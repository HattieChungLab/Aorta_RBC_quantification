#!/usr/bin/env python3
"""Full results of the all-methods job: per-slide red-cell COUNT and AREA for every method variant.

    python report_all_methods.py [--config CONFIG.yaml] [--source all_methods]

Reads outputs/all_methods/<slide>/lesions.csv and outputs/step2/<slide>/summary_measured.json.
Writes outputs/all_methods_report/:
  counts_per_slide.csv   one row per slide, one column per method variant (red-cell count)
  area_per_slide.csv     same, red-cell area in mm2
  slide_info.csv         per slide: stain (red_ref), lesions, lesion area, media area excl. false lumen,
                         QC flags (calibration failed, paler than validated)
  report.md              the above as readable tables + cohort totals per variant
Variant names: <segmenter>+HP or +HP+fill, followed by the input in brackets (red-excess or colour);
Blob and HistoPLUS alone have no input.
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_config, out_dir   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config')
ap.add_argument('--source', default='all_methods')
a = ap.parse_args()
cfg = load_config(a.config)
O = cfg['paths']['output_dir']
rdir = out_dir(cfg, a.source + '_report')
LBL = {'redexcess': 'red-excess', 'rgb': 'colour', '-': ''}

counts, areas, info = OrderedDict(), OrderedDict(), []
variants = []
for d in sorted(glob.glob(os.path.join(O, a.source, '*', 'lesions.csv'))):
    sid = os.path.basename(os.path.dirname(d))
    c, ar = OrderedDict(), OrderedDict()
    lids = set()
    for r in csv.DictReader(open(d)):
        v = r['method'] + (f" ({LBL[r['input']]})" if r['input'] != '-' else '')
        if v not in variants:
            variants.append(v)
        c[v] = c.get(v, 0) + int(r['count'])
        ar[v] = ar.get(v, 0.0) + float(r['area_um2'])
        lids.add(r['lesion_id'])
    counts[sid], areas[sid] = c, ar
    sj = os.path.join(O, 'step2', sid, 'summary_measured.json')
    s2 = json.load(open(sj)) if os.path.exists(sj) else {}
    cal = (s2.get('parameters') or {}).get('calibration')
    flags = []
    if not cal:
        flags.append('Blob calibration failed')
    elif cal['red_ref'] < cfg['qc']['pale_red_ref']:
        flags.append('paler than validated')
    info.append(dict(slide_id=sid, red_ref=cal['red_ref'] if cal else '', lesions_done=f"{len(lids)}/{s2.get('n_lesions', '?')}",
                     lesion_area_mm2=s2.get('total_lesion_area_mm2', ''), media_excl_FL_mm2=s2.get('media_area_excl_FL_mm2', ''),
                     flags='; '.join(flags)))
if not counts:
    raise SystemExit(f'no results under {a.source}/')


def write(name, table, fmt):
    with open(os.path.join(rdir, name), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['slide_id'] + variants)
        for sid, row in table.items():
            w.writerow([sid] + [fmt(row.get(v, 0)) for v in variants])


write('counts_per_slide.csv', counts, lambda x: int(x))
write('area_per_slide.csv', areas, lambda x: round(x / 1e6, 5))
with open(os.path.join(rdir, 'slide_info.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, list(info[0].keys())); w.writeheader(); w.writerows(info)

L = ['# All methods, all slides: red-cell count and area inside the lesions', '',
     '"+HP" = outlines HistoPLUS labels red blood cell; "+HP+fill" = plus HistoPLUS\'s own red-cell outlines where '
     'nothing else was found. Input in brackets: the image given to the segmenter (red-excess or colour).', '',
     '## Cohort totals per method (all slides)', '', '| method | red cells | red-cell area (mm²) |', '|---|---|---|']
for v in variants:
    L.append(f"| {v} | {sum(t.get(v, 0) for t in counts.values()):,} | {sum(t.get(v, 0) for t in areas.values()) / 1e6:.3f} |")
L += ['', '## Slides', '', '| slide | red_ref | lesions | lesion area mm² | media excl. FL mm² | flags |', '|---|---|---|---|---|---|']
L += [f"| {r['slide_id']} | {r['red_ref'] or 'failed'} | {r['lesions_done']} | {r['lesion_area_mm2']} | {r['media_excl_FL_mm2']} | {r['flags']} |"
      for r in info]
for title, table, fmt in (('Red-cell COUNT per slide', counts, lambda x: f'{int(x):,}'),
                          ('Red-cell AREA per slide (mm²)', areas, lambda x: f'{x / 1e6:.4f}')):
    L += ['', f'## {title}', '', '| method | ' + ' | '.join(table) + ' |', '|---|' + '---|' * len(table)]
    for v in variants:
        L.append(f'| {v} | ' + ' | '.join(fmt(table[s].get(v, 0)) for s in table) + ' |')
L += ['', 'Pictures: outputs/all_methods/<slide>/dense_all_methods.jpg (every method on the densest 100 µm of each slide).']
with open(os.path.join(rdir, 'report.md'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print(f'{len(counts)} slides, {len(variants)} method variants -> {rdir}')
