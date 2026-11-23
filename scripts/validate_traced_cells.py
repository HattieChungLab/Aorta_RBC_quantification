#!/usr/bin/env python3
"""Validation on the pathologist's traced red cells (a sample of individually outlined cells).

    python validate_traced_cells.py [--limit N]

Slides and Blob calibration values: config `validation` (traced cells = 'RBC' outlines in each slide's
QuPath GeoJSON). For every traced cell, a validation_window_px square centred on it is read at full
resolution and every variant is run on it, with the same code as the full-slide job
(methods.figures.all_label_images): Blob, HistoPLUS alone, 6 segmenters x 2 inputs x {+HP, +HP+fill}.

Scores per slide:
  found          a detected red cell overlaps the traced outline (any overlap)
  counted once   exactly one counted cell centre lies inside the traced outline
                 (0 = missed, 2 or more = split)
Output: outputs/validation/traced_cells.{csv,md}
"""
import argparse
import csv
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np  # noqa: E402
import openslide  # noqa: E402
from skimage.draw import polygon as fill_polygon  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_config, out_dir                    # noqa: E402
from core.slides import get_slide, load_geojson, rbc_cells      # noqa: E402
from methods import detectors as D                              # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config')
ap.add_argument('--limit', type=int, help='first N traced cells per slide (smoke test)')
a = ap.parse_args()
cfg = load_config(a.config)
os.environ.setdefault('MICROSAM_CACHEDIR', os.path.join(os.path.dirname(cfg['paths']['histoplus_weights']), 'microsam'))
from methods import segmenters as S                       # noqa: E402
from methods.figures import all_label_images              # noqa: E402

MPP, WIN = cfg['slide_mpp'], cfg['validation_window_px']
segs = S.build(cfg)
hp = D.build_histoplus(cfg)
names = D.class_names(hp)

res, t0 = {}, time.time()        # (method, input) -> {slide: dict(found, once, missed, split, n)}
for v in cfg['validation']:
    D.set_blob_thresholds(v['red_ref'], v['sat_ref'])
    s = get_slide(cfg, slide_id=v['slide_id'])
    sl = openslide.OpenSlide(s['svs'])
    cells = rbc_cells(load_geojson(s['geojson']), MPP)[:a.limit]
    acc = {}
    for i, c in enumerate(cells):
        x0, y0 = int(round(c['cx'])) - WIN // 2, int(round(c['cy'])) - WIN // 2
        rgb = np.asarray(sl.read_region((x0, y0), 0, (WIN, WIN)).convert('RGB'))
        r = np.asarray(c['ring'], float)
        rr, cc = fill_polygon(r[:, 1] - y0 - 0.5, r[:, 0] - x0 - 0.5, (WIN, WIN))
        hand = np.zeros((WIN, WIN), bool)
        hand[rr, cc] = True
        for key, (lab, cen) in all_label_images(rgb, MPP, segs, hp, names, cfg, centres=True).items():
            iy, ix = np.round(cen[:, 0]).astype(int), np.round(cen[:, 1]).astype(int)
            ok = (iy >= 0) & (iy < WIN) & (ix >= 0) & (ix < WIN)
            acc.setdefault(key, []).append((bool((lab[hand] > 0).any()), int(hand[iy[ok], ix[ok]].sum())))
        if i % 10 == 0:
            print(f"{v['slide_id']} {i + 1}/{len(cells)}  [{time.time() - t0:.0f} s]", flush=True)
    for key, lst in acc.items():
        ov = np.array([x for x, _ in lst]); n = np.array([y for _, y in lst])
        res.setdefault(key, {})[v['slide_id']] = dict(n=len(lst), found=ov.mean(), once=(n == 1).mean(),
                                                      missed=(n == 0).mean(), split=(n >= 2).mean())

sids = [v['slide_id'] for v in cfg['validation']]
batch = {v['slide_id']: v['batch'] for v in cfg['validation']}
pc = lambda x: f'{100 * x:.0f}%'
rows, L = [], ['# Traced cells: found and counted once', '',
               'Found = a detected red cell overlaps the traced outline. Counted once = exactly one counted cell centre '
               'inside the outline (missed = none, split = two or more).', '',
               '| Method | Input | ' + ' | '.join(f'{s} ({batch[s]}, {res[next(iter(res))][s]["n"]} cells): found | '
                                                 f'{s}: counted once (missed / split)' for s in sids) + ' |',
               '|---|---|' + '---|---|' * len(sids)]
for (meth, inp), d in res.items():
    cols = []
    for s in sids:
        x = d[s]
        cols += [pc(x['found']), f"{pc(x['once'])} ({pc(x['missed'])} / {pc(x['split'])})"]
        rows.append(dict(method=meth, input=inp, slide=s, **{k: round(float(val), 4) for k, val in x.items()}))
    L.append(f'| {meth} | {inp} | ' + ' | '.join(cols) + ' |')
vdir = out_dir(cfg, 'validation')
with open(os.path.join(vdir, 'traced_cells.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(vdir, 'traced_cells.md'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('\n'.join(L))
