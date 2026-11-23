#!/usr/bin/env python3
"""ONE job: study numbers (step 2) + all 26 method variants inside every lesion of every slide.

    python run_all_methods.py [--slides ID ...] [--lesion-ids 1,2] [--workers 8]

--slides ID ...  only these slides; their step 2 is always recomputed (use after annotations change;
                 move the slide's old outputs/all_methods/<slide>/ away first, or finished lesions are skipped)

Phase 1  step 2 (Blob study numbers) recomputed for every slide, in parallel on the job's CPUs
         (measure_lesions.py, unchanged). Skipped on a resubmit once outputs/all_methods/step2_done.txt exists.
Phase 2  per slide, per lesion (square regions of region_px + margin, as run_dl_lesions.py), all variants:
           Blob                        collaborator's detector, thresholds from this slide's step 2
           HistoPLUS alone             HistoPLUS's own red-cell outlines
           <segmenter>+HP              6 segmenters x 2 inputs (red-excess, colour); outlines kept if
           <segmenter>+HP+fill         HistoPLUS labels them red cell; +fill adds HistoPLUS red cells
                                       that touch nothing already kept
         A cell counts if its centre is in the region core and inside the lesion outline.
         Per slide one picture: all variants on the densest 100 um of its biggest lesion.
Output   outputs/all_methods/<slide>/lesions.csv   long format: lesion_id, method, input, count, area_um2
         outputs/all_methods/<slide>/dense_all_methods.jpg
Resumable: lesions already in a slide's lesions.csv are skipped. --lesion-ids writes to all_methods_smoke/.
Then: analysis/report_all_methods.py
"""
import argparse
import csv
import os
import subprocess
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

warnings.filterwarnings('ignore')
import cv2  # noqa: E402
import numpy as np  # noqa: E402
import openslide  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from core.config import load_config, out_dir                      # noqa: E402
from core.slides import discover, load_geojson, lesions_by_name    # noqa: E402
from methods import detectors as D                                # noqa: E402
from methods import segmenters as S                               # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config')
ap.add_argument('--slides', nargs='*')
ap.add_argument('--lesion-ids', help='comma-separated lesion ids only (smoke test; separate output folder)')
ap.add_argument('--workers', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
ap.add_argument('--skip-step2', action='store_true')
a = ap.parse_args()
cfg = load_config(a.config)
os.environ.setdefault('MICROSAM_CACHEDIR', os.path.join(os.path.dirname(cfg['paths']['histoplus_weights']), 'microsam'))
import measure_lesions as ML   # noqa: E402

MPP, R, M = cfg['slide_mpp'], cfg['dl_lesions']['region_px'], cfg['dl_lesions']['margin_px']
RBC, CAP = cfg['histoplus']['rbc_class'], ML.Q['cell_area_cap_um2']
O = cfg['paths']['output_dir']
slides = [s for s in discover(cfg) if not a.slides or s['slide_id'] in a.slides]
root = out_dir(cfg, 'all_methods' + ('_smoke' if a.lesion_ids else ''))
t0 = time.time()
say = lambda m: print(f'[{(time.time() - t0) / 60:6.1f} min] {m}', flush=True)

# ------------------------------------------------------------------ phase 1: step 2
done_flag = os.path.join(root, 'step2_done.txt')
# step 2 runs on: every slide the first time; afterwards only slides named with --slides (always redone,
# e.g. after their annotations changed) plus any slide that has no step-2 output yet
s2 = [] if a.skip_step2 else slides if not os.path.exists(done_flag) or a.slides else \
    [s for s in slides if not os.path.exists(os.path.join(O, 'step2', s['slide_id'], 'summary_measured.json'))]
if s2:
    say(f'step 2 on {len(s2)} slides, {a.workers} at a time')

    def step2(s):
        r = subprocess.run([sys.executable, '-W', 'ignore', os.path.join(HERE, 'scripts', 'run_measure.py'),
                            '--slide-id', s['slide_id'], '--redo'] + (['--config', a.config] if a.config else []),
                           capture_output=True, text=True)
        return s['slide_id'], r.returncode, r.stdout[-400:] + r.stderr[-400:]
    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(step2, s2))
    bad = [(sid, msg) for sid, rc, msg in res if rc]
    for sid, msg in bad:
        say(f'step 2 FAILED {sid}: {msg}')
    if bad:
        raise SystemExit('step 2 failed on some slides - see above')
    if not a.slides:
        open(done_flag, 'w').write(time.strftime('%Y-%m-%d %H:%M\n'))
    say('step 2 done')

# ------------------------------------------------------------------ phase 2: all methods
segs = S.build(cfg)
hp = D.build_histoplus(cfg)
names = D.class_names(hp)
say(f'models ready: {list(segs)}')


def all_variants(rgb, mpp):
    """{(method, input): [(cy, cx, n_pixels or area_um2-for-blob)]} plus blob labels."""
    out = {}
    blab, peaks = D.blob_labels(rgb, mpp)
    bcount = np.bincount(blab.ravel())
    out[('Blob', '-')] = [(float(py), float(px), int(bcount[blab[py, px]]) if blab[py, px] > 0 else 0) for py, px in peaks]
    tmap = D.histoplus_typemap(hp, rgb, cfg)
    hpi = D._histoplus_instances(hp, rgb, cfg)
    hpm = S.hp_rbc_masks(hpi, rgb.shape[:2], RBC)
    out[('HistoPLUS alone', '-')] = [(yy.mean() + y0, xx.mean() + x0, int(m.sum()))
                                     for m, y0, x0 in hpm for yy, xx in [np.nonzero(m)]]
    for name, fn in segs.items():
        for inp in S.INPUTS:
            for ver, cells in S.red_cells(fn(rgb, mpp, inp), tmap, names, hpm, cfg).items():
                out[(name + ver, inp)] = cells
    return out, blab


FIELDS = ['slide_id', 'lesion_id', 'qupath_name', 'method', 'input', 'count', 'area_um2']
for si, s in enumerate(slides, start=1):
    sid = s['slide_id']
    odir = out_dir(cfg, os.path.relpath(root, O), sid)
    les_csv = os.path.join(odir, 'lesions.csv')
    done = {int(r['lesion_id']) for r in csv.DictReader(open(les_csv))} if os.path.exists(les_csv) else set()
    D.set_blob_thresholds_from_step2(os.path.join(O, 'step2', sid, 'summary_measured.json'))
    sl = openslide.OpenSlide(s['svs'])
    W0, H0 = sl.dimensions
    lesions = lesions_by_name(load_geojson(s['geojson']))
    wanted = {int(x) for x in a.lesion_ids.split(',')} if a.lesion_ids else None
    say(f'slide {si}/{len(slides)} {sid}: {len(lesions)} lesions, {len(done)} already done')
    blob_pts = []
    for i, (lname, polys) in enumerate(lesions, start=1):
        if i in done or (wanted is not None and i not in wanted):
            continue
        tl = time.time()
        pts = np.array([q for p in polys for ring in p for q in ring], float)
        X0, Y0 = np.floor(pts.min(0) - 10 / MPP).astype(int)
        X1, Y1 = np.ceil(pts.max(0) + 10 / MPP).astype(int)
        X0, Y0, X1, Y1 = max(0, X0), max(0, Y0), min(W0, X1), min(H0, Y1)
        tot = {}
        for ry in range(Y0, Y1, R):
            for rx in range(X0, X1, R):
                cx1, cy1 = min(rx + R, X1), min(ry + R, Y1)
                ox, oy = max(0, rx - M), max(0, ry - M)
                ex, ey = min(W0, cx1 + M), min(H0, cy1 + M)
                rgb = np.asarray(sl.read_region((ox, oy), 0, (ex - ox, ey - oy)).convert('RGB'))
                inside = ML.rasterise_local(polys, rgb.shape[:2], lambda x, y: (np.asarray(x) - ox, np.asarray(y) - oy))
                if not inside.any():
                    continue
                ok = np.zeros(inside.shape, bool)
                ok[ry - oy:cy1 - oy, rx - ox:cx1 - ox] = True
                ok &= inside
                res, blab = all_variants(rgb, MPP)
                bin_ = np.bincount((blab * inside).ravel())        # Blob area inside the outline, as step 2
                for key, cells in res.items():
                    c, ar = tot.get(key, (0, 0.0))
                    for cy, cx, n in cells:
                        iy, ix = int(round(cy)), int(round(cx))
                        if not (0 <= iy < ok.shape[0] and 0 <= ix < ok.shape[1] and ok[iy, ix]):
                            continue
                        c += 1
                        if key[0] == 'Blob':
                            lid = blab[iy, ix]
                            ar += min(bin_[lid] * MPP * MPP, CAP) if lid > 0 else 0.0
                            blob_pts.append((cx + ox, cy + oy, i))
                        else:
                            ar += n * MPP * MPP
                    tot[key] = (c, ar)
        rows = [dict(slide_id=sid, lesion_id=i, qupath_name=lname, method=k[0], input=k[1], count=c, area_um2=round(ar, 1))
                for k, (c, ar) in tot.items()]
        new = not os.path.exists(les_csv)
        with open(les_csv, 'a', newline='') as fh:           # written after EVERY lesion (resumable)
            w = csv.DictWriter(fh, FIELDS)
            if new:
                w.writeheader()
            w.writerows(rows)
        b = tot.get(('Blob', '-'), (0, 0))[0]
        cp = tot.get(('cpsam_v2+HP+fill', 'redexcess'), (0, 0))[0]
        say(f'  {sid} lesion {i}/{len(lesions)}: Blob {b}, cpsam_v2+HP+fill(red-excess) {cp}  [{time.time() - tl:.0f} s]')

    # ---- one picture per slide: all variants on the densest 100 um of the lesion with most Blob cells
    if blob_pts and not os.path.exists(os.path.join(odir, 'dense_all_methods.jpg')):
        P = np.array(blob_pts)
        top = np.bincount(P[:, 2].astype(int)).argmax()
        Q = P[P[:, 2] == top][:, :2]
        cx, cy = Q[np.argmax([len(x) for x in cKDTree(Q).query_ball_point(Q, 10 / MPP)])]
        Wd, PAD = int(100 / MPP), 200
        big = np.asarray(sl.read_region((int(cx) - Wd // 2 - PAD, int(cy) - Wd // 2 - PAD), 0,
                                        (Wd + 2 * PAD, Wd + 2 * PAD)).convert('RGB'))
        res, _ = all_variants(big, MPP)
        panels = []
        for (meth, inp), cells in res.items():
            p = big[PAD:PAD + Wd, PAD:PAD + Wd].copy()
            n = 0
            for y, x, _ in cells:
                if PAD <= y < PAD + Wd and PAD <= x < PAD + Wd:
                    n += 1
                    cv2.circle(p, (int(x) - PAD, int(y) - PAD), 3, (0, 255, 255) if meth != 'Blob' else (0, 200, 0), -1)
            bar = np.full((22, Wd, 3), 255, np.uint8)
            cv2.putText(bar, f"{meth}{'' if inp == '-' else ' (' + inp + ')'}: {n}", (4, 16), cv2.FONT_HERSHEY_SIMPLEX,
                        0.42, (0, 0, 0), 1, cv2.LINE_AA)
            panels.append(np.vstack([bar, p]))
        cols = 6
        while len(panels) % cols:
            panels.append(np.full_like(panels[0], 255))
        grid = np.vstack([np.hstack(panels[k:k + cols]) for k in range(0, len(panels), cols)])
        Image.fromarray(grid).save(os.path.join(odir, 'dense_all_methods.jpg'), quality=88)
say('ALL DONE')
