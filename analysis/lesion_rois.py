#!/usr/bin/env python3
"""Pictures of whole lesions with all 26 variants, to check counts by eye (no ground truth).

    python lesion_rois.py --slide S07 [--lesions 3 15 16]

One figure per lesion, same layout as the validation figures: legend row, Blob, HistoPLUS alone,
then one row per segmenter (red-excess +HP, red-excess +HP+fill, colour +HP, colour +HP+fill).
Each panel = the lesion's box exactly as the full-slide job processed it (outline extent + 10 um),
pathologist's outline in green, detected red cells outlined (light blue = counted, white = not counted).
The methods run on the SAME regions as scripts/run_all_methods.py (2048 px squares + 64 px margin) with
the same counting rule (centre in the region core and inside the outline), so the counts equal
outputs/all_methods/<slide>/lesions.csv. Lesions are named by their QuPath name; the pipeline's
lesion_id (alphabetical order of QuPath names, as in the original code) is shown too.
Output: outputs/figures/lesion_rois_<slide>/<QuPath name>.jpg
"""
import argparse
import csv
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import cv2  # noqa: E402
import numpy as np  # noqa: E402
import openslide  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from skimage.segmentation import find_boundaries  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_config, out_dir                          # noqa: E402
from core.slides import get_slide, load_geojson, lesions_by_name       # noqa: E402
from methods import detectors as D                                    # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config')
ap.add_argument('--slide', required=True)
ap.add_argument('--lesions', type=int, nargs='*', help='lesion numbers (default: all)')
ap.add_argument('--panel-px', type=int, default=800, help='longest side of each panel')
a = ap.parse_args()
cfg = load_config(a.config)
os.environ.setdefault('MICROSAM_CACHEDIR', os.path.join(os.path.dirname(cfg['paths']['histoplus_weights']), 'microsam'))
import measure_lesions as ML   # noqa: E402
from methods import segmenters as S   # noqa: E402
from methods.figures import COL, INPUT_NAMES, all_label_images, layout, method_title, text_block  # noqa: E402

MPP, O = cfg['slide_mpp'], cfg['paths']['output_dir']
R, M = cfg['dl_lesions']['region_px'], cfg['dl_lesions']['margin_px']
odir = out_dir(cfg, 'figures', f'lesion_rois_{a.slide}')

s = get_slide(cfg, slide_id=a.slide)
sl = openslide.OpenSlide(s['svs'])
W0, H0 = sl.dimensions
lesions = lesions_by_name(load_geojson(s['geojson']))
D.set_blob_thresholds_from_step2(os.path.join(O, 'step2', a.slide, 'summary_measured.json'))
job = {}                                    # counts from the full-slide job, for comparison
jcsv = os.path.join(O, 'all_methods', a.slide, 'lesions.csv')
if os.path.exists(jcsv):
    for r in csv.DictReader(open(jcsv)):
        job[(int(r['lesion_id']), r['method'], r['input'])] = int(r['count'])

segs = S.build(cfg)
hp = D.build_histoplus(cfg)
names = D.class_names(hp)
ITEMS = [('line', COL['lesion'], "Pathologist's lesion outline"),
         ('ring', COL['marked'], 'Red cell counted (centre inside the lesion)'),
         ('ring', COL['fp'], 'Red cell outside the lesion (not counted)')]

for i, (lname, polys) in enumerate(lesions, start=1):
    if a.lesions and i not in a.lesions:
        continue
    # the lesion box and regions exactly as run_all_methods.py
    pts = np.array([q for p in polys for ring in p for q in ring], float)
    X0, Y0 = np.floor(pts.min(0) - 10 / MPP).astype(int)
    X1, Y1 = np.ceil(pts.max(0) + 10 / MPP).astype(int)
    X0, Y0, X1, Y1 = max(0, X0), max(0, Y0), min(W0, X1), min(H0, Y1)
    h, w = Y1 - Y0, X1 - X0
    img = np.asarray(sl.read_region((X0, Y0), 0, (w, h)).convert('RGB'))
    les = ML.rasterise_local(polys, (h, w), lambda x, y: (np.asarray(x) - X0, np.asarray(y) - Y0))
    canv, counted, nxt = {}, {}, {}
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
            core = (slice(ry - oy, cy1 - oy), slice(rx - ox, cx1 - ox))
            dst = (slice(ry - Y0, cy1 - Y0), slice(rx - X0, cx1 - X0))
            for key, (lab, cen) in all_label_images(rgb, MPP, segs, hp, names, cfg, centres=True).items():
                c = canv.setdefault(key, np.zeros((h, w), np.int64))
                off = nxt.get(key, 0)
                sub = lab[core].astype(np.int64)
                c[dst] = np.where(sub > 0, sub + off, 0)
                nxt[key] = off + int(lab.max()) + 1
                ids = set()
                for cy, cx in cen:                       # the job's counting rule
                    iy, ix = int(round(cy)), int(round(cx))
                    if 0 <= iy < ok.shape[0] and 0 <= ix < ok.shape[1] and ok[iy, ix]:
                        ids.add(int(lab[iy, ix]) + off if lab[iy, ix] > 0 else -1)
                        counted.setdefault(key, []).append(int(lab[iy, ix]) + off if lab[iy, ix] > 0 else -1)
    sc = a.panel_px / max(h, w)
    size = (int(round(w * sc)), int(round(h * sc)))
    base = cv2.resize(img, size, interpolation=cv2.INTER_AREA if sc < 1 else cv2.INTER_CUBIC)
    les_b = cv2.resize(les.astype(np.uint8), size, interpolation=cv2.INTER_NEAREST)
    les_edge = cv2.dilate(find_boundaries(les_b, mode='inner').astype(np.uint8), np.ones((3, 3))).astype(bool)
    PW = a.panel_px                                    # every panel column has the same width
    panels, line = {}, []
    for key, c in canv.items():
        ids = counted.get(key, [])
        n = len(ids)
        big = cv2.resize(c.astype(np.int32), size, interpolation=cv2.INTER_NEAREST)
        bnd = find_boundaries(big, mode='inner')
        if sc > 1:
            bnd = cv2.dilate(bnd.astype(np.uint8), np.ones((2, 2))).astype(bool)
        isin = np.isin(big, [k for k in ids if k > 0])
        p = base.copy()
        p[bnd & ~isin] = COL['fp']
        p[bnd & isin] = COL['marked']
        p[les_edge] = COL['lesion']
        p = np.pad(p, ((0, 0), ((PW - p.shape[1]) // 2, PW - p.shape[1] - (PW - p.shape[1]) // 2), (0, 0)),
                   constant_values=255)
        meth, inp = key
        jn = job.get((i, meth, inp))
        head = text_block(PW, [(method_title(meth, inp), 30, True), (INPUT_NAMES[inp], 26, False),
                               (f'red cells counted in the lesion: {n}', 26, False)])
        panels[key] = np.vstack([head, p])
        line.append((meth, inp, n, jn))
    area = les.sum() * MPP * MPP
    grid = layout(panels, list(segs), f'{a.slide}: {lname}',
                  f'pipeline lesion_id {i} · lesion area {area:,.0f} µm² · box {w * MPP:.0f} × {h * MPP:.0f} µm '
                  '· each panel is one method', items=ITEMS)
    Image.fromarray(grid).save(os.path.join(odir, f"{lname.replace(' ', '_') or f'lesion_{i:02d}'}.jpg"), quality=85)
    diff = [(m, inp, n, jn) for m, inp, n, jn in line if jn is not None and n != jn]
    ms = [n for m, inp, n, _ in line if m == 'microsam+HP+fill' and inp == 'redexcess'][0]
    print(f'{lname} (lesion_id {i}): Blob {line[0][2]}, micro-SAM+HP+fill (red-excess) {ms}; '
          f'{len(diff)} of {len(line)} variants differ from lesions.csv' + (f' {diff[:3]}' if diff else ''), flush=True)
print('figures in', odir)
