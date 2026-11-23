#!/usr/bin/env python3
"""Validation on the collaborator's FULLY counted fields: every red cell marked, so false positives count.

    python validate_complete_fields.py [--locate-only]

Input: paths.complete_fields_dir - screenshots of 100 um squares of dense clot (panels of
dense_all_methods.jpg) with every red cell marked green. File name = <real slide name>..._<n marks>.jpg
(private folder); slides are matched through the private slide-code map and named by code from here on.

1. Locate. Field positions are read from outputs/validation/complete_fields/fields.csv when present.
   Otherwise (first run): the square is recomputed with the rule run_all_methods.py used to pick it (densest
   10 um neighbourhood of Blob cells in the lesion with most Blob cells; outputs/dl_lesions/<real id>/
   cells_blob.csv), the screenshot is template-matched (scale + offset) against the clean slide around it
   (stop if score < 0.8), and the result is saved to fields.csv.
   Marks -> slide pixels. Saves <code>_marks.csv and an alignment picture (marks on the clean slide square).
2. Run all 26 variants (Blob, HistoPLUS alone, 6 segmenters x 2 inputs x {+HP, +HP+fill}) on the square
   plus a 200 px margin, with the same code as the full-slide job, at full resolution (0.2524 um/px).
3. Score. A detected red cell belongs to the field if any part of its outline is inside the square.
     found            % of marks inside a detected red-cell outline
     false positives  detected red cells in the field with no mark inside their outline
   Blob: one label per counted peak. Per field and pooled.
Output: outputs/validation/complete_fields.{csv,md}, outputs/validation/complete_fields/<slide>_marks.csv,
        outputs/figures/complete_field_<slide>.jpg (all 26 variants), complete_field_alignment.jpg
"""
import argparse
import csv
import glob
import os
import re
import sys
import warnings

warnings.filterwarnings('ignore')
import cv2  # noqa: E402
import numpy as np  # noqa: E402
import openslide  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from skimage.segmentation import find_boundaries  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from core.config import load_config, out_dir         # noqa: E402
from methods import detectors as D                   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config')
ap.add_argument('--locate-only', action='store_true', help='only step 1 (marks + alignment picture)')
a = ap.parse_args()
cfg = load_config(a.config)
os.environ.setdefault('MICROSAM_CACHEDIR', os.path.join(os.path.dirname(cfg['paths']['histoplus_weights']), 'microsam'))

MPP, RBC, O = cfg['slide_mpp'], cfg['histoplus']['rbc_class'], cfg['paths']['output_dir']
WD, PAD, SEARCH = int(100 / MPP), 200, 150          # square side and margin exactly as run_all_methods.py
VDIR, FDIR = out_dir(cfg, 'validation'), out_dir(cfg, 'figures')
MDIR = out_dir(cfg, 'validation', 'complete_fields')


from methods.figures import COL, INPUT_NAMES, all_label_images as _all_labels, layout, method_title, text_block  # noqa: E402


# ------------------------------------------------------------------ 1. locate
def field_slide(fname):
    """Match a screenshot file name to a slide (its real name, alphanumerics before '_hand'); returns the slide."""
    key = re.sub(r'[^0-9A-Z]', '', fname.upper().split('HAND')[0])
    from core.slides import discover
    norm = lambda x: re.sub(r'[^0-9A-Z]', '', x.upper())
    hits = [s for s in discover(cfg) if norm(s['real_slide_id']).startswith(key) or key.startswith(norm(s['real_slide_id']))]
    hits = sorted(hits, key=lambda s: len(s['real_slide_id']))
    if not hits:
        sys.exit(f'no slide matches {fname}')
    return hits[0]


def green_marks(rgb):
    """Centres (y, x) of the green marks, and a mask of them."""
    r, g, b = (rgb[..., i].astype(int) for i in range(3))
    m = (g > 150) & (g - r > 40) & (g - b > 40)
    lab, n = ndi.label(m)
    sz = ndi.sum(m, lab, range(1, n + 1))
    keep = [i + 1 for i in range(n) if sz[i] >= 0.3 * np.median(sz)]
    return np.array(ndi.center_of_mass(m, lab, keep), float).reshape(-1, 2), m


def dense_corner(real_id):
    """Top-left slide pixel of the dense square, same rule as run_all_methods.py (first run only)."""
    rows = list(csv.DictReader(open(os.path.join(O, 'dl_lesions', real_id, 'cells_blob.csv'))))
    P = np.array([[float(r['x']), float(r['y']), int(r['lesion_id'])] for r in rows])
    top = np.bincount(P[:, 2].astype(int)).argmax()
    Q = P[P[:, 2] == top][:, :2]
    cx, cy = Q[np.argmax([len(x) for x in cKDTree(Q).query_ball_point(Q, 10 / MPP)])]
    return int(cx) - WD // 2, int(cy) - WD // 2


FIELDS_CSV = os.path.join(MDIR, 'fields.csv')
saved = {r['slide_id']: r for r in csv.DictReader(open(FIELDS_CSV))} if os.path.exists(FIELDS_CSV) else {}
fields = []
for path in sorted(glob.glob(os.path.join(cfg['paths']['complete_fields_dir'], '*.jpg'))):
    fname = os.path.basename(path)
    slide = field_slide(fname)
    sid = slide['slide_id']
    n_expected = int(re.findall(r'(\d+)\.jpg$', fname)[0])
    shot = np.asarray(Image.open(path).convert('RGB'))
    marks, mmask = green_marks(shot)
    if len(marks) != n_expected:
        sys.exit(f'{fname}: found {len(marks)} marks, file name says {n_expected}')
    sl = openslide.OpenSlide(slide['svs'])
    if sid in saved:
        r = saved[sid]
        X0, Y0, sx, sy, zoom, score = (int(r['square_x']), int(r['square_y']), int(r['shot_x']), int(r['shot_y']),
                                       float(r['zoom']), float(r['match_score']))
        mk = np.column_stack([sy + marks[:, 0] / zoom, sx + marks[:, 1] / zoom])
        print(f'{sid}: {len(mk)} marks, position from fields.csv', flush=True)
        fields.append(dict(sid=sid, sl=sl, X0=X0, Y0=Y0, sx=sx, sy=sy, marks=mk, score=score, zoom=zoom))
        continue
    X0, Y0 = dense_corner(slide['real_slide_id'])
    reg = np.asarray(sl.read_region((X0 - SEARCH, Y0 - SEARCH), 0, (WD + 2 * SEARCH,) * 2).convert('RGB'))
    clean = cv2.inpaint(shot, cv2.dilate(mmask.astype(np.uint8), np.ones((5, 5)), iterations=2), 5, cv2.INPAINT_TELEA)
    sg = cv2.cvtColor(reg, cv2.COLOR_RGB2GRAY).astype(np.float32)
    tg = cv2.cvtColor(clean, cv2.COLOR_RGB2GRAY).astype(np.float32)
    best = (-1,)
    for z in np.linspace(0.8, 5.0, 421):                 # screenshot px per slide px
        w, h = int(round(shot.shape[1] / z)), int(round(shot.shape[0] / z))
        if w >= sg.shape[1] or h >= sg.shape[0]:
            continue
        res = cv2.matchTemplate(sg, cv2.resize(tg, (w, h), interpolation=cv2.INTER_AREA), cv2.TM_CCOEFF_NORMED)
        _, mx, _, loc = cv2.minMaxLoc(res)
        if mx > best[0]:
            best = (mx, z, loc)
    score, zoom, (lx, ly) = best
    if score < 0.8:
        sys.exit(f'{fname}: screenshot not found on the slide (match score {score:.2f})')
    sx, sy = X0 - SEARCH + lx, Y0 - SEARCH + ly          # slide pixel of the screenshot's top-left
    mk = np.column_stack([sy + marks[:, 0] / zoom, sx + marks[:, 1] / zoom])   # (y, x) slide px
    print(f'{sid}: {len(mk)} marks, match {score:.3f}, zoom {zoom:.2f}, '
          f'screenshot corner vs square {sx - X0:+d},{sy - Y0:+d} px', flush=True)
    fields.append(dict(sid=sid, sl=sl, X0=X0, Y0=Y0, sx=sx, sy=sy, marks=mk, score=score, zoom=zoom))

with open(FIELDS_CSV, 'w', newline='') as fh:        # field positions, so later runs need no re-locating
    w = csv.writer(fh)
    w.writerow(['slide_id', 'n_marks', 'square_x', 'square_y', 'shot_x', 'shot_y', 'zoom', 'match_score'])
    for f in fields:
        w.writerow([f['sid'], len(f['marks']), f['X0'], f['Y0'], f['sx'], f['sy'], round(f['zoom'], 4), round(f['score'], 3)])
for f in fields:
    with open(os.path.join(MDIR, f"{f['sid']}_marks.csv"), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['x', 'y'])
        w.writerows([[round(x, 1), round(y, 1)] for y, x in f['marks']])
fields.sort(key=lambda f: f['sid'])

# alignment picture: marks on the clean slide square (x2), one panel per field
pan = []
for f in fields:
    sq = np.asarray(f['sl'].read_region((f['X0'], f['Y0']), 0, (WD, WD)).convert('RGB')).copy()
    sq = cv2.resize(sq, (2 * WD, 2 * WD), interpolation=cv2.INTER_CUBIC)
    for y, x in f['marks']:
        cv2.circle(sq, (int(2 * (x - f['X0'])), int(2 * (y - f['Y0']))), 6, (0, 255, 0), -1)
    head = text_block(2 * WD, [(f['sid'], 34, True), (f"{len(f['marks'])} marks (green), placed on the clean slide", 26, False)])
    pan.append(np.vstack([head, sq]))
    pan.append(np.full((pan[-1].shape[0], 30, 3), 255, np.uint8))
Image.fromarray(np.hstack(pan[:-1])).save(os.path.join(FDIR, 'complete_field_alignment.jpg'), quality=92)
print('alignment picture:', os.path.join(FDIR, 'complete_field_alignment.jpg'), flush=True)
if a.locate_only:
    sys.exit(0)

# ------------------------------------------------------------------ 2. run all 26 variants
from methods import segmenters as S   # noqa: E402

segs = S.build(cfg)
hp = D.build_histoplus(cfg)
names = D.class_names(hp)


def all_label_images(rgb):
    return _all_labels(rgb, MPP, segs, hp, names, cfg)


def score(lab, marks_local):
    """marks_local (y, x) in lab coordinates; the square is lab[PAD:PAD+WD, PAD:PAD+WD]."""
    inside = np.zeros(lab.shape, bool)
    inside[PAD:PAD + WD, PAD:PAD + WD] = True
    in_field = set(np.unique(lab[inside & (lab > 0)]).tolist())          # any part of the outline inside
    iy = np.clip(np.round(marks_local[:, 0]).astype(int), 0, lab.shape[0] - 1)
    ix = np.clip(np.round(marks_local[:, 1]).astype(int), 0, lab.shape[1] - 1)
    hit = lab[iy, ix]
    marked = set(hit[hit > 0].tolist())
    fp = sorted(in_field - marked)
    return dict(n_marks=len(marks_local), found=int((hit > 0).sum()), detections=len(in_field),
                false_pos=len(fp)), hit > 0, fp


pct = lambda a_, b_: f'{100 * a_ / b_:.0f}%' if b_ else '-'
results = {}
for f in fields:
    x0, y0 = f['X0'] - PAD, f['Y0'] - PAD
    rgb = np.asarray(f['sl'].read_region((x0, y0), 0, (WD + 2 * PAD,) * 2).convert('RGB'))
    D.set_blob_thresholds_from_step2(os.path.join(O, 'step2', f['sid'], 'summary_measured.json'))
    loc = f['marks'] - [y0, x0]
    panels = {}
    for key, lab in all_label_images(rgb).items():
        s, found, fp = score(lab, loc)
        results.setdefault(key, {})[f['sid']] = s
        sub = lab[PAD:PAD + WD, PAD:PAD + WD]
        p = cv2.resize(rgb[PAD:PAD + WD, PAD:PAD + WD], (2 * WD, 2 * WD), interpolation=cv2.INTER_CUBIC)
        big = cv2.resize(sub, (2 * WD, 2 * WD), interpolation=cv2.INTER_NEAREST)
        bnd = cv2.dilate(find_boundaries(big, mode='inner').astype(np.uint8), np.ones((2, 2))).astype(bool)
        isfp = np.isin(big, fp)
        p[bnd & ~isfp] = COL['marked']
        p[bnd & isfp] = COL['fp']
        for (yy, xx), ok in zip(loc - PAD, found):
            cv2.circle(p, (int(2 * xx), int(2 * yy)), 7, COL['found'] if ok else COL['missed'], -1)
            cv2.circle(p, (int(2 * xx), int(2 * yy)), 7, (0, 0, 0), 1)
        meth, inp = key
        head = text_block(2 * WD, [
            (method_title(meth, inp), 34, True), (INPUT_NAMES[inp], 28, False),
            (f"found {s['found']}/{s['n_marks']} ({pct(s['found'], s['n_marks'])})   ·   "
             f"false positives {s['false_pos']}/{s['detections']} ({pct(s['false_pos'], s['detections'])})", 28, False)])
        panels[key] = np.vstack([head, p])
    grid = layout(panels, list(segs), f"{f['sid']}: 100 µm field of dense clot",
                  f"{len(loc)} red cells marked by the collaborator · each panel is one method")
    Image.fromarray(grid).save(os.path.join(FDIR, f"complete_field_{f['sid']}.jpg"), quality=85)
    print(f['sid'], 'scored', flush=True)

# ------------------------------------------------------------------ 3. report
sids = [f['sid'] for f in fields]
pct = lambda a_, b_: f'{100 * a_ / b_:.0f}%' if b_ else '-'
rows = []
L = ['# Fully counted fields: found and false positives', '',
     'Three 100 µm squares of dense clot (one per slide) with every red cell marked.',
     '**Found** = % of marks inside a detected red-cell outline. **False positives** = detected red cells '
     'in the field (any part inside the square) with no mark inside, as number / all detections in the field.',
     '',
     '| Method | Input | ' + ' | '.join(f'{s} ({fields[i]["marks"].shape[0]} marks): found | false pos' for i, s in enumerate(sids))
     + ' | All fields: found | All: false pos |',
     '|---|---|' + '---|---|' * (len(sids) + 1)]
for (meth, inp), d in results.items():
    cells = []
    T = dict(n_marks=0, found=0, detections=0, false_pos=0)
    for s in sids:
        x = d[s]
        for k in T:
            T[k] += x[k]
        cells += [pct(x['found'], x['n_marks']), f"{x['false_pos']}/{x['detections']} ({pct(x['false_pos'], x['detections'])})"]
        rows.append(dict(method=meth, input=inp, field=s, **x))
    cells += [pct(T['found'], T['n_marks']), f"{T['false_pos']}/{T['detections']} ({pct(T['false_pos'], T['detections'])})"]
    rows.append(dict(method=meth, input=inp, field='all', **T))
    L.append(f"| {meth} | {inp} | " + ' | '.join(cells) + ' |')
with open(os.path.join(VDIR, 'complete_fields.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(VDIR, 'complete_fields.md'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('\n'.join(L))
