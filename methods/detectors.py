"""The two red-cell methods, as label images (one integer id per cell, 0 = background).

  Blob               the collaborator's RBC-specific detector: measure_lesions.rbc_segment,
                     imported unchanged from the frozen handover code.
  Cellpose+HistoPLUS Cellpose draws cell outlines; HistoPLUS labels every pixel with a cell
                     type; each Cellpose cell takes the majority HistoPLUS label of its pixels.
                     A cell is an RBC if that majority label is 'Red blood cell'.

Cellpose can be fed two inputs (both are tested in validation):
  'rgb'        the colour H&E image
  'redexcess'  greyscale 'how much redder than green/blue' = R - (G+B)/2, the image the Blob
               detector works on. (The collaborator found Cellpose sees no RBCs on colour H&E but
               works on this; development_record_addendum.md F.2.)
"""
import json

import numpy as np
from PIL import Image


# ------------------------------------------------------------------ Blob (CPU)
def set_blob_thresholds(red_ref, sat_ref):
    """Calibrated thresholds exactly as check_detector.set_thresholds / measure_lesions.main."""
    from check_detector import set_thresholds
    set_thresholds(red_ref, sat_ref)


def set_blob_thresholds_from_step2(summary_json):
    """Use the thresholds step 2 actually used for this slide (incl. its fallback when the
    false-lumen calibration failed), so Blob here is the study's detector exactly."""
    import measure_lesions as ML
    with open(summary_json) as fh:
        p = json.load(fh)['parameters']
    for k in ('_peak_thr', '_blob_thr', '_mask_red_thr', '_mask_sat_thr'):
        ML.Q[k] = p[k]
    return p.get('calibration')


def blob_labels(rgb, mpp):
    """(label image, peaks (n,2) as (y, x)). One label per counted cell centre."""
    import measure_lesions as ML
    return ML.rbc_segment(rgb, mpp)


# ------------------------------------------------------------------ Cellpose (GPU)
def redexcess_gray(rgb):
    """R - (G+B)/2 as 8-bit greyscale (unchanged from handover scripts/cp_recall.py)."""
    t = rgb.astype(np.float32) / 255.
    r = np.clip(t[..., 0] - (t[..., 1] + t[..., 2]) / 2., 0, None)
    return (r / (r.max() or 1.) * 255).astype(np.uint8)


def build_cellpose(cfg, gpu=True):
    from cellpose import models
    return models.CellposeModel(gpu=gpu, pretrained_model=cfg['cellpose']['pretrained_model'])


def cellpose_labels(model, rgb, mpp, cfg, input_mode):
    c = cfg['cellpose']
    if input_mode == 'redexcess':
        x, kw = redexcess_gray(rgb), {}
    elif input_mode == 'rgb':
        x, kw = rgb, dict(channel_axis=-1)
    else:
        raise ValueError(f'unknown cellpose input {input_mode!r}')
    return model.eval(x, diameter=c['diameter_um'] / mpp, flow_threshold=c['flow_threshold'],
                      cellprob_threshold=c['cellprob_threshold'], min_size=c['min_size'],
                      batch_size=c.get('batch_size', 8), **kw)[0].astype(np.int32)


# ------------------------------------------------------------------ HistoPLUS (GPU)
def build_histoplus(cfg):
    from histoplus.helpers.segmentor.cellvit import CellViTSegmentor
    from histoplus.helpers.hub import histoplus_cellvit_segmentor_40x as M
    h = cfg['histoplus']
    return CellViTSegmentor.from_weights_path(
        weights_path=cfg['paths']['histoplus_weights'], cell_type_mapping=M.class_mapping, mpp=M.mpp,
        segmentor_name=M.name, mean=M.mean, std=M.std, mixed_precision=bool(h['mixed_precision']),
        output_layers=M.output_layers, train_image_size=M.train_image_size, inference_image_size=h['tile'])


def class_names(seg):
    m = seg.class_mapping
    return [m[i] for i in sorted(m)]


def tile_origins(n, tile, overlap):
    """Top-left offsets covering length n, never negative (handover run_histoplus_fields.py)."""
    if n <= tile:
        return [0]
    stride = tile - overlap
    xs = list(range(0, max(n - tile, 0) + 1, stride))
    if xs[-1] != n - tile:
        xs.append(n - tile)
    return xs


def tile_cores(n, tile, overlap):
    """{tile origin: (lo, hi)}: the part of [0, n) each tile owns, so every pixel/cell belongs to
    exactly ONE tile. Boundaries sit at the middle of each pair's overlap. (The earlier rule -
    'origin + overlap/2' - double-claimed a strip when the last tile was shifted back to fit the
    image, so cells there were counted twice.)"""
    o = tile_origins(n, tile, overlap)
    bounds = [0] + [(o[k + 1] + o[k] + tile) // 2 for k in range(len(o) - 1)] + [n]
    return {o[k]: (bounds[k], bounds[k + 1]) for k in range(len(o))}


def histoplus_typemap(seg, rgb, cfg):
    """HistoPLUS per-pixel cell-type map (class ids, uint8) for an RGB array of any size.
    Tiled; each pixel is taken from the tile whose core (tile minus half the overlap on
    interior sides) contains it, so tile edges never decide a label. Reflection padding
    (for arrays smaller than a tile) is cropped off."""
    import torch
    h = cfg['histoplus']
    T, OV, bs = h['tile'], h['overlap'], h['batch_size']
    H, W = rgb.shape[:2]
    padH, padW = max(0, T - H), max(0, T - W)
    if padH or padW:
        rgb = np.pad(rgb, ((0, padH), (0, padW), (0, 0)), mode='reflect')
    Hp, Wp = rgb.shape[:2]
    ycore, xcore = tile_cores(Hp, T, OV), tile_cores(Wp, T, OV)
    out = np.zeros((Hp, Wp), np.uint8)
    origins = [(y0, x0) for y0 in tile_origins(Hp, T, OV) for x0 in tile_origins(Wp, T, OV)]
    for b in range(0, len(origins), bs):
        chunk = origins[b:b + bs]
        batch = torch.stack([seg.transform(Image.fromarray(rgb[y0:y0 + T, x0:x0 + T])) for y0, x0 in chunk])
        with torch.no_grad():
            tp = seg.forward(batch)['tp'].cpu().numpy().astype(np.uint8)
        if tp.shape[1:] != (T, T):
            raise RuntimeError(f'HistoPLUS type map is {tp.shape[1:]}, expected {(T, T)}')
        for (y0, x0), t in zip(chunk, tp):
            cy0, cy1 = ycore[y0]
            cx0, cx1 = xcore[x0]
            out[cy0:cy1, cx0:cx1] = t[cy0 - y0:cy1 - y0, cx0 - x0:cx1 - x0]
    return out[:H, :W]


def _histoplus_instances(seg, rgb, cfg):
    """HistoPLUS's own detections over an RGB array of any size: [(y, x, class, contour_xy)].
    Tiled (tile/overlap from config); a cell is kept once (tile core) and never from reflection
    padding. HistoPLUS centroids and contours are (x, y); contours are shifted to array coords."""
    import torch
    h = cfg['histoplus']
    T, OV = h['tile'], h['overlap']
    H, W = rgb.shape[:2]
    padH, padW = max(0, T - H), max(0, T - W)
    if padH or padW:
        rgb = np.pad(rgb, ((0, padH), (0, padW), (0, 0)), mode='reflect')
    Hp, Wp = rgb.shape[:2]
    ycore, xcore = tile_cores(Hp, T, OV), tile_cores(Wp, T, OV)
    post = seg.get_postprocess_fn()
    out = []
    for y0 in tile_origins(Hp, T, OV):
        for x0 in tile_origins(Wp, T, OV):
            with torch.no_grad():
                o = seg.forward(seg.transform(Image.fromarray(rgb[y0:y0 + T, x0:x0 + T])).unsqueeze(0))
            o = {k: (v.cpu().numpy() if k == 'tp' else v.float().cpu().numpy()) for k, v in o.items()}
            pr = post(o)[0]
            cy0, cy1 = ycore[y0]
            cx0, cx1 = xcore[x0]
            for cen, ct, cont in zip(pr.centroids, pr.cell_types, pr.contours):
                gy, gx = cen[1] + y0, cen[0] + x0
                if cy0 <= gy < cy1 and cx0 <= gx < cx1 and gy < H and gx < W:
                    out.append((gy, gx, str(ct), np.asarray(cont, float) + [x0, y0]))
    return out


def histoplus_rbc_centres(seg, rgb, cfg):
    """HistoPLUS used on its own: centres (y, x) of the cells it detects AND calls RBC."""
    rbc = cfg['histoplus']['rbc_class']
    return np.array([(y, x) for y, x, ct, _ in _histoplus_instances(seg, rgb, cfg) if ct == rbc], float).reshape(-1, 2)


def histoplus_rbc_labels(seg, rgb, cfg):
    """HistoPLUS used on its own, as a label image of the cells it calls RBC (0 = background)."""
    import cv2
    rbc = cfg['histoplus']['rbc_class']
    lab = np.zeros(rgb.shape[:2], np.int32)
    k = 0
    for _, _, ct, cont in _histoplus_instances(seg, rgb, cfg):
        if ct == rbc and len(cont) >= 3:
            k += 1
            cv2.fillPoly(lab, [np.round(cont).astype(np.int32).reshape(-1, 1, 2)], int(k))
    return lab


def cellpose_histoplus_rbc_centres(cp, seg, rgb, mpp, cfg):
    """Cellpose+HistoPLUS: centres (y, x) of Cellpose cells whose majority HistoPLUS label is RBC."""
    from scipy import ndimage as ndi
    lab = cellpose_labels(cp, rgb, mpp, cfg, cfg['cellpose']['lesion_input'])
    cls = classify_cells(lab, histoplus_typemap(seg, rgb, cfg), class_names(seg), cfg['histoplus']['rbc_class'])
    ids = [i for i, (k, _) in cls.items() if k == cfg['histoplus']['rbc_class']]
    if not ids:
        return np.zeros((0, 2))
    return np.array(ndi.center_of_mass(np.ones_like(lab), lab, ids), float).reshape(-1, 2)


def combined_rbc(cp, seg, rgb, mpp, cfg, fill=None):
    """CP(0)+HP + HPfill - the combined method:
      1 Cellpose outlines cells on the red-excess image
      2 each Cellpose cell takes the majority HistoPLUS pixel label; kept if 'Red blood cell'
      3 HistoPLUS's own red-cell outlines are added where Cellpose has no red cell, if the added
        cell is less than fill['max_covered'] covered by Cellpose red cells, at least
        fill['min_area_um2'] in area, and its centre more than fill['min_dist_um'] from every
        Cellpose red-cell centre (rules in config: combined.*)
    Returns (labels, cells): labels = int image of the kept red cells (0 = none);
    cells = list of dict(label, y, x, area_um2, source='cellpose'|'histoplus').
    Third value: the Cellpose cells HistoPLUS did NOT call RBC, as
    dict(cells=[dict(label, y, x, area_um2, hp_class)], labels=<Cellpose label image>)."""
    from scipy import ndimage as ndi
    fill = fill or cfg['combined']
    rbc = cfg['histoplus']['rbc_class']
    names = class_names(seg)
    clab = cellpose_labels(cp, rgb, mpp, cfg, cfg['cellpose']['lesion_input'])
    tmap = histoplus_typemap(seg, rgb, cfg)
    cls = classify_cells(clab, tmap, names, rbc)
    out = np.zeros(clab.shape, np.int32)
    cells, rejected = [], []
    ids = [int(i) for i in np.unique(clab[clab > 0])]
    if ids:
        cen = ndi.center_of_mass(np.ones_like(clab), clab, ids)
        areas = np.bincount(clab.ravel())[ids]
        for lid, (cy, cx), n in zip(ids, cen, areas):
            k = cls.get(lid, ('unlabelled', 0.0))[0]
            d = dict(label=lid, y=float(cy), x=float(cx), area_um2=float(n * mpp * mpp))
            if k == rbc:
                out[clab == lid] = lid
                cells.append(dict(d, source='cellpose'))
            else:
                rejected.append(dict(d, hp_class=k))
    cpc = np.array([[c['y'], c['x']] for c in cells], float).reshape(-1, 2)
    import cv2
    nxt = int(clab.max()) + 1
    H, W = clab.shape
    for _, _, ct, cont in _histoplus_instances(seg, rgb, cfg):
        if ct != rbc or len(cont) < 3:
            continue
        q = np.round(cont).astype(np.int32)
        # rasterise in the cell's own bounding box (a full-image mask per cell is far too slow)
        bx0, by0 = max(0, q[:, 0].min()), max(0, q[:, 1].min())
        bx1, by1 = min(W, q[:, 0].max() + 1), min(H, q[:, 1].max() + 1)
        if bx1 <= bx0 or by1 <= by0:
            continue
        m = np.zeros((by1 - by0, bx1 - bx0), np.uint8)
        cv2.fillPoly(m, [(q - [bx0, by0]).reshape(-1, 1, 2)], 1)
        m = m.astype(bool)
        if not m.any():
            continue
        sub = out[by0:by1, bx0:bx1]
        area = float(m.sum() * mpp * mpp)
        yy, xx = np.nonzero(m)
        cy, cx = yy.mean() + by0, xx.mean() + bx0
        dmin = float(np.min(np.hypot(cpc[:, 0] - cy, cpc[:, 1] - cx)) * mpp) if len(cpc) else 1e9
        if (sub[m] > 0).mean() < fill['max_covered'] and area >= fill['min_area_um2'] and dmin > fill['min_dist_um']:
            sub[m & (sub == 0)] = nxt
            cells.append(dict(label=nxt, y=float(cy), x=float(cx), area_um2=area, source='histoplus'))
            nxt += 1
    return out, cells, dict(cells=rejected, labels=clab)


def classify_cells(labels, typemap, names, rbc_class):
    """Majority HistoPLUS label inside each cell (Background ignored).
    Returns {cell_id: (class_name, rbc_fraction)}; 'unlabelled' if every pixel is Background."""
    K = len(names)
    m = labels > 0
    if not m.any():
        return {}
    n_lab = int(labels.max()) + 1
    counts = np.bincount(labels[m].astype(np.int64) * K + typemap[m], minlength=n_lab * K).reshape(n_lab, K)
    counts[:, 0] = 0                                   # class 0 = Background
    rbc = names.index(rbc_class)
    out = {}
    for i in np.unique(labels[m]):
        c = counts[i]
        tot = c.sum()
        out[int(i)] = ('unlabelled', 0.0) if tot == 0 else (names[int(c.argmax())], float(c[rbc] / tot))
    return out
