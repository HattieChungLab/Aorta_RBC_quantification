"""General-purpose cell segmenters, each usable on two inputs, plus the HistoPLUS red-cell steps.

Segmenters (none has a red-cell class of its own):
  cpsam_v2, cpsam, cpdino, cpdino-vitb   Cellpose models
  instanseg   InstanSeg: 'fluorescence_nuclei_and_cells' (cells channel) on red-excess,
              'brightfield_nuclei' on the colour image
  microsam    micro-SAM 'vit_b_lm' automatic instance segmentation, tiled in 448 px tiles
Inputs:
  redexcess   greyscale R - (G+B)/2 (red cells bright)        rgb   the plain colour H&E
Red-cell steps (shared by every segmenter):
  +HP         keep an outline if its majority HistoPLUS pixel label is 'Red blood cell'
  +HP+fill    also add HistoPLUS's own red-cell outlines that touch nothing already kept
"""
import numpy as np

from methods import detectors as D

ALL = ['cpsam_v2', 'cpsam', 'cpdino', 'cpdino-vitb', 'instanseg', 'microsam']
INPUTS = ['redexcess', 'rgb']


def build(cfg, names=ALL):
    """{segmenter: fn(rgb, mpp, input) -> int label image}"""
    c0 = cfg['cellpose']
    segs = {}
    for m in [n for n in names if n.startswith('cp')]:
        from cellpose import models
        mdl = models.CellposeModel(gpu=True, pretrained_model=m)

        def run_cp(rgb, mpp, inp, mdl=mdl):
            img, kw = (D.redexcess_gray(rgb), {}) if inp == 'redexcess' else (rgb, dict(channel_axis=-1))
            return mdl.eval(img, diameter=c0['diameter_um'] / mpp, flow_threshold=c0['flow_threshold'],
                            cellprob_threshold=c0['cellprob_threshold'], min_size=c0['min_size'],
                            batch_size=c0['batch_size'], **kw)[0].astype(np.int32)
        segs[m] = run_cp
    if 'instanseg' in names:
        from instanseg import InstanSeg
        ins = {'redexcess': InstanSeg('fluorescence_nuclei_and_cells', verbosity=0),
               'rgb': InstanSeg('brightfield_nuclei', verbosity=0)}

        def run_ins(rgb, mpp, inp):
            img = D.redexcess_gray(rgb)[..., None] if inp == 'redexcess' else rgb
            out = ins[inp].eval_medium_image(img, pixel_size=mpp, tile_size=512)
            lab = out[0] if isinstance(out, (tuple, list)) else out
            lab = np.squeeze(lab.detach().cpu().numpy() if hasattr(lab, 'detach') else np.asarray(lab))
            if lab.ndim == 3:                       # (nuclei, cells) -> cells
                lab = lab[-1]
            return lab.astype(np.int32)
        segs['instanseg'] = run_ins
    if 'microsam' in names:
        from micro_sam.automatic_segmentation import get_predictor_and_segmenter, automatic_instance_segmentation
        pred, sg = get_predictor_and_segmenter(model_type='vit_b_lm', device='cuda', is_tiled=True)

        def run_ms(rgb, mpp, inp):
            img = D.redexcess_gray(rgb) if inp == 'redexcess' else rgb
            return np.asarray(automatic_instance_segmentation(
                predictor=pred, segmenter=sg, input_path=img, ndim=2, tile_shape=(448, 448), halo=(64, 64),
                verbose=False)).astype(np.int32)
        segs['microsam'] = run_ms
    return segs


def rasterise(cont, shape):
    """HistoPLUS contour (x, y) -> (bool mask in its bounding box, y0, x0) or None."""
    import cv2
    q = np.round(cont).astype(np.int32)
    H, W = shape
    x0, y0 = max(0, q[:, 0].min()), max(0, q[:, 1].min())
    x1, y1 = min(W, q[:, 0].max() + 1), min(H, q[:, 1].max() + 1)
    if x1 <= x0 or y1 <= y0:
        return None
    m = np.zeros((y1 - y0, x1 - x0), np.uint8)
    cv2.fillPoly(m, [(q - [x0, y0]).reshape(-1, 1, 2)], 1)
    m = m.astype(bool)
    return (m, y0, x0) if m.any() else None


def hp_rbc_masks(hpi, shape, rbc):
    """HistoPLUS's own red-cell outlines as [(mask, y0, x0)]."""
    out = []
    for _, _, ct, cont in hpi:
        if ct == rbc and len(cont) >= 3:
            r = rasterise(cont, shape)
            if r is not None:
                out.append(r)
    return out


def cells_of(lab, ids):
    """[(label, cy, cx, n_pixels)] for the given labels."""
    if not len(ids):
        return []
    from scipy import ndimage as ndi
    ids = list(ids)
    cen = ndi.center_of_mass(np.ones_like(lab), lab, ids)
    n = np.bincount(lab.ravel())[ids]
    return [(int(i), float(c[0]), float(c[1]), int(k)) for i, c, k in zip(ids, cen, n)]


def red_cells(lab, tmap, names, hp_masks, cfg, labels=False):
    """From one segmenter's label image -> {'+HP': cells, '+HP+fill': cells};
    a cell = (cy, cx, n_pixels). With labels=True also returns {'+HP': label image, '+HP+fill': label image}
    (fill cells get new label ids above the segmenter's)."""
    rbc = cfg['histoplus']['rbc_class']
    cls = D.classify_cells(lab, tmap, names, rbc)
    keep = [l for l, (k, _) in cls.items() if k == rbc]
    hp = [(cy, cx, n) for _, cy, cx, n in cells_of(lab, keep)]
    occ = np.isin(lab, keep)
    fill = list(hp)
    if labels:
        lab_hp = np.where(occ, lab, 0).astype(np.int32)
        lab_fill, nxt = lab_hp.copy(), int(lab_hp.max()) + 1
    rule, mpp2 = cfg['combined'], cfg['slide_mpp'] ** 2
    for m, y0, x0 in hp_masks:
        sub = occ[y0:y0 + m.shape[0], x0:x0 + m.shape[1]]
        # same rule as the combined method: add only if (almost) nothing there yet and big enough
        if sub[m].mean() < rule['max_covered'] and m.sum() * mpp2 >= rule['min_area_um2']:
            sub[m] = True
            yy, xx = np.nonzero(m)
            fill.append((yy.mean() + y0, xx.mean() + x0, int(m.sum())))
            if labels:
                lab_fill[y0:y0 + m.shape[0], x0:x0 + m.shape[1]][m] = nxt
                nxt += 1
    if labels:
        return {'+HP': hp, '+HP+fill': fill}, {'+HP': lab_hp, '+HP+fill': lab_fill}
    return {'+HP': hp, '+HP+fill': fill}
