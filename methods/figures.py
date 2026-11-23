"""Shared figure helpers: fonts, centred titles, legend row, method names, all-variant label images,
and the row layout (legend, Blob, HistoPLUS alone, one row per segmenter x 4 variants)."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

from methods import detectors as D
from methods import segmenters as S

FONT = font_manager.findfont('DejaVu Sans')
FONT_B = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight='bold'))
font = lambda size, bold=False: ImageFont.truetype(FONT_B if bold else FONT, size)
COL = dict(found=(255, 220, 0), missed=(230, 0, 0), marked=(0, 200, 255), fp=(255, 255, 255),
           lesion=(0, 220, 0))
SEG_NAMES = {'cpsam_v2': 'Cellpose v2', 'cpsam': 'Cellpose (cpsam)', 'cpdino': 'Cellpose-DINO',
             'cpdino-vitb': 'Cellpose-DINO ViT-B', 'instanseg': 'InstanSeg', 'microsam': 'micro-SAM'}
INPUT_NAMES = {'redexcess': 'red-excess image', 'rgb': 'colour image', '-': 'colour image'}
VALIDATION_LEGEND = [('dot', COL['found'], 'Marked cell, found'), ('dot', COL['missed'], 'Marked cell, missed'),
                     ('ring', COL['marked'], 'Detected red cell with a mark'),
                     ('ring', COL['fp'], 'False positive (no mark inside)')]


def text_block(width, lines, pad=10):
    """White block with centred text lines [(text, size, bold)]."""
    h = pad + sum(sz + 8 for _, sz, _ in lines) + pad
    im = Image.new('RGB', (width, h), 'white')
    d = ImageDraw.Draw(im)
    y = pad
    for t, sz, b in lines:
        k = sz
        while k > 10 and d.textlength(t, font=font(k, b)) > width - 2 * pad:
            k -= 1                                   # shrink until the line fits the panel
        f = font(k, b)
        d.text(((width - d.textlength(t, font=f)) / 2, y + (sz - k) / 2), t, fill='black', font=f)
        y += sz + 8
    return np.asarray(im)


def legend_row(width, title, subtitle, items=VALIDATION_LEGEND):
    """Whole first row: title + colour key with large swatches. items: [(dot|ring|line, colour, text)]."""
    im = Image.new('RGB', (width, 230), 'white')
    d = ImageDraw.Draw(im)
    for t, y, f in [(title, 14, font(46, True)), (subtitle, 74, font(30))]:
        d.text(((width - d.textlength(t, font=f)) / 2, y), t, fill='black', font=f)
    f = font(32)
    widths = [70 + d.textlength(t, font=f) for _, _, t in items]
    gap = (width - sum(widths)) / (len(items) + 1)
    x, y = gap, 160
    for (kind, c, t), w in zip(items, widths):
        if kind == 'dot':
            d.ellipse([x + 8, y - 16, x + 40, y + 16], fill=c, outline='black')
        elif kind == 'line':
            d.rectangle([x, y - 22, x + 48, y + 22], fill=(200, 120, 140))
            d.line([x + 2, y, x + 46, y], fill=c, width=7)
        else:
            d.rectangle([x, y - 22, x + 48, y + 22], fill=(200, 120, 140))
            d.ellipse([x + 4, y - 18, x + 44, y + 18], outline=c, width=5)
        d.text((x + 64, y - 20), t, fill='black', font=f)
        x += w + gap
    return np.asarray(im)


def method_title(meth, inp):
    if meth in ('Blob', 'HistoPLUS alone'):
        return meth + (' (original detector, no AI)' if meth == 'Blob' else '')
    seg, ver = meth.split('+', 1)
    return f"{SEG_NAMES.get(seg, seg)} + HistoPLUS{' + fill' if ver.endswith('fill') else ''}"


def all_label_images(rgb, mpp, segs, hp, names, cfg, centres=False):
    """{(method, input): label image of red cells} for all 26 variants; one label per detected red cell.
    Same steps as scripts/run_all_methods.py. With centres=True: {key: (label image, centres (n,2) y,x)},
    where the centres are the points the full-slide job counts (Blob peaks; otherwise cell centroids)."""
    rbc = cfg['histoplus']['rbc_class']
    out, cen = {}, {}
    blab, peaks = D.blob_labels(rgb, mpp)
    out[('Blob', '-')] = blab.astype(np.int32)
    cen[('Blob', '-')] = np.asarray(peaks, float).reshape(-1, 2)
    tmap = D.histoplus_typemap(hp, rgb, cfg)
    hpi = D._histoplus_instances(hp, rgb, cfg)
    hpm = S.hp_rbc_masks(hpi, rgb.shape[:2], rbc)
    hl = np.zeros(rgb.shape[:2], np.int32)
    for k, (m, y0, x0) in enumerate(hpm, start=1):
        hl[y0:y0 + m.shape[0], x0:x0 + m.shape[1]][m] = k
    out[('HistoPLUS alone', '-')] = hl
    cen[('HistoPLUS alone', '-')] = np.array([(yy.mean() + y0, xx.mean() + x0) for m, y0, x0 in hpm
                                              for yy, xx in [np.nonzero(m)]], float).reshape(-1, 2)
    for name, fn in segs.items():
        for inp in S.INPUTS:
            cells, labs = S.red_cells(fn(rgb, mpp, inp), tmap, names, hpm, cfg, labels=True)
            for ver, lab in labs.items():
                out[(name + ver, inp)] = lab
                cen[(name + ver, inp)] = np.array([(cy, cx) for cy, cx, _ in cells[ver]], float).reshape(-1, 2)
    return {k: (v, cen[k]) for k, v in out.items()} if centres else out


def layout(panels, seg_names, title, subtitle, items=VALIDATION_LEGEND, gap=30):
    """panels {(method, input): image}; rows = legend, Blob, HistoPLUS alone, then per segmenter
    (red-excess +HP, red-excess +HP+fill, colour +HP, colour +HP+fill). Panels padded to equal size."""
    H = max(p.shape[0] for p in panels.values())
    W = max(p.shape[1] for p in panels.values())
    pad = lambda p: np.pad(p, ((0, H - p.shape[0]), (0, W - p.shape[1]), (0, 0)), constant_values=255)
    blank = np.full((H, W, 3), 255, np.uint8)
    rows = [[panels[('Blob', '-')]], [panels[('HistoPLUS alone', '-')]]]
    for seg in seg_names:
        rows.append([panels[(seg + v, i)] for i in S.INPUTS for v in ('+HP', '+HP+fill')])
    ncol = max(len(r) for r in rows)
    width = ncol * W + (ncol - 1) * gap
    out = [legend_row(width, title, subtitle, items)]
    for r in rows:
        r = [pad(p) for p in r] + [blank] * (ncol - len(r))
        line = []
        for k, p in enumerate(r):
            line.append(p)
            if k < ncol - 1:
                line.append(np.full((H, gap, 3), 255, np.uint8))
        out += [np.full((gap, width, 3), 255, np.uint8), np.hstack(line)]
    return np.vstack(out)
