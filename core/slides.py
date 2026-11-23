"""Slide discovery and QuPath GeoJSON reading.

Each slide lives in its own folder under data_dir with exactly one .svs and one .geojson.
macOS resource-fork files ('._*') are ignored. Folder names may contain spaces (one has a
trailing space), so paths are always handled as whole strings, never split.

Slides are identified by neutral codes (S01, S02, ...) everywhere in outputs. The private,
gitignored file paths.slide_codes (columns code, real_slide_id) maps each code to the slide's real
name; real names never appear in output folder or file names.
"""
import csv
import glob
import json
import os
import re

import numpy as np


def slide_id_of(svs_path):
    """Same safe name run_all.py used for its Results/ folders (the real name; replaced by its code in discover())."""
    base = os.path.splitext(os.path.basename(svs_path))[0]
    return re.sub(r'[^A-Za-z0-9]+', '_', base).strip('_')


def _visible(paths):
    return sorted(p for p in paths if not os.path.basename(p).startswith('._'))


def slide_codes(cfg):
    """{real_slide_id: code} from the private map (paths.slide_codes)."""
    path = cfg['paths'].get('slide_codes')
    if not path or not os.path.exists(path):
        raise SystemExit(f'slide-code map not found: {path!r} (private file; see metadata.md)')
    with open(path) as fh:
        return {r['real_slide_id']: r['code'] for r in csv.DictReader(fh)}


def discover(cfg, real_ids=False):
    """All slides as dicts (slide_id, real_slide_id, folder, svs, geojson). slide_id is the neutral
    code from the private map; real_ids=True skips the map (used only to build it)."""
    root = cfg['paths']['data_dir']
    skip = set(cfg.get('skip_folders', []))
    slides, problems = [], []
    for folder in sorted(glob.glob(os.path.join(glob.escape(root), '*'))):
        if not os.path.isdir(folder) or os.path.basename(folder) in skip:
            continue
        svs = _visible(glob.glob(os.path.join(glob.escape(folder), '*.svs')))
        gjs = _visible(glob.glob(os.path.join(glob.escape(folder), '*.geojson')))
        if len(svs) != 1 or len(gjs) != 1:
            problems.append(f'{folder!r}: {len(svs)} svs, {len(gjs)} geojson')
            continue
        slides.append(dict(slide_id=slide_id_of(svs[0]), folder=folder, svs=svs[0], geojson=gjs[0]))
    if problems:
        raise SystemExit('Slide folders that are not exactly 1 svs + 1 geojson:\n  ' + '\n  '.join(problems))
    ids = [s['slide_id'] for s in slides]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise SystemExit(f'duplicate slide ids: {dup}')
    for s in slides:
        s['real_slide_id'] = s['slide_id']
    if real_ids:
        return slides
    codes = slide_codes(cfg)
    missing = [s['real_slide_id'] for s in slides if s['real_slide_id'] not in codes]
    if missing:
        raise SystemExit(f'{len(missing)} slide(s) have no code in {cfg["paths"]["slide_codes"]!r}; add them there')
    for s in slides:
        s['slide_id'] = codes[s['real_slide_id']]
    return sorted(slides, key=lambda s: s['slide_id'])


def scrub_text(text, s):
    """Replace a slide's real names (svs / geojson file names, real id) in text with its code."""
    for real in sorted({os.path.basename(s['svs']), os.path.basename(s['geojson']),
                        os.path.splitext(os.path.basename(s['svs']))[0],
                        os.path.splitext(os.path.basename(s['geojson']))[0], s['real_slide_id']}, key=len, reverse=True):
        if real:
            text = text.replace(real, s['slide_id'])
    return text


def scrub_files(folder, s, exts=('.csv', '.json', '.md', '.txt')):
    """scrub_text on every text file under folder (in place). Returns the number of files changed."""
    n = 0
    for root, _, files in os.walk(folder):
        for f in files:
            if f.endswith(exts):
                p = os.path.join(root, f)
                t = open(p, encoding='utf-8', errors='surrogateescape').read()
                u = scrub_text(t, s)
                if u != t:
                    open(p, 'w', encoding='utf-8', errors='surrogateescape').write(u)
                    n += 1
    return n


def write_manifest(slides, path):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, ['index', 'slide_id', 'svs', 'geojson', 'folder'])
        w.writeheader()
        for i, s in enumerate(slides):
            w.writerow(dict(index=i, **s))


def read_manifest(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def get_slide(cfg, slide_id=None, index=None):
    slides = discover(cfg)
    if index is not None:
        return slides[int(index)]
    for s in slides:
        if s['slide_id'] == slide_id:
            return s
    raise SystemExit(f'slide {slide_id!r} not found; known: {[s["slide_id"] for s in slides]}')


# ---------------------------------------------------------------- GeoJSON
def _class_of(feature):
    c = (feature.get('properties') or {}).get('classification')
    return c.get('name') if isinstance(c, dict) else c


def lesions_by_name(geo):
    """Lesions exactly as measure_lesions.py (main, lines 206-222) defines them:
    class or name 'HaemorrhagicLesion' (spaces/case ignored); outlines sharing a QuPath
    name (case-insensitive) are ONE lesion. Returns [(name, [polygon, ...])] sorted by name,
    so lesion index i here == lesion_id i in lesions_measured.csv."""
    names = ('haemorrhagiclesion', 'hemorrhagiclesion', 'lesion')
    found = []
    for ft in geo.get('features', []):
        pr = ft.get('properties', {})
        c = ((pr.get('classification') or {}).get('name', '') or pr.get('name', '')).replace(' ', '').lower()
        if c in names:
            g = ft['geometry']
            polys = [g['coordinates']] if g['type'] == 'Polygon' else g['coordinates'] if g['type'] == 'MultiPolygon' else []
            if polys:
                found.append((pr.get('name', '').strip(), polys))
    byname = {}
    for name, polys in found:
        key = name.lower() or f'unnamed {len(byname) + 1}'
        byname.setdefault(key, [name, []])[1].extend(polys)
    out = [(v[0] or k, v[1]) for k, v in byname.items()]
    out.sort(key=lambda t: t[0].lower())
    return out


def _shoelace(ring):
    c = np.asarray(ring, float)
    return abs(np.dot(c[:, 0], np.roll(c[:, 1], 1)) - np.dot(c[:, 1], np.roll(c[:, 0], 1))) / 2


def rbc_cells(geo, mpp, cls='RBC'):
    """Individually outlined red cells -> [{cx, cy, area_um2}] in slide px.
    Same rule as handover scripts/make_cells.py: Polygon/MultiPolygon outer rings only,
    rings of area <= 1 um2 dropped (slivers); LineStrings ignored."""
    out = []
    for f in geo.get('features', []):
        if _class_of(f) != cls:
            continue
        g = f['geometry']
        if g['type'] == 'Polygon':
            rings = [g['coordinates'][0]]
        elif g['type'] == 'MultiPolygon':
            rings = [p[0] for p in g['coordinates']]
        else:
            continue
        for ring in rings:
            r = np.asarray(ring, float)
            area = _shoelace(r) * mpp * mpp
            if area > 1.0:
                out.append(dict(cx=float(r[:, 0].mean()), cy=float(r[:, 1].mean()), area_um2=float(area),
                                ring=r.tolist()))
    return out


def load_geojson(path):
    with open(path) as fh:
        return json.load(fh)
