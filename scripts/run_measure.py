#!/usr/bin/env python3
"""STEP 2 (the study's numbers) for ONE slide: runs the validated measure_lesions.py,
unchanged, as a subprocess - exactly what run_all.py did, minus its '_reviewed' file-name
matching (each slide folder here holds exactly one geojson).

    python run_measure.py --index N [--config CONFIG.yaml] [--redo]

Output: outputs/step2/<slide_id>/  (lesions_measured.csv, summary_measured.json/.csv, review/)
A slide whose summary_measured.json exists is skipped unless --redo.
"""
import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_config, out_dir   # noqa: E402
from core.slides import get_slide, scrub_files               # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--index', type=int)
ap.add_argument('--slide-id')
ap.add_argument('--config')
ap.add_argument('--redo', action='store_true')
a = ap.parse_args()
cfg = load_config(a.config)
s = get_slide(cfg, slide_id=a.slide_id, index=a.index)
out = out_dir(cfg, 'step2', s['slide_id'])
if os.path.exists(os.path.join(out, 'summary_measured.json')) and not a.redo:
    print(f'[{s["slide_id"]}] already done, skipped (use --redo)')
    sys.exit(0)

step2 = os.path.join(cfg['paths']['handover_dir'], 'measure_lesions.py')
print(f'[{s["slide_id"]}] svs={s["svs"]!r}\n  geojson={s["geojson"]!r}\n  out={out!r}', flush=True)
t = time.time()
r = subprocess.run([sys.executable, '-W', 'ignore', step2, s['svs'], '--geojson', s['geojson'], '--out', out])
if r.returncode == 0:
    scrub_files(out, s)          # outputs name the slide by its code only, not its file names
print(f'[{s["slide_id"]}] step 2 {"done" if r.returncode == 0 else "FAILED"} in {time.time() - t:.0f} s', flush=True)
sys.exit(r.returncode)
