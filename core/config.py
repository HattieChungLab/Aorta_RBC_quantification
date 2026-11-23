"""Config loading. The YAML path is the only thing scripts need."""
import os
import sys

import yaml

PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CONFIG = os.path.join(PIPELINE_DIR, 'config', 'base.yaml')


def load_config(path=None):
    with open(path or DEFAULT_CONFIG) as fh:
        cfg = yaml.safe_load(fh)
    # the validated handover code is imported, never copied or edited
    hd = cfg['paths']['handover_dir']
    if hd not in sys.path:
        sys.path.insert(0, hd)
    os.environ.setdefault('CELLPOSE_LOCAL_MODELS_PATH', cfg['paths']['cellpose_models_dir'])
    return cfg


def out_dir(cfg, *parts):
    p = os.path.join(cfg['paths']['output_dir'], *parts)
    os.makedirs(p, exist_ok=True)
    return p
