"""Explicit artifact writes; nothing is saved during import or model creation."""
from pathlib import Path
import json
import numpy as np


def write_json(path,value):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2)+"\n")
    temporary.replace(path)


def save_band_data(path,model,currents):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,kx=model.kx,ky=model.ky,bands=model.bands(currents),currents=currents)
