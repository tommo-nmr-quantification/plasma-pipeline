"""Verify the shipped file manifest before running the frozen models."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from nmr_cohort import core


def verify():
    manifest=core.load_json(core.ROOT/'MANIFEST.json')
    for name,expected in manifest['files'].items():
        path=core.ROOT/name
        if path.is_symlink() or core.digest(path)!=expected:
            raise ValueError('Missing or altered release file: '+name)
    axis=core.read_axis(core.ROOT/'reference/BinTable_Header.tsv')
    _,models=core.load_models(core.ROOT/'config/models.json',axis)
    print('Verified %d files, %d models, %d spectral bins'%(len(manifest['files']),len(models),len(axis)))


if __name__=='__main__':verify()
