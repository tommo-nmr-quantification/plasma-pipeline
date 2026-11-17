"""Refresh source-file checksums after intentional repository edits."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = {'nmr_cohort', 'config', 'reference', 'weights', 'docs',
               'examples', 'tests', 'scripts', '.github'}
TOP_FILES = {'README.md', 'README.ja.md', 'requirements.txt', 'environment.yml',
             '.gitignore', 'LICENSE', 'CITATION.cff'}


def update():
    previous = json.loads((ROOT/'MANIFEST.json').read_text())
    files = {}
    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT)
        if '__pycache__' in rel.parts or path.name == '.DS_Store' or path.suffix in {'.pyc', '.pyo'}:
            continue
        selected = (len(rel.parts) == 1 and rel.name in TOP_FILES) or (len(rel.parts) > 1 and rel.parts[0] in SOURCE_DIRS)
        if not selected or not path.is_file():
            continue
        if path.is_symlink():
            raise ValueError('Source files must not be symlinks: ' + str(rel))
        files[str(rel)] = hashlib.sha256(path.read_bytes()).hexdigest()
    previous['files'] = files
    (ROOT/'MANIFEST.json').write_text(json.dumps(previous, indent=2)+'\n')
    print('Updated manifest for %d source files' % len(files))


if __name__ == '__main__':
    update()
