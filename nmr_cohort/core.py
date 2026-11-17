"""Input validation, explicit spectral coordinates, and scalar features."""
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def axis_hash(axis):
    return hashlib.sha256(('\n'.join(float(x).hex() for x in axis) + '\n').encode()).hexdigest()


def read_axis(path):
    tokens = Path(path).read_text().split()
    if tokens and tokens[0] == 'exID':
        tokens = tokens[1:]
    axis = tuple(float(x) for x in tokens)
    if (len(axis) < 2 or not all(math.isfinite(x) for x in axis) or
            not all(a > b for a, b in zip(axis, axis[1:]))):
        raise ValueError('Expected a finite, strictly descending ppm axis')
    return axis


def indices(axis, roi, selector):
    hi, lo = roi['high'], roi['low']
    if not (math.isfinite(hi) and math.isfinite(lo) and axis[0] >= hi > lo >= axis[-1]):
        raise ValueError('ROI must be within the model axis and high > low')
    if selector == 'inclusive':
        chosen = [i for i, p in enumerate(axis) if lo <= p <= hi]
        if not chosen:
            raise ValueError('No bins within ROI')
        return chosen[0], chosen[-1]
    if selector == 'bracketed':
        left = next((i - 1 for i, p in enumerate(axis) if p <= hi), None)
        right = next((i for i, p in enumerate(axis) if p < lo), None)
    elif selector == 'mlp_bracketed':
        left = max(i for i, p in enumerate(axis) if p >= hi)
        right = next((i for i, p in enumerate(axis) if p <= lo), None)
    else:
        raise ValueError('Unknown ROI selector: ' + str(selector))
    if left is None or right is None or left < 0 or left > right:
        raise ValueError('ROI cannot be bracketed by the model axis')
    return left, right


def serial_sum(values):
    total = 0.0
    for value in values:
        total += value
    return total


def feature_plan(model, axis, override=None):
    if override is not None:
        if model['kind'] != 'linear':
            raise ValueError('ROI overrides apply only to scalar linear models')
        a, b = indices(axis, override['roi_ppm'], override.get('selector', model['selector']))
    else:
        a, b = model['bin_start'], model['bin_end']
    ref = indices(axis, model['reference_roi'], model['selector']) if model['feature'] == 'max_minus_reference' else None
    return a, b, ref


def scalar_feature(values, model, plan):
    a, b, ref = plan
    selected = values[a:b + 1]
    kind = model['feature']
    if kind == 'max':
        return max(selected)
    if kind == 'sum':
        return serial_sum(selected)
    if kind == 'mean':
        return serial_sum(selected) / len(selected)
    if kind == 'max_minus_reference':
        c, d = ref
        return max(selected) - max(values[c:d + 1]) / model['reference_divisor']
    raise ValueError('Unsupported scalar feature')


def sample_id(value, mode):
    ident = value if mode == 'full' else value.split('.')[0].split('_')[-1]
    if not ident:
        raise ValueError('Empty sample ID')
    return ident


def spectra(path, axis, fmt, id_mode='full'):
    skip = {'raw8': 8, 'header': 1, 'none': 0}[fmt]
    seen = set()
    with Path(path).open(newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for _ in range(skip):
            row = next(reader, None)
            if row is None:
                raise ValueError('Input shorter than declared header/preamble')
            if fmt == 'header' and axis_hash([float(x) for x in row[1:]]) != axis_hash(axis):
                raise ValueError('BIN header does not match the supplied model axis')
        for row in reader:
            if len(row) != len(axis) + 1:
                raise ValueError('BIN row must contain ID plus %d spectral values' % len(axis))
            ident = sample_id(row[0], id_mode)
            if ident in seen:
                raise ValueError('Duplicate/colliding sample ID: ' + ident)
            seen.add(ident)
            values = [float(x) for x in row[1:]]
            if not all(math.isfinite(x) for x in values):
                raise ValueError('Non-finite spectral value: ' + ident)
            yield ident, values
    if not seen:
        raise ValueError('No spectra')


def read_table(path, required):
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or not set(required) <= set(fields):
            raise ValueError('Missing or duplicate table columns')
        rows = list(reader)
    seen = set()
    if not rows:
        raise ValueError('Empty input table')
    for row in rows:
        if None in row or any(x is None for x in row.values()):
            raise ValueError('Inconsistent table width')
        ident = row['SampleID']
        if not ident or ident in seen:
            raise ValueError('Empty or duplicate sample ID')
        seen.add(ident)
    return fields, rows


def number(text):
    if text in ('', 'NA'):
        return None
    value = float(text)
    if not math.isfinite(value):
        raise ValueError('Use NA for missing values; infinity and NaN are not accepted')
    return value


def text_number(value):
    return 'NA' if value is None else format(value, '.17g')


def load_models(path, axis):
    cfg = load_json(path)
    if (cfg.get('schema_version') != 1 or cfg.get('unit') != 'uM' or
            cfg['axis_bin_n'] != len(axis) or cfg['axis_sha256'] != axis_hash(axis)):
        raise ValueError('Model schema, unit or axis mismatch')
    models = cfg['models']
    if [m['tcn'] for m in models] != ['TCN%06d' % k for k in range(1, 51)]:
        raise ValueError('The catalogue must contain exactly TCN000001–TCN000050 in order')
    for m in models:
        a, b = m['bin_start'], m['bin_end']
        if not 0 <= a <= b < len(axis) or b - a + 1 != m['bin_n']:
            raise ValueError('Invalid frozen feature offsets')
        if (axis[a], axis[b]) != (m['actual_ppm_high'], m['actual_ppm_low']):
            raise ValueError('Frozen ROI coordinates do not match axis')
        if m['kind'] == 'adaptive_linear':
            if m['tcn'] != 'TCN000042' or (a, b) != (25295, 25537):
                raise ValueError('Unsupported adaptive-feature model')
        elif indices(axis, m['roi_ppm'], m['selector']) != (a, b):
            raise ValueError('Selector and frozen feature offsets differ')
        if m['kind'] in ('linear', 'adaptive_linear'):
            if not all(math.isfinite(m[k]) for k in ('slope', 'intercept')):
                raise ValueError('Non-finite coefficient')
        elif m['kind'] == 'mlp':
            weights = ROOT / m['weights_file']
            if digest(weights) != m['weights_sha256']:
                raise ValueError('MLP weight checksum mismatch')
            with np.load(weights, allow_pickle=False) as z:
                if set(z.files) != {'W1', 'b1', 'W2', 'b2'}:
                    raise ValueError('Expected four numeric MLP arrays')
                arrays = {k: z[k].copy() for k in z.files}
            n, h = m['input_dim'], m['hidden_dim']
            if n != b - a + 1:
                raise ValueError('MLP input dimension mismatch')
            for key, shape in {'W1': (n, h), 'b1': (h,), 'W2': (h, 1), 'b2': (1,)}.items():
                if arrays[key].shape != shape or arrays[key].dtype.name != 'float64' or not np.isfinite(arrays[key]).all():
                    raise ValueError('Invalid MLP array: ' + key)
            m['_weights'] = arrays
        else:
            raise ValueError('Unsupported model kind')
    return cfg, models
