"""Apply the frozen 50-metabolite catalogue; no model training at inference."""
from collections import Counter
import csv
import math
from pathlib import Path
import platform
import warnings
import numpy as np
import pandas as pd
import scipy
from . import core, tcn42


def group_plan(args, models, axis):
    membership = {}
    if args.reference_column is not None:
        if args.batch_map or args.groups:
            raise ValueError('Choose --reference-column OR --batch-map with --groups')
        groups = {'single_batch': {'tcn42_reference_column': args.reference_column}}
    else:
        if not args.batch_map or not args.groups:
            raise ValueError('TCN000042 requires --reference-column or --batch-map with --groups')
        config = core.load_json(args.groups)
        if config.get('schema_version') != 1 or config.get('axis_sha256') != core.axis_hash(axis):
            raise ValueError('Group configuration schema/axis mismatch')
        groups = config['groups']
        _, rows = core.read_table(args.batch_map, ['SampleID', 'BatchID'])
        membership = {r['SampleID']: r['BatchID'] for r in rows}
        if any(b not in groups for b in membership.values()):
            raise ValueError('Batch map contains a group without settings')
    lookup = {m['tcn']: m for m in models}
    m42 = lookup['TCN000042']
    start, stop = m42['bin_start'] + 1, m42['bin_end'] + 2
    plans = {}
    if not groups:
        raise ValueError('No group settings supplied')
    for batch, cfg in groups.items():
        if not batch:
            raise ValueError('Empty batch label')
        tcn42.validate_reference(cfg['tcn42_reference_column'], start, stop)
        overrides = cfg.get('roi_overrides', {})
        if set(overrides) - set(lookup):
            raise ValueError('Unknown metabolite in ROI overrides')
        plans[batch] = {m['tcn']: core.feature_plan(m, axis, overrides.get(m['tcn']))
                        for m in models if m['kind'] != 'adaptive_linear'}
    return membership, groups, plans


def predict(args):
    if args.batch_size < 1:
        raise ValueError('Batch size must be positive')
    axis = core.read_axis(args.axis)
    config, models = core.load_models(args.models, axis)
    membership, groups, plans = group_plan(args, models, axis)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    info = dict(status='running', release=config['release'], unit='uM',
        input_sha256=core.digest(args.bin), models_sha256=core.digest(args.models),
        axis_sha256=core.axis_hash(axis), input_format=args.format, id_mode=args.id_mode,
        batch_size=args.batch_size, on_fit_error=args.on_fit_error,
        batch_map_sha256=core.digest(args.batch_map) if args.batch_map else None,
        groups_sha256=core.digest(args.groups) if args.groups else None,
        reference_column=args.reference_column, python=platform.python_version(),
        numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        weights_sha256={m['tcn']: m['weights_sha256'] for m in models if m['kind']=='mlp'})
    n = 0
    def group(ident):
        if not membership:
            return 'single_batch'
        if ident not in membership:
            raise ValueError('Missing batch assignment for ' + ident)
        return membership[ident]
    try:
        # The v6 amplitude fit retains the original pandas C float parser.
        # Other models use float64 Python conversion and serial summation.
        m42 = next(m for m in models if m['tcn'] == 'TCN000042')
        start, stop = m42['bin_start'] + 1, m42['bin_end'] + 2
        dmg = {}; modes = Counter()
        with (out/'tcn42_diagnostics.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(['SampleID','BatchID','ReferenceColumn','Mode','Status','Height','FallbackReason','Warnings'])
            for raw_id, series in tcn42.spectra_roi(args.bin, axis, args.format, start, stop):
                ident = core.sample_id(raw_id, args.id_mode)
                if ident in dmg:
                    raise ValueError('Duplicate/colliding sample ID: ' + ident)
                batch = group(ident); ref = groups[batch]['tcn42_reference_column']
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    try:
                        result = tcn42.extract_height(series, ref)
                        if result['height'] is not None and not math.isfinite(result['height']):
                            raise ValueError('Non-finite adjusted peak height')
                    except (ValueError, RuntimeError) as exc:
                        if args.on_fit_error == 'raise':
                            raise ValueError('TCN000042 fit failed for ' + ident + ': ' + str(exc)) from exc
                        result = dict(height=None, mode='unresolved', status='fit_error', fallback_reason=str(exc))
                height = result['height']
                value = None if height is None else m42['slope'] * height + m42['intercept']
                if value is not None and not math.isfinite(value):
                    raise ValueError('Non-finite TCN000042 concentration')
                dmg[ident] = (height, value); modes[result['mode']] += 1
                w.writerow([ident,batch,ref,result['mode'],result['status'],core.text_number(height),
                            result['fallback_reason'],'; '.join(str(x.message) for x in caught)])
        with (out/'raw.tsv').open('x') as raw, (out/'features.tsv').open('x') as feat:
            rw = csv.writer(raw, delimiter='\t', lineterminator='\n')
            fw = csv.writer(feat, delimiter='\t', lineterminator='\n')
            rw.writerow(['SampleID'] + [m['tcn'] for m in models])
            scalar_models = [m for m in models if m['kind'] != 'mlp']
            fw.writerow(['SampleID'] + [m['tcn'] for m in scalar_models])
            def emit(batch):
                predictions = {}; features = {}
                for m in models:
                    ident = m['tcn']
                    if m['kind'] == 'mlp':
                        a, b = m['bin_start'], m['bin_end']; weights = m['_weights']
                        x = np.array([v[a:b+1] for _, v in batch], dtype='float64')
                        values = (np.tanh(x @ weights['W1'] + weights['b1']) @ weights['W2'] + weights['b2']).reshape(-1)
                        if not np.isfinite(values).all():
                            raise ValueError('Non-finite MLP prediction')
                        predictions[ident] = list(values)
                    elif m['kind'] == 'adaptive_linear':
                        features[ident] = [dmg[s][0] for s, _ in batch]
                        predictions[ident] = [dmg[s][1] for s, _ in batch]
                    else:
                        x = [core.scalar_feature(v, m, plans[group(s)][ident]) for s, v in batch]
                        values = [m['slope'] * v + m['intercept'] for v in x]
                        if not all(math.isfinite(v) for v in values):
                            raise ValueError('Non-finite linear prediction')
                        features[ident] = x; predictions[ident] = values
                for i, (sample, _) in enumerate(batch):
                    rw.writerow([sample] + [core.text_number(predictions[m['tcn']][i]) for m in models])
                    fw.writerow([sample] + [core.text_number(features[m['tcn']][i]) for m in scalar_models])
            part = []
            for record in core.spectra(args.bin, axis, args.format, args.id_mode):
                part.append(record); n += 1
                if len(part) == args.batch_size:
                    emit(part); part = []
            if part:
                emit(part)
        if n != len(dmg):
            raise ValueError('TCN000042 and standard-model sample counts differ')
        info.update(status='completed',sample_n=n,tcn42_modes=dict(modes),
                    output_sha256={p.name:core.digest(p) for p in out.iterdir() if p.is_file()})
    except Exception as exc:
        info.update(status='failed',sample_n=n,error_type=type(exc).__name__,error=str(exc))
        core.write_json(out/'run.json',info)
        raise
    core.write_json(out/'run.json',info)
    return info
