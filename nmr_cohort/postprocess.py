"""Optional batch correction, then reference 10SD exclusion, then low cutoffs."""
from collections import Counter
import csv
import math
from pathlib import Path
import statistics
from . import core


def join_batches(rows, path):
    _, table = core.read_table(path, ['SampleID', 'BatchID'])
    mapping = {r['SampleID']: r['BatchID'] for r in table}
    if any(not r['BatchID'] for r in table):
        raise ValueError('Empty batch label')
    if any(r['SampleID'] not in mapping for r in rows):
        raise ValueError('Missing batch assignment')
    return mapping


def calibrate(args):
    cfg = core.load_json(args.rules)
    if cfg.get('schema_version') != 1 or cfg.get('unit') != 'uM' or args.min_n < 2:
        raise ValueError('Invalid batch calibration schema/unit/minimum N')
    rules = cfg['rules']
    if not rules or len({r['tcn'] for r in rules}) != len(rules):
        raise ValueError('Expected one calibration rule per metabolite')
    _, rows = core.read_table(args.concentrations, ['SampleID'] + [r['tcn'] for r in rules])
    mapping = join_batches(rows, args.batch_map)
    frozen = []
    for rule in rules:
        ref = rule['reference_batch']
        batches = [ref] + rule['target_batches']
        unchanged = rule.get('unchanged_batches', [])
        if len(set(batches + unchanged)) != len(batches + unchanged):
            raise ValueError('Batch groups must be disjoint')
        if set(mapping[r['SampleID']] for r in rows) - set(batches + unchanged):
            raise ValueError('Undeclared batch in calibration population')
        values = {b: [] for b in batches}; missing = Counter()
        for row in rows:
            batch = mapping[row['SampleID']]
            if batch not in values:
                continue
            v = core.number(row[rule['tcn']])
            if v is None:
                missing[batch] += 1
            else:
                values[batch].append(v)
        if any(len(x) < args.min_n for x in values.values()):
            raise ValueError('Insufficient calibration observations')
        medians = {b: statistics.median(values[b]) for b in batches}
        frozen.append(dict(tcn=rule['tcn'], reference_batch=ref, unchanged_batches=unchanged,
            batches={b:dict(median=medians[b],n=len(values[b]),missing_n=missing[b],
                           offset=medians[ref]-medians[b]) for b in batches}))
    result = dict(schema_version=1,unit='uM',input_stage=args.input_stage,
        method='additive_median',minimum_n=args.min_n,
        formula='corrected = input + median(reference) - median(batch)',
        source_sha256=dict(concentrations=core.digest(args.concentrations),batch_map=core.digest(args.batch_map),rules=core.digest(args.rules)),
        rules=frozen)
    core.write_json(args.out,result)
    return result


def reference_statistics(path, fields):
    """Read fixed reference means/SDs; never estimate them from the target set."""
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != ['TCN', 'mean_uM', 'sd_uM']:
            raise ValueError('Reference columns must be TCN, mean_uM, sd_uM')
        stats = {}
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError('Invalid reference table width')
            tcn = row['TCN']
            if not tcn or tcn in stats:
                raise ValueError('Empty or duplicate TCN in reference statistics')
            mean, sd = float(row['mean_uM']), float(row['sd_uM'])
            if not math.isfinite(mean) or not math.isfinite(sd) or sd < 0:
                raise ValueError('Reference mean/SD must be finite; SD must be nonnegative')
            stats[tcn] = (mean, sd)
    if set(fields) - {'SampleID'} - set(stats):
        raise ValueError('A reference mean and SD are required for every metabolite')
    return stats


def checked_mask_rules(path, fields):
    cfg = core.load_json(path)
    if cfg.get('schema_version') != 1 or cfg.get('unit') != 'uM':
        raise ValueError('Reporting-rule schema/unit mismatch')
    rules = cfg['rules']
    for r in rules:
        if r['tcn'] not in fields or not r.get('reason'):
            raise ValueError('Unknown metabolite or missing reporting-rule reason')
        if r['operator'] != 'le' or not math.isfinite(r['value']):
            raise ValueError('Low-value rules must specify a finite inclusive cutoff (le)')
    if len({r['tcn'] for r in rules}) != len(rules):
        raise ValueError('Duplicate metabolite in low-value rules')
    return rules


def outside_10sd(value, mean, sd):
    return value < mean - 10.0 * sd or value > mean + 10.0 * sd


def write_concentrations(path, fields, rows, values, display=False):
    with Path(path).open('x') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(fields)
        for row in rows:
            v = values[row['SampleID']]
            w.writerow([row[k] if k == 'SampleID' else ('NA' if v[k] is None else
                format(v[k], '.2f') if display else core.text_number(v[k])) for k in fields])


def apply(args):
    fields, rows = core.read_table(args.concentrations, ['SampleID'])
    stats = reference_statistics(args.reference_stats, fields)
    rules = checked_mask_rules(args.cutoffs, fields)
    values = {r['SampleID']:{k:core.number(r[k]) for k in fields if k!='SampleID'} for r in rows}
    input_missing = Counter({k:sum(v[k] is None for v in values.values()) for k in fields if k!='SampleID'})
    offsets = {}; mapping = {}; batch_sha = None
    if args.batch_parameters:
        if not args.batch_map:
            raise ValueError('Batch parameters require --batch-map')
        params = core.load_json(args.batch_parameters)
        if (params.get('schema_version')!=1 or params.get('unit')!='uM' or
            params.get('method')!='additive_median' or params.get('input_stage')!=args.input_stage):
            raise ValueError('Batch parameters do not match schema/unit/input stage')
        mapping = join_batches(rows,args.batch_map)
        for r in params['rules']:
            tcn = r['tcn']
            if tcn not in fields or tcn in offsets:
                raise ValueError('Unknown or duplicate metabolite in batch parameters')
            reference = r['batches'][r['reference_batch']]['median']
            offsets[tcn] = {}
            for b, p in r['batches'].items():
                if not math.isfinite(p['offset']) or not math.isfinite(p['median']) or p['offset'] != reference-p['median']:
                    raise ValueError('Inconsistent frozen median offset')
                offsets[tcn][b] = p['offset']
            for b in r['unchanged_batches']:
                if b in offsets[tcn]:
                    raise ValueError('Unchanged and calibrated groups must be disjoint')
                offsets[tcn][b] = 0.0
            if any(mapping[row['SampleID']] not in offsets[tcn] for row in rows):
                raise ValueError('No frozen offset/unchanged setting for input batch')
        batch_sha = core.digest(args.batch_parameters)
    elif args.batch_map:
        raise ValueError('--batch-map requires --batch-parameters')
    out = Path(args.out); out.mkdir(parents=True,exist_ok=False)
    flags = []; sd_masked = Counter(); low_masked = Counter()
    # Batch correction is part of concentration harmonization, before masking.
    for row in rows:
        sample = row['SampleID']; vals = values[sample]
        for tcn, value in list(vals.items()):
            if value is None:
                continue
            if tcn in offsets:
                value += offsets[tcn][mapping[sample]]
            if not math.isfinite(value):
                raise ValueError('Non-finite post-correction concentration')
            vals[tcn] = value
        # Stage 2: compare with the fixed reference distribution, not this input's SD.
        for tcn, value in list(vals.items()):
            if value is not None and outside_10sd(value, *stats[tcn]):
                flags.append([sample,tcn,'10SD','outside_reference_10sd',core.text_number(value)])
                vals[tcn] = None; sd_masked[tcn] += 1
    write_concentrations(out/'after_10sd.tsv', fields, rows, values)
    # Stage 3: already-missing entries are not evaluated again for the low cutoff.
    for row in rows:
        sample = row['SampleID']; vals = values[sample]
        for rule in rules:
            tcn = rule['tcn']; value = vals[tcn]
            if value is not None and value <= rule['value']:
                flags.append([sample,tcn,'low_cutoff',rule['reason'],core.text_number(value)])
                vals[tcn] = None; low_masked[tcn] += 1
    for name, display in [('concentrations.tsv',False),('display_2dp.tsv',True)]:
        write_concentrations(out/name, fields, rows, values, display)
    with (out/'mask_flags.tsv').open('x') as f:
        w = csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['SampleID','TCN','Stage','Reason','ValueBeforeRule']);w.writerows(flags)
    summary = dict(status='completed',unit='uM',sample_n=len(rows),input_stage=args.input_stage,
        order=['frozen_batch_offsets','reference_10sd','low_cutoffs','display_formatting'],
        input_sha256=core.digest(args.concentrations),cutoffs_sha256=core.digest(args.cutoffs),
        reference_stats_sha256=core.digest(args.reference_stats),sd_multiplier=10,
        batch_parameters_sha256=batch_sha,batch_map_sha256=core.digest(args.batch_map) if args.batch_map else None,
        input_missing_n=dict(input_missing),masked_10sd_n=dict(sd_masked),masked_low_n=dict(low_masked),
        newly_masked_n=dict(sd_masked+low_masked),
        output_sha256={p.name:core.digest(p) for p in out.iterdir() if p.is_file()})
    core.write_json(out/'run.json',summary)
    return summary
