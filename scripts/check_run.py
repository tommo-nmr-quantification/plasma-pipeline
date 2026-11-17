"""Check completion, shape, alignment and recorded output hashes."""
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from nmr_cohort import core


def check(path):
    root=Path(path);r=json.loads((root/'run.json').read_text())
    if r.get('status')!='completed':raise ValueError('Run did not complete')
    for name,h in r['output_sha256'].items():
        if core.digest(root/name)!=h:raise ValueError('Changed output: '+name)
    target='raw.tsv' if (root/'raw.tsv').exists() else 'concentrations.tsv'
    fields,rows=core.read_table(root/target,['SampleID'])
    if fields!=['SampleID']+['TCN%06d'%k for k in range(1,51)]:raise ValueError('Expected all 50 metabolites')
    if len(rows)!=r['sample_n']:raise ValueError('Row count differs from completion record')
    for row in rows:
        for key in fields[1:]:core.number(row[key])
    if (root/'after_10sd.tsv').exists():
        sf,sr=core.read_table(root/'after_10sd.tsv',['SampleID'])
        if sf!=fields or [x['SampleID'] for x in sr]!=[x['SampleID'] for x in rows]:
            raise ValueError('10SD and final tables differ in shape/order')
        for before,after in zip(sr,rows):
            for tcn in fields[1:]:
                if core.number(after[tcn]) is not None and core.number(after[tcn])!=core.number(before[tcn]):
                    raise ValueError('The low-cutoff stage must only retain a value or set it to NA')
        for tcn in fields[1:]:
            missing=sum(core.number(row[tcn]) is None for row in rows)
            expected=r['input_missing_n'].get(tcn,0)+r['masked_10sd_n'].get(tcn,0)+r['masked_low_n'].get(tcn,0)
            if missing!=expected:raise ValueError('Per-stage exclusion counts do not reconcile')
    print('PASS',root,'samples='+str(len(rows)),'metabolites=50')
    if 'tcn42_modes' in r:print('TCN000042 modes:',r['tcn42_modes'])
    if 'masked_10sd_n' in r:print('Newly masked: 10SD=%d, low cutoffs=%d'%(sum(r['masked_10sd_n'].values()),sum(r['masked_low_n'].values())))
    return [r['SampleID'] for r in rows]


if __name__=='__main__':
    if len(sys.argv)<2:raise SystemExit('Usage: python scripts/check_run.py RUN_DIRECTORY [REPORTED_DIRECTORY]')
    ids=[check(p) for p in sys.argv[1:]]
    if any(x!=ids[0] for x in ids):raise ValueError('Sample order differs between stages')
