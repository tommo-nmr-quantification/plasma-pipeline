"""Generate artificial spectra and group settings; contains no cohort observations."""
import csv
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from nmr_cohort import core, tcn42


def make(path):
    out=Path(path);out.mkdir(parents=True,exist_ok=False)
    axis=core.read_axis(core.ROOT/'reference/BinTable_Header.tsv')
    start,stop=tcn42.roi_columns(axis,[2.96,2.895])
    x=np.arange(start,stop);quartet=[start+35,start+62,start+138,start+165]
    references=[start+110,quartet[2],quartet[2]+8,quartet[3],quartet[3]+8,start+110]
    groups={}
    with (out/'spectra.tsv').open('x') as f, (out/'batch_map.tsv').open('x') as g:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');bw=csv.writer(g,delimiter='\t',lineterminator='\n')
        w.writerow(['SampleID']+list(axis));bw.writerow(['SampleID','BatchID'])
        for i,ref in enumerate(references):
            values=np.full(len(axis),0.0005+i*0.00001)
            y=np.full(len(x),0.00001)
            for k in quartet:y+=0.0008*np.exp(-((x-k)**2)/(2*2.7**2))
            y+=0.0005*np.exp(-((x-ref)**2)/(2*2.7**2))
            if i==5:y=np.zeros(len(x))
            values[start-1:stop-1]=y
            sample='SYNTHETIC_%03d'%(i+1);batch='DEMO_%d'%(i+1)
            w.writerow([sample]+[format(v,'.17g') for v in values]);bw.writerow([sample,batch])
            groups[batch]={'tcn42_reference_column':ref}
    core.write_json(out/'groups.json',dict(schema_version=1,axis_sha256=core.axis_hash(axis),groups=groups))
    with (out/'reference_stats.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['TCN','mean_uM','sd_uM'])
        for k in range(1,51):
            # Artificial settings exercise both masking stages; no cohort statistics.
            w.writerow(['TCN%06d'%k,6.0 if k==5 else 50.0,0.05 if k==5 else 1000.0])
    print(out)


if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Usage: python examples/make_synthetic_input.py OUTPUT_DIRECTORY')
    make(sys.argv[1])
