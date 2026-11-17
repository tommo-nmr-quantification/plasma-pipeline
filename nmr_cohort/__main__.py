import argparse
from pathlib import Path
from . import core


def main():
    p=argparse.ArgumentParser(description='Frozen 50-metabolite NMR quantification')
    sub=p.add_subparsers(dest='command',required=True)
    for command in ['run','predict']:
        q=sub.add_parser(command,help='Quantify and filter all 50 metabolites' if command=='run' else 'Quantification only')
        q.add_argument('--bin',required=True)
        q.add_argument('--axis',default=str(core.ROOT/'reference/BinTable_Header.tsv'))
        q.add_argument('--models',default=str(core.ROOT/'config/models.json'))
        q.add_argument('--format',choices=['raw8','header','none'],required=True)
        q.add_argument('--id-mode',choices=['full','legacy'],default='full')
        q.add_argument('--reference-column',type=int)
        q.add_argument('--groups')
        q.add_argument('--batch-map')
        q.add_argument('--batch-size',type=int,default=64)
        q.add_argument('--on-fit-error',choices=['raise','na'],default='raise')
        q.add_argument('--out',required=True)
        if command=='run':
            q.add_argument('--reference-stats',default=str(core.ROOT/'config/reference_stats.tsv'))
            q.add_argument('--cutoffs',default=str(core.ROOT/'config/reporting_cutoffs.json'))
            q.add_argument('--batch-parameters')
    q=sub.add_parser('calibrate-batches',help='Freeze additive median offsets from a designated population')
    q.add_argument('--concentrations',required=True)
    q.add_argument('--batch-map',required=True)
    q.add_argument('--rules',required=True)
    q.add_argument('--min-n',type=int,default=10)
    q.add_argument('--input-stage',default='model_raw')
    q.add_argument('--out',required=True)
    q=sub.add_parser('postprocess',help='Apply batch offsets, reference 10SD exclusion, then low cutoffs')
    q.add_argument('--concentrations',required=True)
    q.add_argument('--reference-stats',default=str(core.ROOT/'config/reference_stats.tsv'))
    q.add_argument('--cutoffs',default=str(core.ROOT/'config/reporting_cutoffs.json'))
    q.add_argument('--batch-parameters')
    q.add_argument('--batch-map')
    q.add_argument('--input-stage',default='model_raw')
    q.add_argument('--out',required=True)
    args=p.parse_args()
    try:
        if args.command=='run':
            from .pipeline import predict
            from .postprocess import apply, reference_statistics, checked_mask_rules
            fields=['SampleID']+['TCN%06d'%k for k in range(1,51)]
            reference_statistics(args.reference_stats,fields)
            checked_mask_rules(args.cutoffs,fields)
            out=Path(args.out)
            if out.exists():raise FileExistsError(out)
            prediction_args=argparse.Namespace(**vars(args));prediction_args.out=str(out/'quantified')
            predict(prediction_args)
            post_args=argparse.Namespace(concentrations=str(out/'quantified/raw.tsv'),
                reference_stats=args.reference_stats,cutoffs=args.cutoffs,batch_parameters=args.batch_parameters,
                batch_map=args.batch_map if args.batch_parameters else None,
                input_stage='model_raw',out=str(out/'filtered'))
            result=apply(post_args)
        elif args.command=='predict':
            from .pipeline import predict
            result=predict(args)
        elif args.command=='calibrate-batches':
            from .postprocess import calibrate
            result=calibrate(args)
        elif args.command=='postprocess':
            from .postprocess import apply
            result=apply(args)
        print('Completed:',args.out)
    except (ValueError,KeyError,OSError) as exc:
        p.exit(2,'Error: '+str(exc)+'\n')


if __name__=='__main__':main()
