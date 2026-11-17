import argparse
from pathlib import Path
import tempfile
import unittest
from nmr_cohort import core, postprocess


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)

    def tearDown(self):self.temp.cleanup()

    def table(self,name,content):
        p=self.root/name;p.write_text(content);return str(p)

    def test_low_cutoff_inclusive_and_sd_strict(self):
        self.assertFalse(postprocess.outside_10sd(30.,10.,2.))
        self.assertFalse(postprocess.outside_10sd(-10.,10.,2.))
        self.assertTrue(postprocess.outside_10sd(30.0001,10.,2.))

    def test_median_frozen_parameters_and_na_propagation(self):
        inp=self.table('cal.tsv','SampleID\tTCN000030\nA1\t10\nA2\t14\nB1\t20\nB2\t24\n')
        mapping=self.table('map.tsv','SampleID\tBatchID\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n')
        rules=self.root/'rules.json';core.write_json(rules,dict(schema_version=1,unit='uM',rules=[dict(tcn='TCN000030',reference_batch='A',target_batches=['B'])]))
        frozen=self.root/'frozen.json'
        result=postprocess.calibrate(argparse.Namespace(rules=str(rules),concentrations=inp,batch_map=mapping,min_n=2,input_stage='model_raw',out=str(frozen)))
        self.assertEqual(result['rules'][0]['batches']['B']['offset'],-10.)
        # The one-observation target subset must use the frozen -10 offset, not align its own median.
        subset=self.table('subset.tsv','SampleID\tTCN000030\nB1\t30\nB2\tNA\n')
        stats=self.table('stats.tsv','TCN\tmean_uM\tsd_uM\nTCN000030\t20\t0.5\n')
        cuts=self.root/'cuts.json';core.write_json(cuts,dict(schema_version=1,unit='uM',rules=[]))
        args=argparse.Namespace(concentrations=subset,cutoffs=str(cuts),batch_parameters=str(frozen),batch_map=mapping,
            out=str(self.root/'out'),input_stage='model_raw',reference_stats=stats)
        postprocess.apply(args)
        _,rows=core.read_table(self.root/'out/concentrations.tsv',['SampleID'])
        self.assertEqual(float(rows[0]['TCN000030']),20.)
        self.assertEqual(rows[1]['TCN000030'],'NA')
        args.input_stage='rounded';args.out=str(self.root/'wrongstage')
        with self.assertRaises(ValueError):postprocess.apply(args)

    def test_mask_flags_counts_and_no_rounding(self):
        inp=self.table('input.tsv','SampleID\tTCN000001\nA\t1\nB\t1.0000001\nC\tNA\n')
        stats=self.table('stats.tsv','TCN\tmean_uM\tsd_uM\nTCN000001\t0\t10\n')
        cfg=self.root/'cut.json';core.write_json(cfg,dict(schema_version=1,unit='uM',rules=[dict(tcn='TCN000001',operator='le',value=1,reason='low')]))
        args=argparse.Namespace(concentrations=inp,cutoffs=str(cfg),batch_parameters=None,batch_map=None,
            out=str(self.root/'out'),input_stage='model_raw',reference_stats=stats)
        r=postprocess.apply(args)
        _,rows=core.read_table(self.root/'out/concentrations.tsv',['SampleID'])
        self.assertEqual([x['TCN000001'] for x in rows],['NA','1.0000001000000001','NA'])
        self.assertEqual(r['input_missing_n'],{'TCN000001':1});self.assertEqual(r['newly_masked_n'],{'TCN000001':1})

    def test_10sd_then_low_cutoff_in_order(self):
        inp=self.table('input.tsv','SampleID\tTCN000001\nA\t-11\nB\t-10\nC\t10\nD\t11\nE\t1\nF\t1.00001\nG\tNA\n')
        stats=self.table('stats.tsv','TCN\tmean_uM\tsd_uM\nTCN000001\t0\t1\n')
        cfg=self.root/'cut.json';core.write_json(cfg,dict(schema_version=1,unit='uM',rules=[dict(tcn='TCN000001',operator='le',value=1,reason='low')]))
        args=argparse.Namespace(concentrations=inp,cutoffs=str(cfg),batch_parameters=None,batch_map=None,
            out=str(self.root/'out'),input_stage='model_raw',reference_stats=stats)
        r=postprocess.apply(args)
        _,stage2=core.read_table(self.root/'out/after_10sd.tsv',['SampleID'])
        self.assertEqual([x['TCN000001'] for x in stage2[:5]],['NA','-10','10','NA','1'])
        _,final=core.read_table(self.root/'out/concentrations.tsv',['SampleID'])
        self.assertEqual([x['TCN000001'] for x in final[:5]],['NA','NA','10','NA','NA'])
        self.assertEqual(r['masked_10sd_n'],{'TCN000001':2})
        self.assertEqual(r['masked_low_n'],{'TCN000001':2})
        self.assertEqual(r['newly_masked_n'],{'TCN000001':4})
        _,flags=core.read_table(self.root/'out/mask_flags.tsv',['SampleID'])
        self.assertEqual({x['SampleID']:x['Stage'] for x in flags},{'A':'10SD','D':'10SD','B':'low_cutoff','E':'low_cutoff'})

    def test_reference_stats_required_for_every_metabolite(self):
        for content in ['TCN\tmean_uM\tsd_uM\nTCN000002\t1\t2\n',
                        'TCN\tmean_uM\tsd_uM\nTCN000001\t1\t-2\n',
                        'TCN\tmean_uM\tsd_uM\nTCN000001\t1\t2\nTCN000001\t1\t2\n']:
            p=self.table('bad.tsv',content)
            with self.assertRaises(ValueError):postprocess.reference_statistics(p,['SampleID','TCN000001'])

    def test_duplicate_ids_and_missing_batch_fail(self):
        duplicate=self.table('duplicate.tsv','SampleID\tTCN000001\nA\t1\nA\t2\n')
        with self.assertRaises(ValueError):core.read_table(duplicate,['SampleID'])
        mapping=self.table('map.tsv','SampleID\tBatchID\nA\tB\n')
        with self.assertRaises(ValueError):postprocess.join_batches([{'SampleID':'Z'}],mapping)


if __name__=='__main__':unittest.main()
