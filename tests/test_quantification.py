import argparse
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from nmr_cohort import core, pipeline, tcn42


def fixture_module():
    spec=importlib.util.spec_from_file_location('synthetic',core.ROOT/'examples/make_synthetic_input.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class QuantificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.root=Path(cls.temp.name)
        fixture_module().make(cls.root/'input')
        cls.axis=core.read_axis(core.ROOT/'reference/BinTable_Header.tsv')
        cls.config,cls.models=core.load_models(core.ROOT/'config/models.json',cls.axis)
        cls.lookup={m['tcn']:m for m in cls.models}
        cls.args=argparse.Namespace(bin=str(cls.root/'input/spectra.tsv'),axis=str(core.ROOT/'reference/BinTable_Header.tsv'),
            models=str(core.ROOT/'config/models.json'),format='header',id_mode='full',batch_size=2,
            reference_column=None,groups=str(cls.root/'input/groups.json'),batch_map=str(cls.root/'input/batch_map.tsv'),
            on_fit_error='raise',out=str(cls.root/'predicted'))
        cls.report=pipeline.predict(cls.args)
        _,cls.predictions=core.read_table(cls.root/'predicted/raw.tsv',['SampleID'])

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_catalogue_complete_and_model_kinds(self):
        self.assertEqual(len(self.models),50)
        self.assertEqual(sum(m['kind']=='linear' for m in self.models),45)
        self.assertEqual(sum(m['kind']=='mlp' for m in self.models),4)
        self.assertEqual(self.lookup['TCN000005']['feature'],'max')
        self.assertEqual(self.lookup['TCN000005']['slope'],17178.39)

    def test_actual_deployed_window_sizes(self):
        expected={38:8,39:13,40:13,41:13,43:13,44:33,45:33}
        for k,n in expected.items():self.assertEqual(self.lookup['TCN%06d'%k]['bin_n'],n)
        self.assertEqual((self.axis[18962],self.axis[18963]),(4.87012,4.66978))

    def test_shared_column_table_matches_all_50_executable_models(self):
        with (core.ROOT/'config/model_columns.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual([r['TCN'] for r in rows],list(self.lookup))
        for r in rows:
            m=self.lookup[r['TCN']]
            self.assertEqual(int(r['TSV_first_column_1based']),m['bin_start']+2)
            self.assertEqual(int(r['TSV_last_column_1based']),m['bin_end']+2)
            self.assertEqual(int(r['BIN_count']),m['bin_n'])
            self.assertEqual(float(r['ppm_first']),self.axis[m['bin_start']])
            self.assertEqual(float(r['ppm_last']),self.axis[m['bin_end']])
            if m['kind']!='mlp':
                self.assertEqual(float(r['Slope']),m['slope'])
                self.assertEqual(float(r['Intercept']),m['intercept'])
            if m.get('reference_roi'):
                a,b=core.indices(self.axis,m['reference_roi'],m['selector'])
                self.assertEqual(int(r['Reference_first_column_1based']),a+2)
                self.assertEqual(int(r['Reference_last_column_1based']),b+2)

    def test_new_peak_height_models_closed_form(self):
        # Artificial test input is constant in each of these six remote ROIs.
        coefs={5:(17178.39,-2.587192),46:(16442.26,-5.319722),47:(9605.018,-4.229689),
               48:(85234.44,-.1383509),49:(12225.08,-6.370426),50:(7910.959,-.4111907)}
        for i,row in enumerate(self.predictions):
            for k,(a,b) in coefs.items():
                self.assertAlmostEqual(float(row['TCN%06d'%k]),a*(.0005+i*.00001)+b,places=12)

    def test_peak_height_endpoints_and_outside_spike(self):
        for k in [5,46,47,48,49,50]:
            model=self.lookup['TCN%06d'%k];a,b=model['bin_start'],model['bin_end']
            vals=[0.0]*len(self.axis);vals[a]=2.;vals[b]=3.;vals[a-1]=100.;vals[b+1]=200.
            self.assertEqual(core.scalar_feature(vals,model,core.feature_plan(model,self.axis)),3.)
            self.assertGreaterEqual(self.axis[b],model['roi_ppm']['low'])
            self.assertLessEqual(self.axis[a],model['roi_ppm']['high'])

    def test_peak_overlap_modes_and_unresolved(self):
        self.assertEqual(self.report['tcn42_modes'],{'isolated':1,'overlap3_subtraction':1,
            'overlap3_partial':1,'overlap4_subtraction':1,'overlap4_partial':1,'unresolved':1})
        _,diag=core.read_table(self.root/'predicted/tcn42_diagnostics.tsv',['SampleID'])
        for row in diag[:-1]:self.assertAlmostEqual(float(row['Height']),.0005,delta=2e-6)
        self.assertEqual(self.predictions[-1]['TCN000042'],'NA')
        self.assertEqual(self.report['sample_n'],6)
        self.assertEqual(len(self.predictions[0]),51)

    def test_mlp_exported_theano_synthetic_probes(self):
        for k in [7,24,32,34]:
            model=self.lookup['TCN%06d'%k];w=model['_weights']
            with np.load(core.ROOT/'tests/fixtures'/('TCN%06d.probe.npz'%k),allow_pickle=False) as data:
                actual=(np.tanh(data['X']@w['W1']+w['b1'])@w['W2']+w['b2']).reshape(-1)
                np.testing.assert_allclose(actual,data['y'],rtol=1e-8,atol=1e-10)

    def test_group_override_is_local(self):
        m=self.lookup['TCN000028'];normal=core.feature_plan(m,self.axis)
        adjusted=core.feature_plan(m,self.axis,dict(roi_ppm=dict(high=3.957,low=3.948)))
        self.assertNotEqual(normal,adjusted)
        self.assertEqual(core.feature_plan(m,self.axis),normal)
        with self.assertRaises(ValueError):core.feature_plan(self.lookup['TCN000007'],self.axis,dict(roi_ppm=dict(high=2.,low=1.)))

    def test_header_axis_width_and_duplicate_ids_rejected(self):
        small=(3.,2.,1.)
        for content in ['ID\t3\t2\t0\nS\t1\t2\t3\n','ID\t3\t2\t1\nS\t1\t2\n',
                        'ID\t3\t2\t1\nS\t1\t2\t3\nS\t1\t2\t3\n',
                        'ID\t3\t2\t1\nS\t1\tnan\t3\n']:
            p=self.root/'invalid.tsv';p.write_text(content)
            with self.assertRaises(ValueError):list(core.spectra(p,small,'header'))

    def test_input_format_choices(self):
        for fmt,head in [('none',''),('raw8','metadata\n'*8),('header','ID\t3\t2\t1\n')]:
            p=self.root/'format.tsv';p.write_text(head+'sample\t1\t2\t3\n')
            self.assertEqual(list(core.spectra(p,(3.,2.,1.),fmt)),[('sample',[1.,2.,3.])])

    def test_reference_convention_and_failed_gaussian(self):
        self.assertEqual(tcn42.roi_columns(self.axis,[2.96,2.895]),(25296,25539))
        with self.assertRaises(ValueError):tcn42.validate_reference(25296,25296,25539)
        with self.assertRaises(ValueError):tcn42.validate_reference(2.93,25296,25539)
        # A partial-overlap fit failure must retain the documented subtraction fallback.
        from unittest.mock import patch
        start,stop=25296,25539;x=np.arange(start,stop);q=[start+35,start+62,start+138,start+165]
        y=np.full(len(x),.00001)
        for k in q:y+=.0008*np.exp(-(x-k)**2/(2*2.7**2))
        ref=q[2]+8;y+=.0005*np.exp(-(x-ref)**2/(2*2.7**2))
        with patch.object(tcn42,'double_gaussian_fit',side_effect=RuntimeError('fixture failure')):
            result=tcn42.extract_height(pd.Series(y,index=x),ref)
        self.assertEqual(result['mode'],'overlap3_subtraction')
        self.assertIn('fixture failure',result['fallback_reason'])

    def test_existing_output_not_overwritten(self):
        with self.assertRaises(FileExistsError):pipeline.predict(self.args)


if __name__=='__main__':unittest.main()
