# -*- coding: utf-8 -*-
"""
Created on Thu Feb 04 09:32:49 2016

@author: abrasaldo.pmb
"""

from pyTETRAD import *
import pickle
import time
import unittest
from pandas.util.testing import assert_frame_equal

class TestPyTetrad(unittest.TestCase):
    
    def setUp(self):
        self.gv = pickle.load(open('test\\gv.p', 'r'))
        self.isim = pickle.load(open('test\\is.p', 'r'))
        self.out = pickle.load(open('test\\out.p', 'r'))
        self.plt = pickle.load(open('test\\plt.p', 'r'))
        self.grid_gv = pickle.load(open('test\\grid_gv.p', 'r'))
        self.grid_is = pickle.load(open('test\\grid_is.p', 'r'))
        self.dat = pickle.load(open('test\\dat.p', 'r'))
        

    def test_gv(self):
        gvfile = 'test\\base.gv'
        gv = TetradGridView(gvfile)
        result = {}
        result['_params'] = gv._params
        result['_fullpos'] = gv._fullpos
        result['times'] = gv.times.tolist()
        gv.read_all_data()
        result_param_table = gv.param_table
        self_param_table = self.gv.pop('param_table')
        self.gv['times'] = self.gv['times'].tolist()
        self.assertEqual(result, self.gv)
        self.assertEqual(len(result_param_table.keys()), len(self_param_table.keys()))
        for key in result_param_table.keys():
            self.assertTrue(self.assert_df_equal(result_param_table[key], self_param_table[key]))
        
        
    def test_out(self):
        outfile = 'test\\base.out'
        out = TetradOut(outfile)
        result = {}
        result['_pos'] = out._pos
        result['_recurpos'] = out._recurpos   
        result['recur_indices'] = out.recur_indices
        result['recur_times'] = out.recur_times
        result_well_table_df = out.read_well_table()
        self_well_table_df = self.out.pop('well_table_df')
        self.assertEqual(result, self.out)
        self.assertTrue(self.assert_df_equal(result_well_table_df, self_well_table_df))
    
    def test_is(self):
        isfile = 'test\\base.is'
        isim = TetradInterSim(isfile)
        result = {}
        result['_fullpos'] = isim._fullpos
        result['_params'] = isim._params
        result_data = isim.read_data()
        self_data = self.isim.pop('data')
        self.assertEqual(result, self.isim)
        self.assertTrue(self.assert_df_equal(result_data, self_data))
    
    def test_grid(self):
        gridfile = 'test\\base.is'
        grid = TetradGrid(gridfile)
        result = {}
        result['filetype'] = grid.filetype
        result['grid_spec0'] = grid.grid_spec()[0]
        grid_spec1 = grid.grid_spec()[1]
        for key in grid_spec1:
            grid_spec1[key] = grid_spec1[key].tolist()
        result['grid_spec1'] = grid_spec1
        result_df = result.pop('grid_spec0')
        self_df = self.grid_is.pop('grid_spec0')
        self.assertEqual(result, self.grid_is)
        self.assertTrue(self.assert_df_equal(result_df, self_df))
        
        gridfile = 'test\\base.gv'
        grid = TetradGrid(gridfile)
        result = {}
        result['filetype'] = grid.filetype
        result['grid_spec0'] = grid.grid_spec()[0]
        grid_spec1 = grid.grid_spec()[1]
        for key in grid_spec1:
            grid_spec1[key] = grid_spec1[key].tolist()
        result['grid_spec1'] = grid_spec1
        result_df = result.pop('grid_spec0')
        self_df = self.grid_gv.pop('grid_spec0')
        self.assertEqual(result, self.grid_gv)
        self.assertTrue(self.assert_df_equal(result_df, self_df))
    
    def test_plot(self):
        pltfile = 'test\\base.plt'
        plt = TetradPlotFile(pltfile)
        result = {}
        result['_pos'] = plt._pos
        result['_params'] = plt._params
        result_all_times = plt.read_all_times()
        self_all_times = self.plt.pop('read_all_times')
        self.assertEqual(result, self.plt)
        self.assertTrue(self.assert_df_equal(result_all_times, self_all_times))
    
    def test_dat(self):
        datfile = 'test\\base.dat'
        dat = TetradInput(datfile)
        result = {}
        result['_params'] = dat._params
        result['_fullpos'] = dat._fullpos
        result_production_rates = dat.get_production_rates()
        self_production_rates = self.dat.pop('production_rates')
        self.assertEqual(result, self.dat)
        self.assertTrue(self.assert_df_equal(result_production_rates, self_production_rates))
        
    def assert_df_equal(self, df1, df2):
        try:
            assert_frame_equal(df1.sort(axis=1), df2.sort(axis=1), check_names=True)
            return True
        except (AssertionError, ValueError, TypeError):
            return False

def GenerateTestPyTetrad():
    start_time = time.time()
    generate_test_gv('test\\base.gv')
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()    
    generate_test_out('test\\base.out')
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    generate_test_is('test\\base.is')
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    generate_test_grid_is('test\\base.is')
    generate_test_grid_gv('test\\base.gv')
    print("--- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    generate_test_plot('test\\base.plt')
    print("--- %s seconds ---" % (time.time() - start_time))
    
    start_time = time.time()    
    generate_test_dat('test\\base.dat')
    print("--- %s seconds ---" % (time.time() - start_time))
    
def generate_test_gv(gvfile):
    gv = TetradGridView(gvfile)
    result = {}
    result['_params'] = gv._params
    result['_fullpos'] = gv._fullpos
    result['times'] = gv.times
    gv.read_all_data()
    result['param_table'] = gv.param_table
    pickle.dump(result, open('test\\gv.p', 'wb'))
    
def generate_test_out(outfile):
    out = TetradOut(outfile)
    result = {}
    result['_pos'] = out._pos
    result['_recurpos'] = out._recurpos   
    result['recur_indices'] = out.recur_indices
    result['recur_times'] = out.recur_times
    result['well_table_df'] = out.read_well_table()
    pickle.dump(result, open('test\\out.p', 'wb'))

def generate_test_is(isfile):
    isim = TetradInterSim(isfile)
    result = {}
    result['_fullpos'] = isim._fullpos
    result['_params'] = isim._params
    result['data'] = isim.read_data()
    pickle.dump(result, open('test\\is.p', 'wb'))

def generate_test_grid_gv(gridfile):
    grid = TetradGrid(gridfile)
    result = {}
    result['filetype'] = grid.filetype
    result['grid_spec0'] = grid.grid_spec()[0]
    grid_spec1 = grid.grid_spec()[1]
    for key in grid_spec1:
        grid_spec1[key] = grid_spec1[key].tolist()
    result['grid_spec1'] = grid_spec1
    pickle.dump(result, open('test\\grid_gv.p', 'wb'))
    
def generate_test_grid_is(gridfile):
    grid = TetradGrid(gridfile)
    result = {}
    result['filetype'] = grid.filetype
    result['grid_spec0'] = grid.grid_spec()[0]
    grid_spec1 = grid.grid_spec()[1]
    for key in grid_spec1:
        grid_spec1[key] = grid_spec1[key].tolist()
    result['grid_spec1'] = grid_spec1
    pickle.dump(result, open('test\\grid_is.p', 'wb'))

def generate_test_plot(pltfile):
    plt = TetradPlotFile(pltfile)
    result = {}
    result['_pos'] = plt._pos
    result['_params'] = plt._params
    result['read_all_times'] = plt.read_all_times()
    pickle.dump(result, open('test\\plt.p', 'wb'))
    
def generate_test_dat(datfile):
    dat = TetradInput(datfile)
    result = {}
    result['_params'] = dat._params
    result['_fullpos'] = dat._fullpos
    result['production_rates'] = dat.get_production_rates()
    pickle.dump(result, open('test\\dat.p', 'wb'))

if __name__ == '__main__':
    unittest.main()
#    GenerateTestPyTetrad()