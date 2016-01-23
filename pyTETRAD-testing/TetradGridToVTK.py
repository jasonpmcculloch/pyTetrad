# -*- coding: utf-8 -*-
"""
JYTA
Last Edit: 1/22/2016
"""

from pyTETRAD import *
import numpy as np 
import pandas as pd


def WriteGridPts(fs, vals, axis): 
    no = vals.shape[0] 
    print >> fs, axis + '_COORDINATES',  no, 'float' 
    for i, v in enumerate(vals): 
        print >> fs, v, 
        if (i % 8 == 7): 
            print >> fs 
    print >> fs 


def WritePointData(fs, name, vals): 
    print >> fs, 'SCALARS',   name, 'float' 
    print >> fs, 'LOOKUP_TABLE default' 
    noX, noY, noZ = vals.shape 
    last = True 
    for k in range(noZ): 
        for j in range(noY): 
            for i in range(noX): 
                print >> fs, vals[i,j,k], 
                if (i % 8 == 7): 
                    print >> fs 
                    last = True 
                else: 
                    last = False 
            if (not last): 
                print >> fs 


def DoVTKLegacy(x, y, z, data, fn, time = None): 
    with open(fn, 'w') as fs: 
        print >> fs, '# vtk DataFile Version 2.0' 
        print >> fs, 'Rectilinear Grid Example' 
        print >> fs, 'ASCII' 
        print >> fs, 'DATASET RECTILINEAR_GRID' 
        print >> fs, 'DIMENSIONS ', x.shape[0], y.shape[0], z.shape[0] 
        if time:
            print >> fs, 'FIELD FieldData 1'
            print >> fs, 'TIME 1 1 double'
            print >> fs, time
            
        WriteGridPts(fs, x, 'X') 
        WriteGridPts(fs, y, 'Y') 
        WriteGridPts(fs, z, 'Z') 

        # Point data 
        names = data.keys() 
        noPts = data[names[0]].size 
        print >> fs, 'CELL_DATA', noPts 
        for name, vals in data.iteritems(): 
            WritePointData(fs, name, vals) 
            


# TETRAD Grid
grid = TetradGrid("base.gv")
centers, dim = grid.grid_spec()

#TETRAD GridView Data
def ProcessGridView():
    gv = TetradGridView("base.gv")
    gv.read_all_data('base_results')

ProcessGridView()
# Dimensions 

nx, ny, nz = dim['dx'].shape[0], dim['dy'].shape[0], dim['dz'].shape[0] 
lx, ly, lz = dim['dx'].sum(), dim['dy'].sum(), dim['dz'].sum()

ncells = nx * ny * nz 
npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Coordinates 
x = np.cumsum(dim['dx'])
x = np.insert(x,0,0)
y = np.cumsum(dim['dy'][::-1])   #inverted y axis in paraview vs tetrad
y = np.insert(y,0,0)
z = np.cumsum(dim['dz'])
z = np.insert(z,0,0)
z = z*-1


def ConvertDataTo3D(data_1d):
    data_3d = np.zeros((nx, ny, nz)) 
    
    for k in range(nz):
        for j in range(ny): 
            for i in range(nx):  
                index = i + nx*j + (nx*ny)*k
                data_3d[i,-j-1,k] = data_1d[index]
                
    return data_3d
    
# Variables 
#pressure = np.zeros((nx, ny, nz)) 
#temperature = np.zeros((nx, ny, nz)) 
block_3d = (np.zeros((nx, ny, nz)))

#df = data.loc[:,'TFR']
#data_array = np.array(df.tolist())
block_1d = np.array(centers.loc[:,'Block'].tolist())
block_3d = ConvertDataTo3D(block_1d)

        
#DoVTKLegacy(x, y, z, {'T': temperature, 'Block': block})

results_files = ["base_results_P.csv", "base_results_Sg.csv", 'base_results_T.csv']
param_names = ['P', 'SG', 'T']
isfile = 'base.is'
is_time, is_index = (0,0)

def TetradResultsToVTK(results_files, param_names, vtk_file_name, x, y, z, block_3d, skip = 1):
    
    data_dict = {'Block':block_3d}
    df = pd.read_csv(results_files[0])
    times = df.columns.tolist()[:-1][::skip]
    
    if isfile:
        isfo = TetradInterSim(isfile)
        is_df = isfo.read_data()
        is_params = is_df.columns.tolist()
        if 'Sg' in is_params:
            is_params[is_params.index('Sg')] = 'SG'
            is_df.columns = is_params
        for is_param in is_params:
            data_1d = np.array(is_df[is_param].tolist())
            data_dict[is_param] = ConvertDataTo3D(data_1d)
            
#        DoVTKLegacy(x, y, z, data_dict, vtk_file_name + '.vtk.' + is_index, time=time)
        
    
    for step, time in enumerate(times):
        print time
        for i, results_csv in enumerate(results_files):
            df = pd.DataFrame()
            df = pd.read_csv(results_csv)
            data_1d = np.array(df.loc[:,time])
            data_3d = ConvertDataTo3D(data_1d)
            data_dict[param_names[i]] = data_3d
            
        DoVTKLegacy(x, y, z, data_dict, vtk_file_name + '.vtk.' + str(step), time=time)
        
TetradResultsToVTK(results_files, param_names, 'base_results', x, y, z, block_3d)
