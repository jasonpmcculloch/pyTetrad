# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
grid = TetradGrid("T2-R22.GV")
centers, dim = grid.grid_spec()
data = pd.read_csv('Grid.csv')

#TETRAD GridView Data
def ProcessGridView():
    gv = TetradGridView("T2-R22.GV")
    gv.read_all_data('T2-R22_results')

# Dimensions 

nx, ny, nz = dim['dx'].shape[0], dim['dy'].shape[0], dim['dz'].shape[0] 
lx, ly, lz = dim['dx'].sum(), dim['dy'].sum(), dim['dz'].sum()

ncells = nx * ny * nz 
npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Coordinates 
x = np.cumsum(dim['dx'])
x = np.insert(x,0,0)
y = np.cumsum(dim['dy'])
y = np.insert(y,0,0)
z = np.cumsum(dim['dz'])
z = np.insert(z,0,0)
z = z*-1

# Variables 

#pressure = np.zeros((nx, ny, nz)) 
#temperature = np.zeros((nx, ny, nz)) 
block_3d = (np.zeros((nx, ny, nz)))

#df = data.loc[:,'TFR']
#data_array = np.array(df.tolist())
block_1d = np.array(data.loc[:,'Block'].tolist())

for k in range(nz):
    for j in range(ny): 
        for i in range(nx): 
            index = i + nx*j + (nx*ny)*k
#            print i,j,k,index
#            temperature[i,-j,k] = data_array[index]
            block_3d[i,-j,k] = block_1d[index]

        
#DoVTKLegacy(x, y, z, {'T': temperature, 'Block': block})

results_files = ["T2-R22_results_PFR.csv", "T2-R22_results_SGFR.csv", 'T2-R22_results_TFR.csv']
param_names = ['PFR', 'SGFR', 'TFR']

def TetradResultsToVTK(results_files, param_names, vtk_file_name, x, y, z, block_3d, skip = 1):
    
    data_dict = {'Block':block_3d}
    
    df = pd.read_csv(results_files[0])
    times = df.columns.tolist()[:-1][::skip]
    
        
    
    for step, time in enumerate(times):
        print time
        for i, results_csv in enumerate(results_files):
            df = pd.DataFrame()
            df = pd.read_csv(results_csv)
            data_1d = np.array(df.loc[:,time])
            data_3d = ConvertDataTo3D(data_1d)
            data_dict[param_names[i]] = data_3d
            
        DoVTKLegacy(x, y, z, data_dict, vtk_file_name + '.vtk.' + str(step), time=time)
        
        
def ConvertDataTo3D(data_1d):
    data_3d = np.zeros((nx, ny, nz)) 
    
    for k in range(nz):
        for j in range(ny): 
            for i in range(nx):  
                index = i + nx*j + (nx*ny)*k
                data_3d[i,-j,k] = data_1d[index]
                
    return data_3d
                
TetradResultsToVTK(results_files, param_names, 'T2-R22_results', x, y, z, block_3d, 30)