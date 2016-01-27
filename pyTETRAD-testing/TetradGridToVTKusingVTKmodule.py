# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:22:05 2016

@author: alcaraz.jt
"""

from vtk import vtkRectilinearGrid, vtkDoubleArray, vtkRectilinearGridWriter, vtkStringArray, vtkXMLRectilinearGridWriter
from pyTETRAD import TetradGrid, TetradGridView, TetradInterSim
import numpy as np
import pandas as pd

def ConvertDataTo3D(data_1d, nx, ny, nz):
    data_3d = np.zeros((nx, ny, nz)) 
    if all(isinstance(d,str) for d in data_1d):
        data_3d=data_3d.astype(str)
    for k in range(nz):
        for j in range(ny): 
            for i in range(nx):  
                index = i + nx*j + (nx*ny)*k
                data_3d[i,-j-1,k] = data_1d[index]
                
    return data_3d

def SetVtkGrid(x,y,z):
    """Set up the vtk rectilinear grid using x, y, z data.
    
    Parameters:
        x, y, z -- the points in the x, y and z directions respectively
    Returns:
        grid -- a vtkRectilinearGrid object
    """
    grid = vtkRectilinearGrid();
    grid.SetDimensions(len(x),len(y),len(z));
     
    xArray = vtkDoubleArray();
    for xCoord in x: xArray.InsertNextValue(xCoord)
     
    yArray = vtkDoubleArray();
    for yCoord in y: yArray.InsertNextValue(yCoord)
    
    zArray = vtkDoubleArray();
    for zCoord in z: zArray.InsertNextValue(zCoord)
     
    grid.SetXCoordinates(xArray);
    grid.SetYCoordinates(yArray);
    grid.SetZCoordinates(zArray);
    
    print "There are " + str(grid.GetNumberOfPoints()) + " points.";
    print "There are " + str(grid.GetNumberOfCells()) + " cells.";

    return grid

def AddScalarData(grid, param_name, tetrad_data):
    """Add scalar/string data to a vtk grid.
    
    Parameters:
        grid -- vtk grid where the scalar data is added 
        param_name -- the name of the scalar dataset
        tetrad_data -- the data extracted from an IS or GV file
        nx, ny, nz -- the number of cells in the x, y and z directions respectively
    """
    nx, ny, nz = [i-1 for i in grid.GetDimensions()]    
    
    data_3d = (np.zeros((nx, ny, nz)))
    data_3d = ConvertDataTo3D(tetrad_data, nx, ny, nz)
    if all(isinstance(d,str) for d in tetrad_data):
        dataArray = vtkStringArray()
    else:
        dataArray = vtkDoubleArray()
    dataArray.SetName(param_name)
    for d in data_3d.flatten('F'): dataArray.InsertNextValue(d)
    grid.GetCellData().AddArray(dataArray)

def WriteVtkFile(grid, filename, time=None, xml=False, binary=False):
    """Writes a vtk Grid object to a file.
    
    Parameters:
        grid -- vtk Grid Object
        filename -- the filename for the vtk file
        xml -- write the vtk in xml format
        binary -- writes in binary format; xml should be True
    """
    if binary: xml=True
    
    if not xml:
        #write to .vtk file
        writer = vtkRectilinearGridWriter();
        writer.SetFileName(filename)
        if binary:
            print 
        writer.SetInputData(grid)
        writer.Write()
        
        if time:
            #manually add the time data
            vtk_file = open(filename, 'r+').read().split('\n')
            line_index = 0
            line = ''
            _pos = []
            while 'CELL_DATA' not in line:
                line = vtk_file[line_index]
                line_index += 1
            vtk_file.insert(line_index-1,'FIELD FieldData 1\nTime 1 1 double\n' + str(time))
            with open(filename, 'w+') as fo:
                fo.writelines('\n'.join(vtk_file))
    
    else:
        if time:
            #insert time data; this has a bug on the vtk legacy writer, manually insert time field instead for vtk legacy files.
            timeArray = vtkDoubleArray()
            timeArray.SetName('Time')
#            print time, type(time)
            timeArray.InsertNextValue(time)
            grid.GetFieldData().AddArray(timeArray)
    
        writer = vtkXMLRectilinearGridWriter();   
        if binary:
            writer.SetDataModeToBinary()
        else: writer.SetDataModeToAscii()
        print filename
        writer.SetFileName(filename)
        writer.SetInputData(grid)
        writer.Write()

#process TETRAD GridView Data
def ProcessGridView(gv_filename, csv_filename_base):
    """Reads a GridView file and extracts the parameter data to csv files. 
    
    Parameters:
        gv_filename -- GridView filename to read.
        csv_filename_base -- the filename base for the csv files. 
                             e.g. myFile -> becomes myFile_param1.csv, ..., etc.
    """
    gv = TetradGridView(gv_filename)
    gv.read_all_data(csv_filename_base)


def TetradResultsToVtk(grid, results_files, param_names, vtk_filename_base, skip = 1, isfile = None, xml = False, binary = False):
    """Writes the GridView results from the .csv files to the vtkGrid object and writes the vtk files.
    Also reads InterSim files if available.
    
    Parameters:
        grid -- vtk Grid Object
        results_files -- a list containing the filesnames of the .csv files to read from
        param_names -- the equivalent parameter names for each of the results_files
        skip -- time steps to skip, default = 1 (means include everything)
        isfile -- the intersim file to read, default = None
        xml -- writes the vtk files in .vtr XML format. This also writes a .pvd file.
        binary -- writes the .vtk or .vtr files in binary format; xml should be True when setting bin to True.
    """    
    nx, ny, nz = [i-1 for i in grid.GetDimensions()]        
    
    #insert block number data
    block_1d = np.array(range(1, grid.GetNumberOfCells()+1))
    AddScalarData(grid, 'Block', block_1d)
    
    df = pd.read_csv(results_files[0])
    times = df.columns.tolist()[:-1][::skip]
    vtk_filenames = []    
    
    if binary: xml=True    
    
    if xml:
        file_extension = '.vtr'
    else: file_extension = '.vtk'    
    
    offset = 0
    #read the InterSim file
    if isfile:
        offset = 1
        isfo = TetradInterSim(isfile)
        is_df = isfo.read_data()
        is_params = is_df.columns.tolist()
        #Note: For IS files, the parameter names are in 'Sentence case'. GV parameters are 'ALLCAPS'. 
#        if 'Sg' in is_params:
#            is_params[is_params.index('Sg')] = 'SG'
#            is_df.columns = is_params
        for is_param in is_params:
            data_1d = np.array(is_df[is_param].tolist())
            AddScalarData(grid, is_param, data_1d)
    
    results_dict={}
    for i, results_csv in enumerate(results_files):
        results_dict[param_names[i]] = pd.read_csv(results_csv)
        
    for step, time in enumerate(times):
        print time
        for param, df in results_dict.iteritems():
            data_1d = np.array(df.loc[:,time])
            AddScalarData(grid, param, data_1d)
        
        vtk_filenames += [vtk_filename_base + '_' +str(step+offset) + file_extension]
        WriteVtkFile(grid, vtk_filenames[-1], time=float(time), xml=xml, binary=binary)
    
    if xml:
        write_pvd_file(vtk_filename_base + '.pvd', times, vtk_filenames)
        
def write_pvd_file(pvd_filename, times, files):
    """Writes a .pvd file relating the .vtr files
    
    Parameters:
        pvd_filename -- name of the .pvd file
        times -- list of times
        files -- list of .vtr filenames
    """ 
    with open(pvd_filename, 'w+') as fo:
        print >> fo, '<?xml version="1.0"?>'
        print >> fo, '<VTKFile type="Collection" version="0.1">'
        print >> fo, '  <Collection>'
        for i, time in enumerate(times):
            print >> fo, '    <DataSet timestep="{0}" file="{1}"/>'.format(time, files[i])
        print >> fo, '  </Collection>'
        print >> fo, '</VTKFile>'
            

# TETRAD Grid
grid = TetradGrid('T2-R22.GV')
centers, dim = grid.grid_spec()

# Dimensions 
nx, ny, nz = dim['dx'].shape[0], dim['dy'].shape[0], dim['dz'].shape[0] 
lx, ly, lz = dim['dx'].sum(), dim['dy'].sum(), dim['dz'].sum()

ncells = nx * ny * nz 
npoints = (nx + 1) * (ny + 1) * (nz + 1) 

# Coordinates 
x = np.cumsum(dim['dx'])
x = np.insert(x,0,0)
y = np.cumsum(dim['dy'][::-1])
y = np.insert(y,0,0)
z = np.cumsum(dim['dz'])
z = z*-1
z = np.insert(z,0,0)

#Set up the vtk grid
grid = SetVtkGrid(x,y,z)

#insert scalar data
block_1d = np.array(range(1, grid.GetNumberOfCells()+1))
AddScalarData(grid, 'Block', block_1d)

#insert string data; this has a bug in VTK 6.2.0, should update to 6.3.0 if available
#block_names = ['This is block {}.'.format(str(b)) for b in block_1d]
#AddScalarData(grid, 'Block_Name', block_names)

results_files = ["T2-R22_results_PFR.csv", "T2-R22_results_SGFR.csv", 'T2-R22_results_TFR.csv']
param_names = ['PFR', 'SGFR', 'TFR']
isfile = 'T2-R22A.IS'

#process GridView
#ProcessGridView('base.gv', 'base')

#write the vtk files from GridView and InterSim
TetradResultsToVtk(grid, results_files, param_names, 'T2-R22_results', skip=20, isfile=isfile, xml=True, binary=True)
