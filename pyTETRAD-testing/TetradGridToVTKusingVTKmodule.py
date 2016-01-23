# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 11:22:05 2016

@author: alcaraz.jt
"""

from vtk import vtkRectilinearGrid, vtkDoubleArray, vtkRectilinearGridWriter, vtkStringArray, vtkXMLRectilinearGridWriter
from pyTETRAD import TetradGrid, TetradGridView
import numpy as np

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


# TETRAD Grid
grid = TetradGrid("base.gv")
centers, dim = grid.grid_spec()

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
z = z*-1
z = np.insert(z,0,0)


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

#Set up the vtk grid
grid = SetVtkGrid(x,y,z)

def AddScalarData(grid, param_name, tetrad_data, nx, ny, nz):
    """Add scalar/string data to a vtk grid.
    
    Parameters:
        grid -- vtk grid where the scalar data is added 
        param_name -- the name of the scalar dataset
        tetrad_data -- the data extracted from an IS or GV file
        nx, ny, nz -- the number of cells in the x, y and z directions respectively
    """
    data_3d = (np.zeros((nx, ny, nz)))
    data_3d = ConvertDataTo3D(tetrad_data, nx, ny, nz)
    if all(isinstance(d,str) for d in tetrad_data):
        dataArray = vtkStringArray()
    else:
        dataArray = vtkDoubleArray()
    dataArray.SetName(param_name)
    for d in data_3d.flatten('F'): dataArray.InsertNextValue(d)
    grid.GetCellData().AddArray(dataArray)

#insert scalar data
block_1d = np.array(range(1, grid.GetNumberOfCells()+1))
AddScalarData(grid, 'Block', block_1d, nx, ny, nz)

#insert string data
block_names = ['This is block {}.'.format(str(b)) for b in block_1d]
AddScalarData(grid, 'Block_Name', block_names, nx, ny, nz)



def WriteVtkFile(grid, filename, time=None):
    """Writes a vtk Grid object to a file.
    
    Parameters:
        grid -- vtk Grid Object
        filename -- the filename for the vtk file
    """
    
    #insert time data; this has a bug on the vtk writer, manually insert time field instead.
    #time = vtkDoubleArray()
    #time.SetName('Time')
    #time.InsertNextValue(np.array(1999))
    #grid.GetFieldData().AddArray(time)
    
    #write to .vtk file
    writer = vtkRectilinearGridWriter();
    writer.SetFileName(filename)
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
        vtk_file.insert(line_index-1,'FIELD FieldData 1\nTime 1 1 double\n' + str(1999))
        with open(filename, 'w+') as fo:
            fo.writelines('\n'.join(vtk_file))

    #write to .vtr file; there is also a bug for VTK 6.2.0, update needed to VTK 6.3.0
    #writer = vtkXMLRectilinearGridWriter();
    #writer.SetFileName('Sample Rectilinear.vtr')
    #writer.SetInputData(grid)
    #writer.Write()

#write the vtk file
WriteVtkFile(grid, 'Sample Rectilinear.vtk', time = 1999)

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
    
#process GridView
ProcessGridView('base.gv', 'base')