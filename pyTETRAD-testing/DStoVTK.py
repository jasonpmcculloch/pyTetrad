#-------------------------------------------------------------------------------
# Name:        
# Purpose:     1. Converting NIGBU actual coordinates to model coordinates
#              2. Finding the most probable block given a set of fault intersections
# Author:      alcaraz.jt
#
# Created:     /6/2015
#
# Copyright:   (c) alcaraz.jyt 2015
#-------------------------------------------------------------------------------

from vtk import vtkUnstructuredGrid, vtkDoubleArray, vtkUnstructuredGridWriter, vtkXMLUnstructuredGridWriter, vtkStringArray, \
    vtkPoints, vtkCellArray, vtkPolyLine, VTK_POLY_LINE
from numpy import cos, sin, sqrt
from math import radians
import pandas as pd
from pdb import set_trace
#from mulgrids import *
#from t2data import *
from collections import Counter


def main():
    welltracks_df = convertDeviationSurvey('DS_updated.xlsx')
    makeVTKWellsUsingModule('Wells', welltracks_df, xml = True)

def mainOld():
    global geo, blockmap
    xldata_path = r"d:\Users\alcaraz.jt\Desktop\NIGBU Reservoir Modeling\NIGBU drilling sequence_2013-2031.xlsx"
    xldata_sheet = "Summary Table"
    faultdata_df = pd.read_excel(xldata_path,xldata_sheet, header=1)

    faultdata_df['model_x_from'], faultdata_df['model_y_from'] = transformCoords((faultdata_df.loc[:,'mE'], faultdata_df.loc[:,'mN']))
    faultdata_df['model_x_to'], faultdata_df['model_y_to'] = transformCoords((faultdata_df.loc[:,'mE_to'], faultdata_df.loc[:,'mN_to']))

    fault_intersections = [Segment(Point((faultdata_df.model_x_from[i], faultdata_df.model_y_from[i], faultdata_df.mRSL[i])), 
        Point((faultdata_df.model_x_to[i], faultdata_df.model_y_to[i], faultdata_df.mRSL_to[i]))) for i in faultdata_df.index]

    faultdata_df['Segment_Length'] = [i.length for i in fault_intersections]

    t2_datfile = t2data('NSN07.dat')
    geo, blockmap = t2_datfile.grid.rectgeo() #mulgrid file without wells

    blocks_intersected_list = [findIntersectedBlock(fault_intersections[i]) for i in faultdata_df.index[:2]]

    for i in range(3):
        faultdata_df['block_' + str(i)], faultdata_df['count_' + str(i)] = np.array(blocks_intersected_list)[:,i].T

    faultdata_df.to_excel("Fault Intersection Data.xlsx")
    set_trace()

class Point:
    def __init__(self, (x, y, z)):
        self.x = x
        self.y = y
        self.z = z
        self.coords = (x, y, z)

class Segment:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.length = sqrt((p2.x-p1.x)**2 + (p2.y-p1.y)**2 + (p2.z-p1.z)**2)

    def divide(self, num_segments):
        x_points = np.linspace(self.p1.x, self.p2.x, num_segments)
        y_points = np.linspace(self.p1.y, self.p2.y, num_segments)
        z_points = np.linspace(self.p1.z, self.p2.z, num_segments)
        line_points = [Point(i) for i in zip(x_points,y_points,z_points)]
        return line_points 

def findIntersectedBlock(line, divisions = 100, ranks = 3): #geo and blockmap are mulgraph classes
    points = line.divide(divisions)
    blocks_intersected = [geo.block_name_containing_point([p.x,p.y,p.z]) for p in points]
    count = Counter(blocks_intersected)
    most_probable_blocks = count.most_common(ranks)
    return most_probable_blocks
    
    
def transformCoords(x, y):

    origin_coords = (495658.8468, 1024275.7711) #top left corner coordinates in model
#    origin_coords = (521460.0927, 1003382.334) #bottom right corner coordinates in model
#    origin_coords = (542856.9860, 1029805.2969)
#    origin_coords = (517055.7401, 1050698.7338)
    rotation_rad = radians(51) #degrees of rotation from actual to model
    x = x - origin_coords[0]
    y = y - origin_coords[1]
    return (x*cos(rotation_rad) + y*sin(rotation_rad), 
         - x*sin(rotation_rad) + y*cos(rotation_rad))
         
def convertDeviationSurvey(fname):
    welltracks = pd.read_excel(fname)
    ly = 33200  #Tetrad Grid length in Y direction
    x, y = transformCoords(welltracks['Easting'], welltracks['Northing'])
    welltracks['X'] = x
    welltracks['Y'] = y + ly  #case of tetrad
    return welltracks
    
def makeVTKWells(fname, welltracks_df):
    
    numpoints = welltracks_df.shape[0]
    wells = welltracks_df['Well'].unique().tolist()
    numwells = len(wells)


    with open(fname, 'w') as fo:
        print >> fo, '# vtk DataFile Version 1.0'
        print >> fo, 'Well Tracks'
        print >> fo, 'ASCII\n'
        print >> fo, 'DATASET UNSTRUCTURED_GRID'
        
        print >> fo, 'POINTS {} float'.format(numpoints)
        
        for i in range(numpoints):
            print >> fo, welltracks_df.loc[i,'X'], welltracks_df.loc[i,'Y'], welltracks_df.loc[i,'Elev_mASL']
            
        print >> fo, '\nCELLS {0} {1}'.format(numwells, numwells + numpoints)
        
        for well in wells:
            indices = welltracks_df[welltracks_df['Well']==well].index.tolist()
            print >> fo, len(indices),
            for i in indices:
                print >> fo, i,
            print >> fo
            
        print >> fo, '\nCELL_TYPES {0}'.format(numwells)
        
        for i in range(numwells):
            print >> fo, '4'

def makeVTKWellsUsingModule(fname_base, welltracks_df, xml=False):
    
    numpoints = welltracks_df.shape[0]
    wells = welltracks_df['Well'].unique().tolist()
    numwells = len(wells)

    grid = vtkUnstructuredGrid()
    points = vtkPoints()  
    
    for i in range(numpoints):
        points.InsertNextPoint(welltracks_df.loc[i,'X'], welltracks_df.loc[i,'Y'], welltracks_df.loc[i,'Elev_mASL'])
    
    cells = vtkCellArray()
    wellname = vtkStringArray()
    wellname.SetName('Well')
    
    for well in wells:
        print well
        polyline = vtkPolyLine()
        indices = welltracks_df[welltracks_df['Well']==well].index.tolist()
        for i, j in enumerate(indices):
            polyline.GetPointIds().SetNumberOfIds(len(indices))
            polyline.GetPointIds().SetId(i,j)
            
        cells.InsertNextCell(polyline)
        wellname.InsertNextValue(well)
        
    grid.SetPoints(points)
    grid.SetCells(VTK_POLY_LINE, cells)
    grid.GetCellData().AddArray(wellname)
    
    if xml:
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName('{}.vtu'.format(fname_base))
        writer.SetDataModeToAscii()
        writer.SetInputData(grid)
        writer.Write()
        
    else:
        writer = vtkUnstructuredGridWriter()
        writer.SetFileName('{}.vtk'.format(fname_base))
        writer.SetInputData(grid)
        writer.Write()
    
if __name__ == '__main__':
    main()