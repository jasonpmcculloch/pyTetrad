#Last Edited: 1/22/2016

import numpy as np
import re
import pandas as pd

class TetradFile(file):

    def read_lines_back(self, numlines=1):
        """Reads the line numlines(default is 1) before the current cursor position. Returns the given line and
        sets the cursor to the original position"""
        orig = self.tell()
        lines = 0
        c = ''
        self.seek(self.tell() - 1)
        while lines<=numlines:
            c = self.read(1)
            if c == "\n":
                lines += 1
                self.seek(self.tell() - 3)
            else: self.seek(self.tell() - 2)
        self.seek(self.tell() + 3)
        linestring = self.readline()
        self.seek(orig)
        return linestring

    def skipto(self,keyword='',start=0, getlinebefore=False):
        """Skips to line starting  with keyword.  keyword can be either a string or a list of strings, in which case
        it skips to a line starting with any of the specified strings.
        Returns the keyword found, or false if it can't find any of them.  The start parameter specifies which
        character in the line is to be considered the first to search."""
        line=''
        if isinstance(keyword,list): keywords=keyword
        else: keywords=[keyword]
        # Check if the passed object is a GridView or Intersim Object
        if isinstance(self, TetradGridView) or isinstance(self, TetradInterSim):
            while not any([line[start:].strip()==kw for kw in keywords]):
                line=self.readline()
                if line=='': return False
            return [kw for kw in keywords if line[start:].strip()==kw][0]
        else:
            while not any([line[start:].startswith(kw) for kw in keywords]):
                line=self.readline()
                if line=='': return False
            return [kw for kw in keywords if line[start:].startswith(kw)][0]

class TetradGridView(TetradFile):
    def __init__(self, filename):
        self.filename = filename
        super(TetradGridView,self).__init__(filename,'r')
        self.setup_pos()
        self.first()

    def setup_pos(self):
        """Sets up _fullpos dict and _pos list for TETRAD gridview files, containing file position of each set of results.
        Also sets up the times array."""
        self.seek(0)
        # set up pos,times
        t=[]
        exclude=["Dx", "Dy", "Top", "Dz", 'ActivBlks', 'Aquifer'] #exclude Top, Dx, Dy, Dz, 'ActivBlks', 'Aquifer'  in keywords
        keywords=[i for i in self.find_parameters() if i.lower() not in map(lambda x:x.lower(), exclude)]
        self._fullpos = {p.strip():[] for p in keywords}
        endfile=False
        while not endfile:
            kwfound=self.skipto(keywords)
            if kwfound:
                self._fullpos[kwfound.strip()]+=[self.tell()]
                self.read_time()
                t.append(self.time)
            else: endfile=True
        self.times=np.unique(np.array(t))
        self._fullpos = {k:v for k,v in self._fullpos.iteritems() if v <> []}
        first_parameter = min(self._fullpos, key=self._fullpos.get)
        self._pos = self._fullpos[first_parameter]
        self._params = self._fullpos.keys()

    def find_parameters(self):
        origin = self.tell()
        self.seek(0)
        self.readline()
        parameters = []
        endfile=False
        while not endfile:
            l = self.readline()
            if l[0].isalpha():
                parameters.append(l.strip())
            elif l.startswith("        -1.0"):
                endfile=True
        self.seek(origin)
        return parameters

    def read_time(self):
        self.time = float(self.read_lines_back(numlines=2).strip())

    def read_parameter(self):
        return self.read_lines_back(numlines=1).strip()

    def read_table(self):
        line = "  "
        values = []
        while (line.startswith("  ") or line.startswith(" -")) and not line.startswith("         0.0") and not line.startswith("        -1.0"):
            line = self.readline()
            values+=(line.strip().split())
        return (np.array(values[:-1]).astype(float))

    def read_all_data(self, filename = ""):
        """Read all the data vs time for each parameter in the GridView object.
        Writes a new GridVew object property called param_table
        
        Parameters:
            filename -- an filename (without the extension) to write the DataFrame on. 
        
        Object Property:
            self.param_table -- a dictionary with parameters as keys; containing DataFrames of data
        """
        self.param_table = {}
        for par in self._fullpos:
            df = pd.DataFrame()
            for pos in self._fullpos[par]:
                self.seek(pos)
                self.read_time()
                t = self.time
                values = self.read_table()
                df[t]=values
            self.param_table[par] = df
            
        if filename:
            for par in self.param_table:
                self.param_table[par].to_csv(filename + "_" + par + ".csv", index = False)

    def read_data(self):
        orig_pos = self.tell()
        self.read_time()
        pos_index = self._pos.index(self.tell())
        df = pd.DataFrame(columns = self._fullpos.keys())
        for param in df.columns:
            self.seek(self._fullpos[param][pos_index])
            df.loc[:,param] = self.read_table()
        self.seek(orig_pos)
        return df

    def first(self):
        self.index = 0
        self.seek(self._pos[self.index])
        self.read_time()

    def next(self):
        if self.index == len(self._pos)-1:
            print "EOF reached, no more results next to current time."
        else:
            self.index += 1
            self.seek(self._pos[self.index])
            self.read_time()

    def last(self):
        self.index = len(self._pos)-1
        self.seek(self._pos[self.index])
        self.read_time()

class TetradOut(TetradFile):
    def __init__(self, filename):
        self.filename = filename
        super(TetradOut,self).__init__(filename,'r')
        self.setup_pos()
        self.first()

    def setup_pos(self):
        """Sets up _pos and _recurpos list for TETRAD output files, containing file position of each set of results.
        Also sets up the times and recur times array in year units."""
        # set up pos,times
        t=[]
        tstep=[]
        keywords=["  TIME STEP"," RECURRENT INPUT"]
        self.seek(0)
        endfile=False
        self._pos=[]
        self._recurpos=[]
        while not endfile:
            kwfound=self.skipto(keywords)
            if kwfound=="  TIME STEP":
                self.read_time()
                if self.time_step in tstep:
                    self._pos.pop()
                    t.pop()
                    tstep.pop()
                self._pos+=[self.tell()]
                t.append(self.time)
                tstep.append(self.time_step)
            elif kwfound==" RECURRENT INPUT" and t<>[]:
                self._recurpos.append(self._pos[-1])
            elif kwfound==False:
                endfile=True
        self.times=t
        self.recur_indices=[self._pos.index(p) for p in self._recurpos]
        self.recur_times=[self.times[i] for i in self.recur_indices]

    def first(self):
        self.index = 0
        self.seek(self._recurpos[self.index])
        self.read_time()
        
    def next(self):
        if self.index == len(self._recurpos)-1:
            print "EOF reached, no more recurrent results next to current time."
        else:
            self.index += 1
            self.seek(self._recurpos[self.index])
            self.read_time()

    def last(self):
        self.index = len(self._recurpos)-1
        self.seek(self._recurpos[self.index])
        self.read_time()

    def read_time(self):
        """Reads the time of the current position in file. The default units for self.time is in years."""
        origin = self.tell()
        timestamp = self.readline().strip().split()
        self.time_step = int(timestamp[0])
        self.time_day = float(timestamp[3])
        self.time = float(timestamp[4])
        self.seek(origin)

    def read_well_table(self):
        """Reads the well table into a pandas DataFrame"""                
        origin = self.tell()
        l = ""
        endtable=False
        columns = ['BLOCK','LAYER','DRAWW',
                   'MASS FLOW STEAM','MASS FLOW WATER',
                   'ENERGY FLOW STEAM','ENERGY FLOW WATER',
                   'CUMULATIVE MASS STEAM','CUMULATIVE MASS WATER','CUMULATIVE MASS TOTAL',
                   'CUMULATIVE ENERGY STEAM','CUMULATIVE ENERGY WATER','CUMULATIVE ENERGY TOTAL',
                   'P','T']
        dtypes = [int]*2+[float]*13
        well_table_df = pd.DataFrame(columns = ['WELL']+columns)
        #read the headers
        self.skipto(" BLOCK  LAYER")
        self.readline()
        self.readline()

        while not endtable: #the whole table for the timestep
            well_table = []
            l=''
            while not l.startswith(' ****') and not l.startswith('\n'): #table per well with the totals
                l=self.readline()
                if l!='  \n':
#                    well_table.append(l.strip().split())
                    curr_line = []
                    curr_line.append(l[0:6])
                    curr_line.append(l[7:13])
                    curr_line.append(l[14:20])
                    curr_line = curr_line + [l[20:][i:i+9] for i in range(0,(len(columns)-3)*9,9)]
                    well_table.append(curr_line)
            if l=='\n': endtable=True      
            if not endtable:
                well_table.pop()
                wellname = well_table[-1][0]
                well_df = pd.DataFrame(well_table[:-1],columns=columns)
                well_df.insert(0,'WELL',wellname)
                well_table_df=well_table_df.append(well_df, ignore_index=True)

        #apply datatypes
        for i,col in enumerate(columns):
            well_table_df.loc[:,col]=well_table_df.loc[:,col].astype(dtypes[i])
        well_table_df.set_index('WELL', inplace=True)
        self.seek(origin)
        
        return well_table_df


class TetradInterSim(TetradFile):
    def __init__(self, filename):
        self.filename = filename
        super(TetradInterSim,self).__init__(filename,'r')
        self.setup_pos()
        
    def setup_pos(self):
        """Sets up _fullpos dict and _pos list for TETRAD intersim files, containing file position of each set of results."""
        self.seek(0)
        # set up pos
        exclude = ["DX", "DY", "Top", "Trans Mods", "Wells"]
        keywords=[i for i in self.find_parameters() if i.lower() not in map(lambda x:x.lower(), exclude)]
        self._fullpos = {p.strip():None for p in keywords}
        endfile=False
        while not endfile:
            kwfound=self.skipto(keywords)
            if kwfound:
                self._fullpos[kwfound.strip()]=self.tell()
            else: endfile=True
        first_parameter = min(self._fullpos, key=self._fullpos.get)
        self._params = self._fullpos.keys()
    
    def find_parameters(self):
        origin = self.tell()
        self.seek(0)
        self.readline()
        parameters = []
        endfile=False
        while not endfile:
            l = self.readline()
            if l:
                if l[0].isalpha():
                    parameters.append(l.strip())
            else:
                endfile=True
        self.seek(origin)
        return parameters
    
    def read_table(self):
        line = "  "
        values = []
        while line.startswith("  ") and not line[0].isalpha():
            line = self.readline()
            if line <> '':
                if not line[0].isalpha():
                    values+=(line.strip().split())
        return (np.array(values).astype(float))
        
    def read_data(self):
        orig_pos = self.tell()
        df = pd.DataFrame(columns = self._fullpos.keys())
        for param in df.columns:
            self.seek(self._fullpos[param])
            df.loc[:,param] = self.read_table()
        self.seek(orig_pos)
        return df
        
class TetradGrid(TetradFile):
    def __init__(self, filename):
        self.filename = filename
        super(TetradGrid,self).__init__(filename,'r')
        self.check_filetype()
        
    def check_filetype(self):
        """Checks filetype and defines read_table method using getattr and setattr"""
        self.seek(0)
        filetype = self.readline().split()[0]
        if filetype == "INTERSIM":
            setattr(self,'read_table',getattr(self, "read_table_intersim"))
        elif filetype == "GRIDVIEW":
            setattr(self,'read_table',getattr(self, "read_table_gridview"))
        else:
            print "Cannot recognize filetype"
        self.filetype = filetype
    
    def read_table_intersim(self):
        line = "  "
        values = []
        while line.startswith("  ") and not line[0].isalpha():
            line = self.readline()
            if not line[0].isalpha():
                values+=(line.strip().split())
        return (np.array(values).astype(float))
    
    def read_table_gridview(self):
        line = "  "
        values = []
        while line.startswith("  ") and not line.startswith("         0.0") and not line.startswith("        -1.0"):
            line = self.readline()
            values+=(line.strip().split())
        return (np.array(values[:-1]).astype(float))
        
    def grid_spec(self):
        """Get the block centers (x, y, z) of each block in a TETRAD grid.
        
        Returns:
            grid_df -- a dataframe containing the block centers (x, y, z) 
                       of each block
            grid_dim -- a dictionary containing lists of the grid dimensions dx, dy and dz
        """
        
        self.seek(0)
        self.readline()
        
        if self.filetype=="INTERSIM":
            self.readline()
            nx = np.array(self.readline(), dtype=int)
            ny = np.array(self.readline(), dtype=int)
            nz = np.array(self.readline(), dtype=int)
            
        elif self.filetype=="GRIDVIEW":
            nx, ny, nz = np.array(self.readline().strip().split(), dtype=int)
            
        self.seek(0)
        self.skipto(["DX","Dx"])
        dx = self.read_table()
        if self.filetype=="INTERSIM":
            dx = dx[:nx]
        x_centers = self.block_centers(dx)
        x_centers = np.tile(x_centers, ny*nz)
        
        self.seek(0)
        self.skipto(["DY","Dy"])
        dy = self.read_table()
        if self.filetype=="INTERSIM":
            dy = dy.reshape(ny,nx).T[0]
        y_centers = self.block_centers(dy)
        y_centers = np.tile(y_centers.reshape(ny,1),nx)
        y_centers = np.reshape(y_centers, nx*ny)
        y_centers = np.tile(y_centers, nz)
        
        self.seek(0)
        self.skipto(["Thickness", "DZ", "Dz"])
        dz = self.read_table()
        dz = dz.reshape(nz, nx*ny).T[0]
        z_centers = self.block_centers(dz)
        z_centers = np.tile(z_centers.reshape(nz,1),nx*ny)
        z_centers = np.reshape(z_centers, nx*ny*nz)
        
        grid_df = pd.DataFrame(np.array((x_centers, y_centers, z_centers)).T, columns = ["X", "Y", "Z"])
        grid_df.loc[:,"Block"] = pd.Series(range(1,nx*ny*nz+1))
        
        grid_dim = {'dx':dx, 'dy':dy, 'dz':dz}     
        
        return grid_df, grid_dim
                   
    def connect_refined_area(self, ref_grid, x_ref_index, y_ref_index, z_ref_index):
        """Connect a refined area grid dataframe to the main grid dataframe.
        
        Input argument(s):
            ref_grid -- the refined grid TetradGrid object.
            ref_index -- a tuple (ix, iy, iz) that indicates the start indices of
                         refinement
        
        Returns:
            grid_combined -- a combined grid and ref_grid dataframe with corrected
                             block centers
    
        """
        grid_df, grid_dim = self.grid_spec()
        ref_grid_df, ref_grid_dim = ref_grid.grid_spec()
        
        x_offset = sum(grid_dim['dx'][:x_ref_index - 1])
        y_offset = sum(grid_dim['dy'][:y_ref_index - 1])
        z_offset = sum(grid_dim['dz'][:z_ref_index - 1])
        block_offset = grid_df.shape[0] * 2 #get the fracture blocks only
        
        ref_grid_df.loc[:,'X'] += x_offset
        ref_grid_df.loc[:,'Y'] += y_offset
        ref_grid_df.loc[:,'Z'] += z_offset
        ref_grid_df.loc[:,'Block'] += block_offset
        
        grid_combined = grid_df.append(ref_grid_df)
        
        return grid_combined
        
    
    def block_centers(self, dx):
        n=len(dx)
        nodes = [0]
        for x in dx:
            nodes.append(x+nodes[-1])
        centers = [(nodes[i+1]+nodes[i])/2 for i in xrange(n)]
        return np.array(centers)
            

          
            