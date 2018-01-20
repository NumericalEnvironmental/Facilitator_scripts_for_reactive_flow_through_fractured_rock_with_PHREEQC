############################################################################
#
# FracProcess.py
#
# python script to write PHREEQC dual-porosity (fracture) flow simulation
# results to 3-D grid for plotting
#
# by Walt McNab
#
############################################################################

from __future__ import print_function
from numpy import *
from pandas import *

refVec = array([1., 0., 0.])        # reference unrotated normal vector

            
### classes ###

              
class Grid:

    def __init__(self):

        # read grid parameters
        lineInput = []        
        inputFile = open('grid_setup.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.x0 = float(lineInput[1][1])                # grid origin
        self.y0 = float(lineInput[1][2])
        self.z0 = float(lineInput[1][3])
        xLength = float(lineInput[2][1])           # length of sides
        yLength = float(lineInput[2][2])
        zLength = float(lineInput[2][3]) 
        self.nx = int(lineInput[3][1])                  # discretization
        self.ny = int(lineInput[3][2])
        self.nz = int(lineInput[3][3])
        self.dx = xLength/self.nx                  # cell size
        self.dy = yLength/self.ny
        self.dz = zLength/self.nz       
        self.xf = self.x0 + xLength             # end of grid along each axis
        self.yf = self.y0 + yLength         
        self.zf = self.z0 + zLength         
        print('Read grid parameters.')

        # read background parameters list and associated quantities
        self.bkgName = []
        self.bkgNum = []
        inputFile = open('background.txt','r')
        for line in inputFile:
            self.bkgName.append(line.split()[0])
            self.bkgNum.append(float(line.split()[1]))            
        inputFile.close()

    def FillGrid(self, fracs):
        # gross fill-in of background grid points
        indxSets = mgrid[0:self.nx, 0:self.ny, 0:self.nz].T.reshape(-1, 3)
        indxCols = transpose(indxSets)
        x = self.x0 + (indxCols[0] + 0.5) * self.dx
        y = self.y0 + (indxCols[1] + 0.5) * self.dy                
        z = self.z0 + (indxCols[2] + 0.5) * self.dz
        df = DataFrame({'x': x, 'y': y, 'z': z})
        # remove points within fracture zones
        df['fracZone'] = 0
        for frac in fracs: df['fracZone'] += frac.WithinFrac(x, y, z)
        df = df[df['fracZone']==0]        
        df.drop(['fracZone'], axis=1, inplace=True)
        # add background quantities
        for i, item in enumerate(self.bkgName):
            df[item] = self.bkgNum[i]
        return df



class Fracture:

    def __init__(self, fileName, nz, nx, ny, dz, dx, dy, w, phi, theta, datum, y0, grid_x_intercept, nzDivide):
        # fracture properties (and simulation results, as data frames)
        self.name = fileName.split('.')[0]                          # fracture name
        self.nz = nz                                                # discretization along fracture
        self.nx = nx                                                # discretization transverse to fracture
        self.ny = ny                                                # extruded discertization, along y-direction
        self.dz = dz
        self.dx = dx
        self.dy = dy                                                # fracture depth (for gridding, independent of PHREEQC)
        self.w = w                                                  # width of advective zone
        self.phi = phi*pi/180.                                      # dip below horizontal (input as degrees)
        self.theta = theta*pi/180.                                  # strike, clockwise from N (input as degrees)
        self.datum = datum                                          # z value at coordinate system rotation point  
        self.y0 = y0                                                # y-axis extrude starting location
        self.grid_x_intercept = grid_x_intercept                    # grid x-axis fracture crossing location
        self.nzDivide = nzDivide                                    # number of gridding subdivisions along advection axis
        self.phreeqc = read_csv(fileName, delim_whitespace=True)    # read associated PHREEQC transport output file
        self.fracVec = self.NormPlane()                             # normal vector to fracture plane
        self.rotMatrix = self.RotationMatrix()  
        self.eqn = self.EqnPlane()                                  # equation of fracture plane

    def WithinFrac(self, x, y, z):
        # determine if point (x, y, z) is within the fracture zone
        interior = (self.DistPtPlane(x, y, z) <= (self.w + self.nx*self.dx))
        return interior

    def DistPtPlane(self, x, y, z):
        # normal distance between point (x, y, z) and fracture plane
        # see College Calculus with Analytic Geometry, 3rd Edition, Protter and Morrey (1982); p. 410
        num = abs(self.eqn[0]*x + self.eqn[1]*y + self.eqn[2]*z + self.eqn[3])
        den = sqrt(self.eqn[0]**2 + self.eqn[1]**2 + self.eqn[2]**2)
        return num/den

    def EqnPlane(self):
        # determine equation of fracture plane in the form Ax + By + Cz + D = 0
        # see http://mathonline.wikidot.com/point-normal-form-of-a-plane
        A = self.fracVec[0]
        B = self.fracVec[1]
        C = self.fracVec[2]
        D = -A*self.grid_x_intercept
        return array([A, B, C, D])

    def NormPlane(self):
        # calculate the vector normal to a plane of strike theta and dip phi in 3-space
        # see http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML1%20FaultNormals.pdf
        nx = sin(self.phi) * cos(self.theta)    
        ny = -sin(self.phi) * sin(self.theta)
        nz = -cos(self.phi)
        return array([nx, ny, nz])
    
    def RotationMatrix(self):
        # methodology: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        v = cross(refVec, self.fracVec)
        c = dot(refVec, self.fracVec)
        vx = array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
        R = identity(3) + vx + linalg.matrix_power(vx, 2) * 1./(1. + c)
        return R      

    def ProcessSimResults(self, groupName, phaseList, initPhases, molesZero, bkgName):
        # clean up PHREEQC results dataframe and process coordinates

        # (1) limit results to last step
        endStep = self.phreeqc['step'].max()                        
        self.phreeqc = self.phreeqc[self.phreeqc['step']==endStep]

        # (2) remove columns pertaining to changing mass
        results = list(self.phreeqc.columns)                
        for i in xrange(len(results)-1, -1, -1):
            if (results[i][0]=='d') and ((results[i][1]=='_') or (results[i][2]=='_')):
                del results[i]
        self.phreeqc = self.phreeqc[results]

        # (3) group together phases/kinetics items, as user-specified
        for i, group in enumerate(groupName):
            self.phreeqc[group] = 0.
            for phase in phaseList[i]:
                self.phreeqc[group] += self.phreeqc[phase]

        # (4) mass losses for user-selected phases.kinetics items
        for i, phase in enumerate(initPhases):
            name = 'delta_' + phase
            self.phreeqc[name] = -(molesZero[i]-self.phreeqc[phase])*sign(self.phreeqc[phase])

        # (5) work out cell coordinates
        self.phreeqc['z'] = self.datum - self.phreeqc['dist_x']
        self.phreeqc['col'] = ((self.phreeqc['soln']-1.5)/self.nz).astype('int') * (self.phreeqc['soln']>self.nz)
        self.phreeqc['x'] = 0.5*self.w*(self.phreeqc['col']==0) + \
            (self.w + (self.phreeqc['col']-0.5)*self.dx)  * (self.phreeqc['col']!=0)
        
        # (6) interpolate at intermediate grid points along (unrotated) z-axis
        for i in xrange(self.nzDivide):
            between = self.phreeqc.copy()     
            between['z'] = self.phreeqc['z'] - 0.5*self.dz/(i+1)
            between = between[between['z'] > self.phreeqc['z'].min()] 
            for name in bkgName:
                between[name] = 0.5*array(self.phreeqc[name][:-(self.nx+1)]) + \
                    0.5*array(self.phreeqc[name][(self.nx+1):])
            self.phreeqc = concat([self.phreeqc, between], axis=0)
            
        # (7) created mirror image half on other side of fracture
        mirror = self.phreeqc.copy()
        mirror['x'] *= -1.
        self.phreeqc = concat([self.phreeqc, mirror], axis=0)

    def GridFracture(self):
        # extrude grid along y-axis
        self.phreeqc['y'] = self.y0 + 0.5*self.dy
        for i in xrange(1, self.ny):
            layer = self.phreeqc.copy()
            layer['y'] = self.y0 + (i+0.5)*self.dy
            if i==1: extrude = layer
            else: extrude = concat([extrude, layer], axis=0)
        self.phreeqc = concat([self.phreeqc, extrude], axis=0)
        # rotate grid in 3-D about origin
        fracPts = transpose([self.phreeqc['x'], self.phreeqc['y'], self.phreeqc['z']])
        rotPts = []          # container for rotated coordinate 3 x 1 vectors 
        for i in xrange(len(fracPts)): rotPts.append(dot(self.rotMatrix, fracPts[i]))
        rotPts = array(rotPts)
        prime = transpose(rotPts)
        self.phreeqc.insert(0, 'zPrime', prime[2])
        self.phreeqc.insert(0, 'yPrime', prime[1])
        self.phreeqc.insert(0, 'xPrime', prime[0])
        # shift along x-axis
        self.phreeqc['xPrime'] += self.grid_x_intercept

    def Cleanup(self, bkgList, grid):
        # alter dataframe to make compatible with background
        self.phreeqc.drop(['x', 'y', 'z'], axis=1, inplace=True)   # remove old coordinates
        self.phreeqc.rename(index=str, columns={'xPrime': 'x', 'yPrime': 'y', 'zPrime': 'z'}, inplace=True)
        self.phreeqc = self.phreeqc[bkgList]
        # clip fracture planes to match background grid
        self.phreeqc = self.phreeqc[(self.phreeqc['x']>=grid.x0) & (self.phreeqc['x']<=grid.xf)]
        self.phreeqc = self.phreeqc[(self.phreeqc['y']>=grid.y0) & (self.phreeqc['y']<=grid.yf)]
        self.phreeqc = self.phreeqc[(self.phreeqc['z']>=grid.z0) & (self.phreeqc['z']<=grid.zf)]


### supporintg functions ###
    
        
def ReadFracs():
    # read fractures files and accompanying PHREEQC output as data frames
    fracture = []
    lineInput = []        
    inputFile = open('frac_setup.txt','r')
    for line in inputFile: lineInput.append(line.split())
    inputFile.close()
    for i in xrange(1, len(lineInput)):
        fileName = lineInput[i][0]
        nz = int(lineInput[i][1])
        nx = int(lineInput[i][2])
        ny = int(lineInput[i][3])        
        dz = float(lineInput[i][4])
        dx = float(lineInput[i][5])
        dy = float(lineInput[i][6])
        w = float(lineInput[i][7])
        phi = float(lineInput[i][8])
        theta = float(lineInput[i][9])
        datum = float(lineInput[i][10])
        y0 = float(lineInput[i][11])
        grid_x_intercept = float(lineInput[i][12])
        nzDivide = int(lineInput[i][13])
        fracture.append(Fracture(fileName, nz, nx, ny, dz, dx, dy, w, phi, theta, datum, y0, grid_x_intercept, nzDivide))
    print('Read fracture attributes.')
    return fracture       
            

def ReadGrouped(fileName):
    # read list of minerals or kinetics items to be lumped together under 'fileName'
    name = fileName.split('.')[0]
    phases = []
    inputFile = open(fileName,'r')
    for line in inputFile: phases.append(line.split()[0])
    inputFile.close()
    print('Read', fileName)          
    return name, phases


def ReadInitial(fileName):
    # read initial mass of active phase to track loss
    phases = []
    molesZero = []
    lineInput = [] 
    inputFile = open(fileName,'r')
    for line in inputFile: lineInput.append(line.split())
    inputFile.close()
    for i in xrange(1, len(lineInput)):
        phases.append(lineInput[i][0])
        molesZero.append(float(lineInput[i][1]))        
    print('Read initial masses.')
    return phases, molesZero

            
### main script ###

    
def FracProcess():
    
    grid = Grid()                   # grid settings
    bkgList = ['x', 'y', 'z']       # list of column headings for final output
    bkgList.extend(grid.bkgName)
    fracture = ReadFracs()          # fracture properties
  
    # read files for lumping mineral phases/linetics items
    groupName = []
    phaseList = []
    groupedPhases = ['sec_Cu_sulfides.txt']
    for group in groupedPhases:
        name, phases = ReadGrouped(group)
        groupName.append(name)
        phaseList.append(phases) 

    # read initial masses       
    initPhases, molesZero = ReadInitial('ref_init.txt')
   
    # populate background grid and return dataframe    
    lattice = grid.FillGrid(fracture)       
    
    # process fracture spaces      
    for frac in fracture:
        print('Processing fracture', frac.name)
        frac.ProcessSimResults(groupName, phaseList, initPhases, molesZero, grid.bkgName)
        frac.GridFracture()
        frac.Cleanup(bkgList, grid)
        lattice = concat([lattice, frac.phreeqc], axis=0)   # append to results lattice
        
    lattice.to_csv('lattice.csv')    
    print('Done.')

    
### run script ###

FracProcess()