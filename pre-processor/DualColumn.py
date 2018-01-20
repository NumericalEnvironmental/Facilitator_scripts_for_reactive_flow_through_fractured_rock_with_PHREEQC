############################################################################
#
# DualColumn.py
#
# python script to setup PHREEQC dual-porosity fracture flow simulation
#
# by Walt McNab
#
############################################################################

from __future__ import print_function
from numpy import *
from scipy import special

phi = 0.3       # default porosity
dy = 1.         # assumed thickness = 1 m


class PropBlock:
    
    def __init__(self, fileName):
        # populate data set (e.g., equilibrium_phases, kinetics)
        self.name = fileName
        self.lineInput = []        
        inputFile = open(fileName,'r')
        for line in inputFile: self.lineInput.append(line)
        inputFile.close()        
        self.bounds = self.lineInput[0].split()
        
    def WriteBlock(self, nx, nz):
        # write appropriately numbered blocks to file
        outputFile = open('blocks_' + self.name,'w')
        outputFile.writelines([self.bounds[0], '\t', self.bounds[1] + '-' + self.bounds[2], '\n'])
        for item in self.lineInput[1:]: outputFile.writelines([item])
        outputFile.writelines(['\n'])        
        for jMatrix in range(1, nx+1):
            outputFile.writelines([self.bounds[0], '\t', str(int(self.bounds[1])+(jMatrix*nz)+1) + '-' + str(int(self.bounds[2])+(jMatrix*nz)+1), '\n'])
            for item in self.lineInput[1:]: outputFile.writelines([item])
            outputFile.writelines(['\n'])
        outputFile.close()        


class Params:

    def __init__(self):
        # read model setup parameters
        lineInput = []        
        inputFile = open('grid_params.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.De = float(lineInput[0][1])                 # effective diffusion coefficient
        self.nz = int(lineInput[1][1])                   # number of advective cells along 1-D column
        self.nx = int(lineInput[2][1])                   # number of diffusive cells transverse to 1-D column
        self.dz = float(lineInput[3][1])                 # length of cells along 1-D column	
        self.dx = float(lineInput[4][1])                 # length of cells transverse to 1-D column
        self.w = float(lineInput[5][1])                 	# length of cells transverse to 1-D column
        self.dt = float(lineInput[6][1])*86400.          # time step size (convert from days to seconds)
        self.nShifts = int(lineInput[7][1])              # number of shifts (nShifts*dt = end time)
        self.alpha = float(lineInput[8][1])              # dispersivity along 1-D column (length)
        self.theta = float(lineInput[9][1]) * pi/180.    # fracture departure from vertical
        self.tEnd = self.dt*self.nShifts                 # implied simulation end time
        print('Read parameters.')
        print('End time (years) = ', self.tEnd/(86400.*365.))    
        q = phi * self.dz/self.dt               # Darcy velocity
        K = q/cos(self.theta)
        print ('Hydraulic conductivity (m/day) = ', K * 86400.)

    def WriteTransport(self):
        # write transport keyword block
        outputFile = open('transport.txt','w')
        outputFile.writelines(['TRANSPORT', '\n'])
        outputFile.writelines(['\t', '-cells', '\t', str(self.nz), '\n'])    
        outputFile.writelines(['\t', '-shifts', '\t', str(self.nShifts), '\n'])    
        outputFile.writelines(['\t', '-time_step', '\t', str(self.dt), '\n'])   
        outputFile.writelines(['\t', '-lengths', '\t', str(self.dz), '\n'])        
        outputFile.writelines(['\t', '-dispersivities', '\t', str(self.alpha), '\n'])
        outputFile.writelines(['\t', '-boundary_conditions', '\t', 'flux', '\t', 'flux', '\n'])
        outputFile.writelines(['\t', '-stagnant', '\t', str(self.nx), '\n'])           
        outputFile.writelines(['\t', '-correct_disp', '\t', 'true', '\n'])  
        outputFile.writelines(['\t', '-print_frequency', '\t', str(self.nShifts), '\n'])     
        outputFile.writelines(['\t', '-punch_frequency', '\t', str(self.nShifts), '\n']) 
        outputFile.close()
        print('Wrote transport keyword block to file.')   
       
    def WriteMix(self):
        # definitions
        dxb = 0.5*self.dx + 0.5*self.w
        A = self.dz * dy
        volFrac = A * self.w
        volMatrix = A * self.dx
        JFrac = self.De * self.dt * A / dxb
        JMatrix = self.De * self.dt * A / self.dx    
        mixStart = [1.-JFrac/volFrac, JFrac/volFrac]        # cell fracture within fracture
        mix2nd = [JFrac/volMatrix, 1.-JFrac/volMatrix-JMatrix/volMatrix, JMatrix/volMatrix]    # border cell: in matrix, but adjacent to fracture
        mixMiddle = [JMatrix/volMatrix, 1.-2*JMatrix/volMatrix, JMatrix/volMatrix]   # interior matrix cells
        mixEnd = [JMatrix/volMatrix, 1.-JMatrix/volMatrix] # end of matrix cell column         
        # write mix keyword block
        outputFile = open('mix.txt','w')    
        for i in range(1, self.nz+1):
            outputFile.writelines(['MIX', '\t', str(i), ';', '\t',  str(i) , '\t',  str(mixStart[0]), ';', '\t',  str(i+self.nz+1), '\t', str(mixStart[1]), '\n'])
            outputFile.writelines(['MIX', '\t', str(i+self.nz+1), ';', '\t',  str(i) , '\t',  str(mix2nd[0]), ';', '\t',  str(i+self.nz+1), '\t', str(mix2nd[1]), ';', '\t',  str(i+2*self.nz+1), '\t', str(mix2nd[2]),'\n'])
            for j in range(2, self.nx):
                outputFile.writelines(['MIX', '\t', str(i+j*self.nz+1), ';', '\t',  str(i+(j-1)*self.nz+1), '\t',  str(mixMiddle[0]), ';', '\t',  str(i+j*self.nz+1), '\t', str(mixMiddle[1]), ';', '\t',  str(i+(j+1)*self.nz+1), '\t', str(mixMiddle[2]),'\n'])
            outputFile.writelines(['MIX', '\t', str(i+self.nx*self.nz+1), ';', '\t', str(i+(self.nx-1)*self.nz+1) , '\t',  str(mixEnd[0]), ';', '\t',  str(i+self.nx*self.nz+1), '\t', str(mixEnd[1]), '\n']) 
            outputFile.writelines(['\n']) 
        outputFile.close()
        x = arange(0.5*self.dx, (self.nx+0.5)*self.dx, self.dx)
        C = special.erfc(x/(2.*sqrt(self.De*self.tEnd)))
        print('Spacing from fracture = ', x)
        print('Normalized concentration at simulation end = ', C)
        print('Wrote mix keyword block to file.')    

    
### main script ###

    
def dualPoros():
    
    params = Params()           # read parameter input
    params.WriteMix()           # write mix keyword block
    params.WriteTransport()     # write transport keyword block
                  
    blockFiles = ['eq_phases_top.txt', 'eq_phases_bottom.txt', 'kinetics_top.txt']
    for block in blockFiles:
        blockObj = PropBlock(block)
        blockObj.WriteBlock(params.nx, params.nz)
    print('Wrote to block files.')
                         
    print('Done.')
    
### run script ###

dualPoros()