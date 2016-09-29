# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 23:22:20 2016

@author: jonas_000
"""

import calfem.core as cfc
import calfem.utils as cfu

import numpy as np
from scipy.sparse import lil_matrix

class Solver:
    def __init__(self, mesh):
        self.mesh = mesh
        self.nDofs = np.size(mesh.dofs)
        self.onInit()
        self.execute();
        
    def onInit(self):
        pass
               
    def execute(self):
        self.assem()
        self.applyBCs()
        self.applyLoads()

        print("Solving system...")
        
        self.a, self.r = cfc.spsolveq(self.K, self.f, self.bc, self.bcVal)
        
        print("Extracting ed...")
        
        self.ed = cfc.extractEldisp(self.mesh.edof, self.a)
        
        
    def assem(self):
        self.K = lil_matrix((self.nDofs, self.nDofs))

        print("Assembling K... ("+str(self.nDofs)+")")
              
        for eltopo, elx, ely in zip(self.mesh.edof, self.mesh.ex, self.mesh.ey):
        
            Ke = self.onCreateKe(elx, ely, self.mesh.shape.elementType)                
            cfc.assem(eltopo, self.K, Ke)
            
    def applyBCs(self):
        print("Applying bc and loads...")
        
        self.bc = np.array([],'i')
        self.bcVal = np.array([],'i')
        
        self.onApplyBCs(self.mesh, self.bc, self.bcVal)
        
    def applyLoads(self):
        self.f = np.zeros([self.nDofs,1])
        
        self.onApplyLoads(self.mesh, self.f)
        
    def calcElementForces(self):
        print("Element forces... ")
        
        for i in range(self.mesh.edof.shape[0]):
            self.onCalcElForce(self.mesh.ex[i,:], self.mesh.ey[i,:], self.ed[i,:])
            
    def onCalcElForce(self, ex, ey, ed):
        pass

    def onCreateKe(self, elx, ely, elementType):
        pass
    
    def onApplyBCs(self, mesh, bc, bcVal):        
        pass
        
    def onApplyLoads(self, mesh, f):
        pass
        



# ---- Calculate elementr stresses and strains ------------------------------

