# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 23:22:20 2016

@author: jonas_000
"""

import calfem.core as cfc
import calfem.utils as cfu
import logging as cflog

import numpy as np
from scipy.sparse import lil_matrix

def error(msg):
    cflog.error(" calfem.solver: "+msg)

def info(msg):
    cflog.info(" calfem.solver: "+msg)


class Results:
    pass

class Solver:
    def __init__(self, mesh):
        self.results = Results()
        self.mesh = mesh
        self.nDofs = np.size(mesh.dofs)
        self.nElements = np.size(self.mesh.edof,0)

        self.bc = np.array([],'i')
        self.bcVal = np.array([],'i')
        self.f = np.zeros([self.nDofs,1])
        
        self.results.elForces = np.zeros([self.nElements, self.onQueryElForceSize()])
        
    def onQueryElForceSize(self):
        return 1
                       
    def execute(self):
        info("Assembling K... ("+str(self.nDofs)+")")
        self.assem()
        
        info("Solving system...")        
        self.results.a, self.results.r = cfc.spsolveq(self.K, self.f, self.bc, self.bcVal)
        
        info("Extracting ed...")        
        self.results.ed = cfc.extractEldisp(self.mesh.edof, self.results.a)
        
        info("Element forces... ")
        self.calcElementForces()
        
        return self.results
        
    def assem(self):
        self.K = lil_matrix((self.nDofs, self.nDofs))
        for eltopo, elx, ely in zip(self.mesh.edof, self.mesh.ex, self.mesh.ey):
            Ke = self.onCreateKe(elx, ely, self.mesh.shape.elementType)                
            cfc.assem(eltopo, self.K, Ke)
            
    def addBC(self, marker, value=0.0, dimension=0):
        self.bc, self.bcVal = cfu.applybc(self.mesh.bdofs, self.bc, self.bcVal, marker, value, dimension)
        
    def addForceTotal(self, marker, value=0.0, dimension=0):
        cfu.applyforcetotal(self.mesh.bdofs, self.f, self.mesh.shape.topId, value, dimension)
          
    def addForce(self, marker, value=0.0, dimension=0):
        cfu.applyforce(self.mesh.bdofs, self.f, self.mesh.shape.topId, value, dimension)
        
    def addForceNode(self, node, value = 0.0, dimension=0):
        cfu.applyforcenode(node, value, dimension)
        
    def addBCNode(self, node, value = 0.0, dimension = 0):
        self.bc, self.bcVal = cfu.applybcnode(node, value, dimension)

    def applyBCs(self):
        self.bc, self.bcVal = self.onApplyBCs(self.mesh, self.bc, self.bcVal)
                
    def calcElementForces(self):
        for i in range(self.mesh.edof.shape[0]):
            elForce = self.onCalcElForce(self.mesh.ex[i,:], self.mesh.ey[i,:], self.results.ed[i,:], self.mesh.shape.elementType)
            if len(elForce)==1:
                self.results.elForces[i,:] = elForce
            else:
                pass
            
            
    def onCalcElForce(self, ex, ey, ed, elementType):
        pass

    def onCreateKe(self, elx, ely, elementType):
        pass
    
    def onApplyBCs(self, mesh, bc, bcVal):        
        pass
        
    def onApplyLoads(self, mesh, f):
        pass
        
class Plan2DSolver(Solver):
        
    def onCreateKe(self, elx, ely, elementType):
        Ke = None
        if self.mesh.shape.elementType == 2:
            Ke = cfc.plante(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
        else:
            Ke = cfc.planqe(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
            
        return Ke
                    
    def onCalcElForce(self, ex, ey, ed, elementType):
        if elementType == 2: 
            es, et = cfc.plants(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
            elMises = np.math.sqrt( pow(es[0,0],2) - es[0,0]*es[0,1] + pow(es[0,1],2) + 3*pow(es[0,2],2) )
        else:
            es, et = cfc.planqs(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
            elMises = np.math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*pow(es[2],2) )
        
        return elMises

class Flow2DSolver(Solver):
        
    def onCreateKe(self, elx, ely, elementType):
        Ke = None
        if self.mesh.shape.elementType == 2:
            Ke = cfc.flw2te(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
        else:
            Ke = cfc.flw2i4e(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
            
        return Ke
                    
    def onCalcElForce(self, ex, ey, ed, elementType):
        es = None
        et = None
        if elementType == 2: 
            es, et = cfc.flw2ts(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
        else:
            es, et, temp = cfc.flw2i4s(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
        
        return [es, et]
