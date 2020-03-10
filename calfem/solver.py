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
        self.n_dofs = np.size(mesh.dofs)
        self.n_elements = np.size(self.mesh.edof,0)

        self.bc = np.array([],'i')
        self.bc_val = np.array([],'i')
        self.f = np.zeros([self.n_dofs,1])
        
        self.results.el_forces = np.zeros([self.n_elements, self.on_query_el_force_size()])
        
    def on_query_el_force_size(self):
        return 1
                       
    def execute(self):
        info("Assembling K... ("+str(self.n_dofs)+")")
        self.assem()
        
        info("Solving system...")        
        self.results.a, self.results.r = cfc.spsolveq(self.K, self.f, self.bc, self.bc_val)
        
        info("Extracting ed...")        
        self.results.ed = cfc.extractEldisp(self.mesh.edof, self.results.a)
        
        info("Element forces... ")
        self.calc_element_forces()
        
        return self.results
        
    def assem(self):
        self.K = lil_matrix((self.n_dofs, self.n_dofs))
        for eltopo, elx, ely in zip(self.mesh.edof, self.mesh.ex, self.mesh.ey):
            Ke = self.on_create_Ke(elx, ely, self.mesh.shape.element_type)                
            cfc.assem(eltopo, self.K, Ke)
            
    def addBC(self, marker, value=0.0, dimension=0):
        self.bc, self.bc_val = cfu.applybc(self.mesh.bdofs, self.bc, self.bc_val, marker, value, dimension)
        
    def addForceTotal(self, marker, value=0.0, dimension=0):
        cfu.applyforcetotal(self.mesh.bdofs, self.f, self.mesh.shape.top_id, value, dimension)
          
    def addForce(self, marker, value=0.0, dimension=0):
        cfu.applyforce(self.mesh.bdofs, self.f, self.mesh.shape.top_id, value, dimension)
        
    def addForceNode(self, node, value = 0.0, dimension=0):
        cfu.applyforcenode(node, value, dimension)
        
    def addBCNode(self, node, value = 0.0, dimension = 0):
        self.bc, self.bc_val = cfu.applybcnode(node, value, dimension)

    def applyBCs(self):
        self.bc, self.bc_val = self.on_apply_bcs(self.mesh, self.bc, self.bcVal)
                
    def calc_element_forces(self):
        for i in range(self.mesh.edof.shape[0]):
            el_force = self.on_calc_el_force(self.mesh.ex[i,:], self.mesh.ey[i,:], self.results.ed[i,:], self.mesh.shape.element_type)
            if el_force!=None:
                self.results.el_forces[i,:] = el_force
            else:
                pass
            
            
    def on_calc_el_force(self, ex, ey, ed, element_type):
        pass

    def on_create_Ke(self, elx, ely, element_type):
        pass
    
    def on_apply_bcs(self, mesh, bc, bcVal):        
        pass
        
    def on_apply_loads(self, mesh, f):
        pass
        
class Plan2DSolver(Solver):
        
    def on_create_Ke(self, elx, ely, element_type):
        Ke = None
        if self.mesh.shape.element_type == 2:
            Ke = cfc.plante(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
        else:
            Ke = cfc.planqe(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
            
        return Ke
                    
    def on_calc_el_force(self, ex, ey, ed, element_type):
        if element_type == 2: 
            es, et = cfc.plants(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
            elMises = np.math.sqrt( pow(es[0,0],2) - es[0,0]*es[0,1] + pow(es[0,1],2) + 3*pow(es[0,2],2) )
        else:
            es, et = cfc.planqs(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
            elMises = np.math.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*pow(es[2],2) )
        
        return elMises

class Flow2DSolver(Solver):
        
    def on_create_Ke(self, elx, ely, element_type):
        Ke = None
        if self.mesh.shape.element_type == 2:
            Ke = cfc.flw2te(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
        else:
            Ke = cfc.flw2i4e(elx, ely, self.mesh.shape.ep, self.mesh.shape.D)
            
        return Ke
                    
    def on_calc_el_force(self, ex, ey, ed, element_type):
        es = None
        et = None
        if element_type == 2: 
            es, et = cfc.flw2ts(ex, ey, self.mesh.shape.D, ed)
        else:
            es, et, temp = cfc.flw2i4s(ex, ey, self.mesh.shape.ep, self.mesh.shape.D, ed)
        
        return [es, et]
