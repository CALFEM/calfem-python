# -*- coding: utf-8 -*-
"""
CALFEM Vedo utility module

Utility functions for Vedo
"""

import numpy as np
from vedo import *
import vtk

def von_mises_3d(es):
    
    sigv = np.zeros(es.shape[0])
    
    for i, es_ip in enumerate(es):
        sigx, sigy, sigz, sigxy, sigyz, sigxz = es_ip
        
        sigv[i] = sqrt(  0.5*((sigx - sigy)**2 + (sigx - sigz)**2 + (sigz - sigx)**2)
                       + 3.0*(sigxy**2 + sigyz**2 + sigxz**2))
        
    return sigv

def sigv_to_hex(es):
    
    n_ip = es.shape[0]
    
    if n_ip == 1:
        return np.ones(8)*es
    elif n_ip == 8:
        return es
    elif n_ip == 27:
        return es[np.array([0, 2, 6, 8, 18, 20, 24, 26])]
    else:
        return None


def ugrid_from_edof_ec(edof, ex, ey, ez, ed=None, dofs_per_node=3, ignore_first=False):
    
    coords, topo, node_dofs, node_displ = convert_to_node_topo(edof, ex, ey, ez, ed, dofs_per_node, ignore_first)

    npoint = coords.shape[0]
    nel = topo.shape[0]
    nnd = topo.shape[1]

    if nnd == 4:
        ct = vtk.VTK_TETRA
    elif nnd == 8:
        ct = vtk.VTK_HEXAHEDRON
    else:
        print("Topology not supported.")

    celltypes = [ct] * nel

    return UGrid([coords, topo, celltypes])

def convert_to_node_topo(edof, ex, ey, ez, ed=None, es=None, dofs_per_node=3, ignore_first=False):
    """
    Routine to convert dof based topology and element coordinates to node based
    topology required for visualisation with VTK and other visualisation frameworks

    :param array edof: element topology [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array ex: element x coordinates [nel x n_nodes]
    :param array ey: element y coordinates [nel x n_nodes]
    :param array ez: element z coordinates [nel x n_nodes]
    :param array n_dofs_per_node: number of dofs per node. (default = 3)
    :param boolean ignore_first: ignore first column of edof. (default = True)
    :return array coords: Array of node coordinates. [n_nodes x 3]
    :return array topo: Node topology. [nel x n_nodes]
    :return array node_dofs: Dofs for each node. [n_nodes x n_dofs_per_node]
    """

    node_hash_coords = {}
    node_hash_numbers = {}
    node_hash_dofs = {}
    node_hash_displ = {}
    node_hash_scalar = {}
    node_hash_count = {}
    el_hash_dofs = []


    nel, cols = edof.shape

    if ignore_first:
        tot_dofs = cols-1
    else:
        tot_dofs = cols

    n_nodes = int(tot_dofs / dofs_per_node)

    # print("n_dofs_per_node =", dofs_per_node)
    # print("cols    =", tot_dofs)
    # print("nel     =", nel)
    # print("n_nodes =", n_nodes)

    if ed is None:
        ed = np.zeros((nel,n_nodes*dofs_per_node))

    if es is None:
        es = np.zeros((nel,n_nodes))

    for elx, ely, elz, eed, ees, dofs in zip(ex, ey, ez, ed, es, edof):

        if ignore_first:
            el_dofs = dofs[1:]
        else:
            el_dofs = dofs

        # 0 1 2  3 4 5  6 7 8  9 12 11 

        el_dof = np.zeros((n_nodes, dofs_per_node), dtype=int)
        el_hash_topo = []

        for i in range(n_nodes):
            el_dof[i] = el_dofs[ (i*dofs_per_node):((i+1)*dofs_per_node) ]
            node_hash_coords[hash(tuple(el_dof[i]))] = [elx[i], ely[i], elz[i]]
            node_hash_numbers[hash(tuple(el_dof[i]))] = -1
            node_hash_dofs[hash(tuple(el_dof[i]))] = el_dof[i]

            if hash(tuple(el_dof[i])) in node_hash_scalar:
                node_hash_scalar[hash(tuple(el_dof[i]))] += ees[i]
            else:
                node_hash_scalar[hash(tuple(el_dof[i]))] = ees[i]

            if hash(tuple(el_dof[i])) in node_hash_count:
                node_hash_count[hash(tuple(el_dof[i]))] += 1
            else:
                node_hash_count[hash(tuple(el_dof[i]))] = 1

            displ_dofs = []
            for j in range(dofs_per_node):
                displ_dofs.append(eed[i*dofs_per_node+j])
            node_hash_displ[hash(tuple(el_dof[i]))] = displ_dofs
            el_hash_topo.append(hash(tuple(el_dof[i])))

        el_hash_dofs.append(el_hash_topo)

    coord_count = 0

    coords = []
    node_dofs = []
    node_displ = []
    node_scalars = []

    for node_hash in node_hash_numbers.keys():
        node_hash_numbers[node_hash] = coord_count
        node_dofs.append(node_hash_dofs[node_hash])
        node_displ.append(node_hash_displ[node_hash])
        node_scalars.append(node_hash_scalar[node_hash]/node_hash_count[node_hash])
        coord_count +=1

        coords.append(node_hash_coords[node_hash])

    topo = []

    for el_hashes in el_hash_dofs:

        el_hash_topo = []

        for el_hash in el_hashes:
            el_hash_topo.append(node_hash_numbers[el_hash])

        topo.append(el_hash_topo)

    return np.asarray(coords), np.asarray(topo), np.asarray(node_dofs), np.asarray(node_displ), np.asarray(node_scalars)