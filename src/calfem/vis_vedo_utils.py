# -*- coding: utf-8 -*-
"""
CALFEM Vedo utility module

Utility functions for Vedo
"""

import numpy as np
from vedo import *
import vtk
import pyvtk
import sys

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# By Jonas Lindemann

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

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
# By Andreas Åmand

def get_coord_from_edof(edof_row,dof,element_type):
    """
    Routine to get element coodinates based on element type and degrees of freedom.

    :param array edof_row: Element topology row [1 x degrees of freedom per element]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param int element_type: Element type [1-6]

    :return array coords: Array of node coordinates for element [n_nodes x 3]
    """
    if element_type == 1 or element_type == 2 or element_type == 5:
        edof_row1,edof_row2 = np.split(edof_row,2)
        coord1 = int(np.where((edof_row1==dof).any(axis=1))[0])
        coord2 = int(np.where((edof_row2==dof).any(axis=1))[0])
        return coord1, coord2
    elif element_type == 3 or element_type == 4:
        edof_row1,edof_row2,edof_row3,edof_row4,edof_row5,edof_row6,edof_row7,edof_row8 = np.split(edof_row,8)
        coord1 = int(np.where(np.any(edof_row1==dof,axis=1))[0])
        coord2 = int(np.where(np.any(edof_row2==dof,axis=1))[0])
        coord3 = int(np.where(np.any(edof_row3==dof,axis=1))[0])
        coord4 = int(np.where(np.any(edof_row4==dof,axis=1))[0])
        coord5 = int(np.where(np.any(edof_row5==dof,axis=1))[0])
        coord6 = int(np.where(np.any(edof_row6==dof,axis=1))[0])
        coord7 = int(np.where(np.any(edof_row7==dof,axis=1))[0])
        coord8 = int(np.where(np.any(edof_row8==dof,axis=1))[0])
        coords = np.array([coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8])
        return coords
    elif element_type == 6:
        edof_row1,edof_row2,edof_row3,edof_row4 = np.split(edof_row,4)
        coord1 = int(np.where(np.any(edof_row1==dof,axis=1))[0])
        coord2 = int(np.where(np.any(edof_row2==dof,axis=1))[0])
        coord3 = int(np.where(np.any(edof_row3==dof,axis=1))[0])
        coord4 = int(np.where(np.any(edof_row4==dof,axis=1))[0])
        coords = np.array([coord1, coord2, coord3, coord4])
        return coords

def get_a_from_coord(coord_row_num,num_of_deformations,a,scale=1):
    """
    Routine to get node displacements based on coordinates.

    :param int coord_row_num: Node coordinate row number
    :param int num_of_deformations: Number of degrees of freedom per node
    :param array a: Global displacement vector [1 x total degrees of freedom]

    :return float dx: Nodal displacement in x-direction
    :return float dy: Nodal displacement in y-direction
    :return float dz: Nodal displacement in z-direction
    """
    dx = a[coord_row_num*num_of_deformations]*scale
    dy = a[coord_row_num*num_of_deformations+1]*scale
    dz = a[coord_row_num*num_of_deformations+2]*scale
    return dx, dy, dz

def get_node_elements(coord,scale,alpha,dof,bcPrescr=None,bc=None,bc_color='red',fPrescr=None,f=None,f_color='blue6',dofs_per_node=None):
    """
    Routine to get node node actors.

    :param array coord: Nodal coordinates [number of nodes x 3]
    :param int scale: Node actor radius
    :param float alpha: Node actor transparency [0-1]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param array bcPrescr: Degrees of freedom with prescribed boundary conditions [number of prescribed boundary contidions x 1]
    :param array bc: Values for prescribed boundary conditions [number of prescribed boundary contidions x 1]
    :param string bc_color: Color for nodes with prescribed boundary conditions
    :param array fPrescr: Degrees of freedom with applied forces [number of applied forces x 1]
    :param array f: Values for forces [number of applied forces x 1]
    :param string f_color: Color for nodes with applied forces
    :param int dofs_per_node: Degrees of freedom per node [1-6]

    :return list nodes: Node actors
    """
    nnode = np.size(coord, axis = 0)
    ncoord = np.size(coord, axis = 1)
    nodes = []

    bc_dict = {}
    indx = 0
    if isinstance(bcPrescr, np.ndarray):
        for i in bcPrescr:
            bc_dict[i] = bc[indx]
            indx += 1

    f_dict = {}
    indx = 0
    if isinstance(fPrescr, np.ndarray):
        for i in fPrescr:
            f_dict[i] = f[indx]
            indx += 1

    for i in range(nnode):

        dofs = dof[i]

        if np.any(np.isin(bcPrescr, dofs, assume_unique=True)) == True:
            color = bc_color
        elif np.any(np.isin(fPrescr, dofs, assume_unique=True)) == True:
            color = f_color
        else:
            color = 'black'

        node = Sphere(c=color).scale(1.5*scale).pos([coord[i,0],coord[i,1],coord[i,2]]).alpha(alpha)

        if np.any(np.isin(bcPrescr, dofs, assume_unique=True)) == True:
            node.name = f"Node nr. {i+1}, DoFs & BCs: ["
            for j in range(dofs_per_node):
                node.name += str(dof[i,j])
                if dof[i,j] in bc_dict:
                    node.name += (': ' + str(bc_dict[dof[i,j]]))
                if j == dofs_per_node-1:
                    node.name += ']'
                else:
                    node.name += ', '

        elif np.any(np.isin(fPrescr, dofs, assume_unique=True)) == True:
            node.name = f"Node nr. {i+1}, DoFs & Forces: ["
            for j in range(dofs_per_node):
                node.name += str(dof[i,j])
                if dof[i,j] in f_dict:
                    node.name += (': ' + str(f_dict[dof[i,j]]))
                if j == dofs_per_node-1:
                    node.name += ']'
                else:
                    node.name += ', '

        else:
            node.name = f"Node nr. {i+1}, DoFs: ["
            for j in range(dofs_per_node):
                node.name += str(dof[i,j])
                if j == dofs_per_node-1:
                    node.name += ']'
                else:
                    node.name += ', '

        nodes.append(node)

    return nodes

def vectors(
    points,
    vectors,
    c="k",
    alpha=1,
    shaftWidth=0.05,
    text=None, vmax=None, vmin=None, cmap='jet', values = None):
    """
    Routine for creating vectors.

    :param list points: Mid point for vector [number of vectors x 3]
    :param list vectors: Vector components [number of vectors x 3]
    :param string dof: Vector color
    :param float alpha: Vector transparancy [0-1]
    :param float shaftWidth: Vector width
    :param list text: Vector values [number of vectors x 1]
    :param float vmax: Maximum vector value for colormapping
    :param float vmin: Minimum vector value for colormapping
    :param string cmap: Vector colormap
    :param list values: [number of vectors x 1]

    :return list cylinders: Vector actors
    """
    if isinstance(points, Points):
        points = points.points()
    else:
        points = np.array(points)
    vectors = np.array(vectors) / 2

    spts = points - vectors
    epts = points + vectors

    npts = np.size(points,0)

    cylinders = []
    for i in range(npts):
        cyl = Cylinder([spts[i],epts[i]],r=shaftWidth*0.01,res=4,c=c)
        cyl.name = text[i]
        cyl.cmap(cmap,input_array=values[i],vmin=vmin,vmax=vmax,on='cells')
        cylinders.append(cyl)

    return cylinders

def check_input(edof,coord,dof,element_type,a=None,values=None,nseg=None):
    """
    Routine for checking input to draw_mesh, draw_displaced_mesh & animation.

    :param array edof: Element topology [number of elements x degrees of freedom per element]
    :param array coord: Nodal coordinates [number of nodes x 3]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param int element_type: Element type [1-6]
    :param array a: Global displacement vector [degrees of freedom x 1]
    :param array values: Scalar values [number of elements x 1 | number of elements x nodes per element | number of nodes x 1]
    :param int nseg: Number of beam segments + 1

    :return int number_of_elements: Number of elements
    :return int number_of_degrees_of_freedom_per_element: Degrees of freesom per element
    :return int number_of_coordinates: Number of coordinates
    :return int number_of_dimensions: Number of dimensions for model [1-3]
    :return int number_of_degrees_of_freedom: Number of degrees of freedom
    :return int degrees_of_freedom_per_node: Degrees of freedom per node
    :return int number_of_displacements: Number of displacements
    :return string val: Types of scalar input ['el_values' / 'nodal_values_by_el' / 'nodal_values']
    """

    if element_type == 1 or element_type == 2 or element_type == 5:
        number_of_nodes_per_element = 2
    elif element_type == 3 or element_type == 4:
        number_of_nodes_per_element = 8
    elif element_type == 6:
        number_of_nodes_per_element = 4

    number_of_elements = np.size(edof, axis=0)
    number_of_degrees_of_freedom_per_element = np.size(edof, axis=1)

    number_of_coordinates = np.size(coord, axis=0)
    number_of_dimensions = np.size(coord, axis=1)

    number_of_degrees_of_freedom = np.size(dof, axis=0)*np.size(dof, axis=1)
    degrees_of_freedom_per_node = np.size(dof, axis=1)

    if element_type == 1:
        element_name = 'Spring'
    elif element_type == 2:
        element_name = 'Bar'
    elif element_type == 3:
        element_name = 'Flow'
    elif element_type == 4:
        element_name = 'Solid'
    elif element_type == 5:
        element_name = 'Beam'
    elif element_type == 6:
        element_name = 'Plate'

    # Checks displacements, degrees of freedom
    if a is not None:
        number_of_displacements = np.size(a, axis=0)
        if number_of_displacements != number_of_degrees_of_freedom:
            print(element_name + " element: Number of displacements (a) and total number of degrees of freedom (Dof) does not correspond")
            sys.exit()

    # Checks scalar values, edof
    if values is not None:
        if element_type == 1 or element_type == 2 or element_type == 5:
            if element_type == 1 or element_type == 2:
                nseg = 1
            number_of_values = np.size(values, axis=0)
            if number_of_values == number_of_elements*nseg:
                val = 'el_values'
            else:
                print(element_name + " element: Invalid number of element-/nodal scalars, please make sure the number of scalar values correspond to total number of elements/nodes")
                sys.exit()
            
        else:
            number_of_values = np.size(values, axis=0)*np.size(values, axis=1)

            if number_of_values == number_of_elements:
                val = 'el_values'
            elif number_of_values == number_of_elements*number_of_nodes_per_element:
                val = 'nodal_values_by_el'
            elif number_of_values == number_of_coordinates:
                val = 'nodal_values'
            else:
                print(element_name + " element: Invalid number of element-/nodal scalars, please make sure the number of scalar values correspond to total number of elements/nodes")
                sys.exit()

    # Checks coord, dof
    if np.size(coord, axis=0) != np.size(dof, axis=0):
        print(element_name + " element: Number of rows in Coord & Dof does not correspond, please check them along with number of nodes")
        sys.exit()

    # Checks edof
    if element_type == 1:
        if number_of_degrees_of_freedom_per_element != 2:
            print(element_name + " element: 2 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 1:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 1 required")
            sys.exit()

    # Checks edof
    elif element_type == 2:
        if number_of_degrees_of_freedom_per_element != 6:
            print(element_name + " element: 6 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 3:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 3 required")
            sys.exit()

    # Checks edof
    elif element_type == 3:
        if number_of_degrees_of_freedom_per_element != 8:
            print(element_name + " element: 8 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 1:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 1 required")
            sys.exit()

    # Checks edof
    elif element_type == 4:
        if number_of_degrees_of_freedom_per_element != 24:
            print(element_name + " element: 24 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 3:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 3 required")
            sys.exit()

    # Checks edof
    elif element_type == 5:
        if number_of_degrees_of_freedom_per_element != 12:
            print(element_name + " element: 12 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 6:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 6 required")
            sys.exit()

    # Checks edof
    elif element_type == 6:
        if number_of_degrees_of_freedom_per_element != 12:
            print(element_name + " element: 12 degrees of freedom per element required")
            sys.exit()
        if degrees_of_freedom_per_node < 3:
            print(element_name + " element: Too few degrees of freedom per node in Dof, at least 3 required")
            sys.exit()

    if a is None and values is None:
        return number_of_elements, \
            number_of_degrees_of_freedom_per_element, \
            number_of_coordinates, \
            number_of_dimensions, \
            number_of_degrees_of_freedom, \
            degrees_of_freedom_per_node,
    elif a is None:
        return number_of_elements, \
            number_of_degrees_of_freedom_per_element, \
            number_of_coordinates, \
            number_of_dimensions, \
            number_of_degrees_of_freedom, \
            degrees_of_freedom_per_node, \
            val
    elif values is None:
        return number_of_elements, \
            number_of_degrees_of_freedom_per_element, \
            number_of_coordinates, \
            number_of_dimensions, \
            number_of_degrees_of_freedom, \
            degrees_of_freedom_per_node, \
            number_of_displacements
    else:
        return number_of_elements, \
            number_of_degrees_of_freedom_per_element, \
            number_of_coordinates, \
            number_of_dimensions, \
            number_of_degrees_of_freedom, \
            degrees_of_freedom_per_node, \
            number_of_displacements, \
            val
            
def convert_nodal_values(edof,topo,dof,values):
    """
    Routine for converting nodal values from element to global and interpolating.

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array topo: element topology [nel x nodes per element]
    :param array dof: Global degrees of freedom [number of nodes x degrees of freedom per node]
    :param array values: Element scalar values [nel x nodes per element]

    :return array nodal_value_array: Global scalar values at nodes
    """
    nel = np.size(edof, axis=0)
    nnode = np.size(dof, axis=0)
    nodal_value_array = np.zeros((nnode,1))
    #print('Number of element values: ', np.size(values, axis=0))
    #print('Number of values per element: ', np.size(values, axis=1))

    topo_num = {}
    for i in range(nel):
        topo_num[i] = tuple(topo[i])

    test = {}
    for el, nodes in topo_num.items():
        it = 0
        points = [0,1,2,3,4,5,6,7]

        for i in nodes:
            test[tuple([el,points[it]])] = i
            it += 1

    test2 = {}
    for i in range(nnode):
        test2[i] = []

    for data, node in test.items():
        test2[node].append(values[data])

    for i in range(nnode):
        nodal_value_array[i] = np.mean(test2[i])

    return nodal_value_array

def convert_el_values(edof,values):
    """
    Routine for converting element values from element to global.

    :param array edof: element topology by degrees of freedom [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array values: Element scalar values [nel x nodes per element]

    :return array el_values: Global scalar values for element
    """
    nel = np.size(edof, axis=0)
    el_values = np.zeros((nel*6,1))

    for i in range(nel):
        el_values[i*6,:] = values[i]
        el_values[i*6+1,:] = values[i]
        el_values[i*6+2,:] = values[i]
        el_values[i*6+3,:] = values[i]
        el_values[i*6+4,:] = values[i]
        el_values[i*6+5,:] = values[i]
    return el_values
