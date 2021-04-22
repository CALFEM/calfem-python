import os
import sys
import tempfile
import shutil
import subprocess
import numpy as np
from calfem.core import createdofs
from calfem.utils import which
import calfem.core as cfc

import logging as cflog


def error(msg):
    """Log error message"""
    cflog.error(msg)


def info(msg):
    """Log information message"""
    cflog.info(msg)


def cmp(a, b):
    return (a > b) ^ (a < b)

# def dofsFromNodes(listOfNodes, dofs):
#        D = []
#        for node in listOfNodes:
#            D.extend(dofs[node])
#        return D


def _offsetIndices(lst, offset=0):
    '''Shifts the indices by offset. 
    Positive offsets move the indices away from 0.
    Negative offsets move the indices towards 0.
    If an index is 0 the offset is added to it.'''
    return [x + cmp(x, -0.5)*offset for x in lst]


def _formatList(lst, offset=0):
    """
    Turns a list of numbers into a corresponding string of comma-separated numbers.
    The parameter offset is a number that is added to the numbers.
    Can be used for turning a list of 0-based indices into a corresponding string
    of comma-separated offset indices. Offsets depend on the sign, i.e. negative
    numbers get the offset subtracted.
    Do not use offsets on lists with negative float values. 
    """
    # Increment the indices by 1. Join list-elements as strings separated by ', '.
    try:
        return ', '.join(map(str, _offsetIndices(lst, offset)))
    except TypeError:
        # If lst is not iterable (causes TypeError), then it is probably an integer.
        return lst+offset


def _insertInSetDict(dictionary, key, values):
    '''inserts values at key in dictionary containing sets. Values may be
    a single value or iterable, in which case each value is inserted'''
    if not key in dictionary:
        dictionary[key] = set()
    try:
        for v in values:
            dictionary[key].add(v)
    # Exception if values is not an iterable - insert values itself instead.
    except TypeError:
        dictionary[key].add(values)


def _insertBoundaryElement(boundaryElements, elementType, marker, nodes):
    """
    Insert an element to the boundaryElements dict.

    Parameters:

        boundaryElements  Dictionary of boundary elements

        elementType       'elm-type' according to GMSH

        marker            Boundary marker

        nodes             List of element nodes, order according to GMSH
    """
    if not marker in boundaryElements:
        boundaryElements[marker] = []
    boundaryElements[marker].append(
        {'elm-type': elementType, 'node-number-list': nodes})


def createGmshMesh(geometry, el_type=2, el_size_factor=1, dofs_per_node=1,
                   gmsh_exec_path=None, clcurv=False,
                   min_size=None, max_size=None, meshing_algorithm=None,
                   additional_options=''):

    meshGen = GmshMeshGenerator(geometry, el_type, el_size_factor, dofs_per_node,
                                gmsh_exec_path, clcurv, min_size, max_size, meshing_algorithm,
                                additional_options)

    return meshGen.create()


createMesh = createGmshMesh
create_mesh = createGmshMesh
mesh = createGmshMesh


class GmshMeshGenerator:
    '''
    Meshes geometry in GeoData objects or geo-files by calling the Gmsh executable.
    This is done when the function create() is called.
    '''

    def __init__(self, geometry, el_type=2, el_size_factor=1, dofs_per_node=1,
                 gmsh_exec_path=None, clcurv=False,
                 min_size=None, max_size=None, meshing_algorithm=None,
                 additional_options='', mesh_dir='', return_boundary_elements=False):
        '''        
        Parameters:

            geometry        GeoData instance or string containing path to .geo-file

            el_type        Integer. Element type and order. 
                           See gmsh manual for details.

            el_size_factor  Float. Factor by which the element sizes are multiplied.

            dofs_per_node    Number of degrees of freedom per node.

            gmsh_exec_path   File path to where the gmsh executable is located.

            clcurv         Set to true to make elements smaller at high curvatures. 
                           (Experimental option according to the gmsh manual)

            min_size        Minimum element size

            max_size        Maximum element size

            meshing_algorithm  String. Select mesh algorithm ('meshadapt', 'del2d',
                              'front2d',  'del3d', 'front3d', ...). 
                              See the gmsh manual for more info.

            return_boundary_elements  Flag for returning dictionary with boundary element
                                    information. Useful for applying loads on boundary.

            additional_options  String containing additional command line args for gmsh.
                               Use this if a gmsh option is not covered by the above 
                               parameters (See section 3.3 in the gmsh manual for a 
                               list of options)):
           '''
        self.geometry = geometry
        self.el_type = el_type
        self.el_size_factor = el_size_factor
        self.dofs_per_node = dofs_per_node
        self.gmsh_exec_path = gmsh_exec_path
        self.clcurv = clcurv
        self.min_size = min_size
        self.max_size = max_size
        self.meshing_algorithm = meshing_algorithm
        self.additional_options = additional_options
        self.mesh_dir = mesh_dir
        self.return_boundary_elements = return_boundary_elements

        # gmsh elements that have rectangle faces
        self._ElementsWithQuadFaces = [3, 5, 10, 12, 16, 17, 92, 93]
        self._2ndOrderElms = [8,  9, 10, 11, 12,
                              13, 14, 16, 17, 18,
                              19]
        self._2dOrderIncompleteElms = [9, 11, 13, 14,
                                       16, 17, 18, 19]
        # Apart from 16 the 2nd orders are totally untested. Only 16 (8-node quad)
        # is implemented in pycalfem though, so it does not matter.

    def create(self, is3D=False):
        '''
        Meshes a surface or volume defined by the geometry in geoData.
        Parameters:
        is3D - Optional parameter that only needs to be set if geometry
               is loaded from a geo-file, i.e. if geoData is a path string.
               Default False.

        Returns:

            coords          Node coordinates

                            [[n0_x, n0_y, n0_z],
                            [   ...           ],
                            [nn_x, nn_y, nn_z]]

            edof            Element topology

                            [[el0_dof1, ..., el0_dofn],
                            [          ...          ],
                            [eln_dof1, ..., eln_dofn]]

            dofs            Node dofs

                            [[n0_dof1, ..., n0_dofn],
                            [         ...         ],
                            [nn_dof1, ..., nn_dofn]]

            bdofs           Boundary dofs. Dictionary containing lists of dofs for
                            each boundary marker. Dictionary key = marker id.

            elementmarkers  List of integer markers. Row i contains the marker of
                            element i. Markers are similar to boundary markers and
                            can be used to identify in which region an element lies.

            boundaryElements  (optional) returned if self.return_boundary_elements is true.
                              Contains dictionary with boundary elements. The keys are markers
                              and the values are lists of elements for that marker.

        Running this function also creates object variables:

            nodesOnCurve    Dictionary containing lists of node-indices. Key is a 
                            curve-ID and the value is a list of indices of all nodes
                            on that curve, including its end points.

            nodesOnSurface  Dictionary containing lists of node-indices. Key is a
                            surface-ID and the value is a list of indices of the nodes
                            on that surface, including its boundary.

            nodesOnVolume   Dictionary containing lists of node-indices. Key is a
                            volume-ID and the value is a list of indices of the nodes
                            in that volume, including its surface. 
        '''
        # Nodes per element for different element types:
        # (taken from Chapter 9, page 89 of the gmsh manual)
        nodesPerElmDict = {1: 2,   2: 3,   3: 4,   4: 4,   5: 8,
                           6: 6,   7: 5,   8: 3,   9: 6,  10: 9,
                           11: 10, 12: 27, 13: 18, 14: 14, 15: 1,
                           16: 8,  17: 20, 18: 15, 19: 13, 20: 9,
                           21: 10, 22: 12, 23: 15, 24: 15, 25: 21,
                           26: 4,  27: 5,  28: 6,  29: 20, 30: 35,
                           31: 56, 92: 64, 93: 125}
        nodesPerElement = nodesPerElmDict[self.el_type]

        # Check for GMSH executable [NOTE]Mostly copied from trimesh2d(). TODO: Test on different systems
        gmshExe = self.gmsh_exec_path
        if gmshExe == None:
            gmshExe = None
            if sys.platform == "win32":
                gmshExe = which("gmsh.exe")
                if gmshExe == None:
                    gmshExe = which("gmsh.bat")
            else:
                gmshExe = which("gmsh")
        else:
            if not os.path.exists(gmshExe):
                gmshExe = os.path.join(
                    os.getcwd(), self.gmsh_exec_path)  # Try relative path
                if not os.path.exists(gmshExe):
                    gmshExe = None  # Relative path didnt work either

        if gmshExe == None:
            raise IOError(
                "Error: Could not find GMSH. Please make sure that the \GMSH executable is available on the search path (PATH).")
        else:
            print("Info    : GMSH -> %s" % gmshExe)

        # Create a temporary directory for GMSH

        oldStyleTempDir = False

        if self.mesh_dir != "":
            tempMeshDir = self.mesh_dir
            if not os.path.exists(tempMeshDir):
                os.mkdir(tempMeshDir)
        else:
            tempMeshDir = tempfile.mkdtemp()

        # If geometry data is given as a .geo file we will just pass it on to gmsh later.
        if type(self.geometry) is str:
            geoFilePath = self.geometry
            # In this case geoData is a path string, so the dimension must be supplied by the user.
            dim = 3 if is3D else 2
            if not os.path.exists(geoFilePath):
                geoFilePath = os.path.join(
                    os.getcwd(), geoFilePath)  # Try relative path
                if not os.path.exists(geoFilePath):
                    raise IOError(
                        "Error: Could not find geo-file " + geoFilePath)
        else:
            # Get the dimension of the model from geoData.
            dim = 3 if self.geometry.is3D else 2

            if oldStyleTempDir:
                if not os.path.exists("./gmshMeshTemp"):
                    os.mkdir("./gmshMeshTemp")
                geoFilePath = os.path.normpath(os.path.join(
                    os.getcwd(), "gmshMeshTemp/tempGeometry.geo"))  # "gmshMeshTemp/tempGeometry.geo"
            else:
                geoFilePath = os.path.normpath(
                    os.path.join(tempMeshDir, 'tempGeometry.geo'))

            self.geofile = open(geoFilePath, "w")  # Create temp geometry file
            self._writeGeoFile()  # Write geoData to file
            self.geofile.close()

        if oldStyleTempDir:
            # Filepath to the msh-file that will be generated.
            mshFileName = os.path.normpath(os.path.join(
                os.getcwd(), 'gmshMeshTemp/meshFile.msh'))
        else:
            mshFileName = os.path.normpath(
                os.path.join(tempMeshDir, 'meshFile.msh'))

        # construct options string:

        options = ""
        options += ' -' + str(dim)
        options += ' -clscale ' + str(self.el_size_factor)  # scale factor
        options += ' -o \"%s\"' % mshFileName
        options += ' -clcurv' if self.clcurv else ''
        options += ' -clmin ' + \
            str(self.min_size) if self.min_size is not None else ''
        options += ' -clmax ' + \
            str(self.max_size) if self.max_size is not None else ''
        options += ' -algo ' + self.meshing_algorithm if self.meshing_algorithm is not None else ''
        options += ' -order 2' if self.el_type in self._2ndOrderElms else ''
        options += ' -format msh22'
        options += ' ' + self.additional_options

        # Execute gmsh

        gmshExe = os.path.normpath(gmshExe)
        info("GMSH binary: "+gmshExe)
        #print('""%s" "%s" %s"' % (gmshExe, geoFilePath, options))
        #os.system('""%s" "%s" %s"' % (gmshExe, geoFilePath, options))
        #retval = os.system(r'"%s" "%s" %s' % (gmshExe, geoFilePath, options))

        output = subprocess.Popen(r'"%s" "%s" %s' % (
            gmshExe, geoFilePath, options), shell=True, stdout=subprocess.PIPE).stdout.read()

        # Read generated msh file:
        # print("Opening msh file " + mshFileName)#TEMP

        mshFile = open(mshFileName, 'r')
        info("Mesh file  : "+mshFileName)

        # print("Reading msh file...")#TEMP
        ln = mshFile.readline()
        while(ln != '$Nodes\n'):  # Read until we find the nodes
            ln = mshFile.readline()
        nbrNodes = int(mshFile.readline())
        allNodes = np.zeros([nbrNodes, dim], 'd')
        for i in range(nbrNodes):
            line = list(map(float, mshFile.readline().split()))
            # Grab the coordinates (1:3 if 2D, 1:4 if 3D)
            allNodes[i, :] = line[1:dim+1]

        while(mshFile.readline() != '$Elements\n'):  # Read until we find the elements
            pass
        # The nbr of elements (including marker elements).
        nbrElements = int(mshFile.readline())
        elements = []
        elementmarkers = []
        # temp dictionary of sets. Key:MarkerID. Value:Set. The sets will be converted to lists.
        bdofs = {}
        boundaryElements = {}
        # nodeOnPoint = {}  #dictionary pointID : nodeNumber
        self.nodesOnCurve = {}  # dictionary lineID  : set of [nodeNumber]
        self.nodesOnSurface = {}  # dictionary surfID  : set of [nodeNumber]
        self.nodesOnVolume = {}  # dictionary volID   : set of [nodeNumber]
        # Read all elements (points, surfaces, etc):
        for i in range(nbrElements):
            line = list(map(int, mshFile.readline().split()))
            eType = line[1]  # second int is the element type.
            nbrTags = line[2]  # Third int is the nbr of tags on this element.
            marker = line[3]  # Fourth int (first tag) is the marker.
            # Fifth int  is the ID of the geometric entity (points, curves, etc) that the element belongs to
            entityID = line[4]
            # The rest after tags are node indices.
            nodes = line[3+nbrTags: len(line)]

            # If the element type is the kind of element we are looking for:
            if(eType == self.el_type):
                # Add the nodes of the elements to the list.
                elements.append(nodes)
                # Add element marker. It is used for keeping track of elements (thickness, heat-production and such)
                elementmarkers.append(marker)
            else:  # If the element is not a "real" element we store its node at marker in bdof instead:
                _insertInSetDict(bdofs, marker, nodes)

                # We also store the full information as 'boundary elements'
                _insertBoundaryElement(boundaryElements, eType, marker, nodes)

            # if eType == 15: #If point. Commmented away because points only make elements if they have non-zero markers, so nodeOnPoint is not very useful.
            #    nodeOnPoint[entityID-1] = nodes[0] #insert node into nodeOnPoint. (ID-1 because we want 0-based indices)

            if eType in [1, 8, 26, 27, 28]:  # If line
                # insert nodes into nodesOnCurve
                _insertInSetDict(self.nodesOnCurve, entityID -
                                 1, _offsetIndices(nodes, -1))
            elif eType in [2, 3, 9, 10, 16, 20, 21, 22, 23, 24, 25]:  # If surfaceelement
                # insert nodes into nodesOnSurface
                _insertInSetDict(self.nodesOnSurface, entityID -
                                 1, _offsetIndices(nodes, -1))
            else:  # if volume element.
                _insertInSetDict(self.nodesOnVolume, entityID -
                                 1, _offsetIndices(nodes, -1))

        elements = np.array(elements)
        for key in bdofs.keys():  # Convert the sets of boundary nodes to lists.
            bdofs[key] = list(bdofs[key])
        for key in self.nodesOnCurve.keys():  # Convert set to list
            self.nodesOnCurve[key] = list(self.nodesOnCurve[key])
        for key in self.nodesOnSurface.keys():  # Convert set to list
            self.nodesOnSurface[key] = list(self.nodesOnSurface[key])
        for key in self.nodesOnVolume.keys():  # Convert set to list
            self.nodesOnVolume[key] = list(self.nodesOnVolume[key])

        mshFile.close()

        # Remove temporary mesh directory if not explicetly specified.

        if self.mesh_dir == "":
            shutil.rmtree(tempMeshDir)

        dofs = createdofs(np.size(allNodes, 0), self.dofs_per_node)

        if self.dofs_per_node > 1:  # This if-chunk copied from pycalfem_utils.py
            self.topo = elements
            expandedElements = np.zeros(
                (np.size(elements, 0), nodesPerElement*self.dofs_per_node), 'i')
            elIdx = 0
            for elementTopo in elements:
                for i in range(nodesPerElement):
                    expandedElements[elIdx, i*self.dofs_per_node:(
                        i*self.dofs_per_node+self.dofs_per_node)] = dofs[elementTopo[i]-1, :]
                elIdx += 1

            for keyID in bdofs.keys():
                bVerts = bdofs[keyID]
                bVertsNew = []
                for i in range(len(bVerts)):
                    for j in range(self.dofs_per_node):
                        bVertsNew.append(dofs[bVerts[i]-1][j])
                bdofs[keyID] = bVertsNew

            if self.return_boundary_elements:
                return allNodes, np.asarray(expandedElements), dofs, bdofs, elementmarkers, boundaryElements
            return allNodes, np.asarray(expandedElements), dofs, bdofs, elementmarkers

        if self.return_boundary_elements:
            return allNodes, elements, dofs, bdofs, elementmarkers, boundaryElements
        return allNodes, elements, dofs, bdofs, elementmarkers

    def _writeGeoFile(self):
        # key is marker, value is a list of point indices (0-based) with that marker
        pointMarkers = {}
        curveMarkers = {}
        surfaceMarkers = {}
        volumeMarkers = {}

        # WRITE POINTS:
        for ID, [coords, elSize, marker] in self.geometry.points.items():
            self.geofile.write("Point(%i) = {%s};\n" % (
                ID+1, _formatList(coords + [elSize])))
            _insertInSetDict(pointMarkers, marker, ID)

        # WRITE CURVES:
        for ID, [curveName, points, marker, elOnCurve, distributionString, distributionVal] in self.geometry.curves.items():
            self.geofile.write("%s(%i) = {%s};\n" % (
                curveName, ID+1, _formatList(points, 1)))

            # Transfinite Line{2} = 20 Using Bump 0.05;
            if elOnCurve != None:
                distribution = "" if distributionString == None else "Using %s %f" % (
                    distributionString, distributionVal)
                self.geofile.write("Transfinite Line{%i} = %i %s;\n" % (
                    ID+1, elOnCurve+1, distribution))
                # +1 on elOnCurve because gmsh actually takes the number of nodes on the curve, not elements on the curve.
            _insertInSetDict(curveMarkers, marker, ID)

        # WRITE SURFACES:
        for ID, [surfName, outerLoop, holes, ID, marker, isStructured] in self.geometry.surfaces.items():
            # First we write line loops for the surface edge and holes (if there are any holes):
            self._writeLineLoop(outerLoop, ID+1)
            holeIDs = []
            for hole, i in zip(holes, range(len(holes))):
                # Create a hopefully unique ID-number for the line loop: Like 10015 or 1540035
                # (If gmsh uses 32-bit integers for IDs then IDs over 214'748 will break)
                holeID = 10000 * (ID+1) + 10 * i + 5
                self._writeLineLoop(hole, holeID)
                holeIDs.append(holeID)
            # Second, we write the surface itself:
            # If we have hole we want to include them in the surface.
            holeString = "" if not holeIDs else ", " + _formatList(holeIDs)
            # Like "Plane Surface(2) = {4, 2, 6, 8}
            self.geofile.write("%s(%i) = {%s%s};\n" % (
                surfName, ID+1, ID+1, holeString))
            # Lastly, we make the surface transfinite if it is a structured surface:
            if isStructured:
                cornerPoints = set()
                # Find the corner points. This is possibly unnecessary since Gmsh can do this automatically.
                for c in outerLoop:
                    curvePoints = self.geometry.curves[c][1]
                    cornerPoints.add(curvePoints[0])
                    cornerPoints.add(curvePoints[-1])
                cornerPoints = list(cornerPoints)
                self.geofile.write("Transfinite Surface{%i} = {%s};\n" % (
                    ID+1, _formatList(cornerPoints, 1)))  # Like Transfinite Surface{1} = {1,2,3,4};
                # Transfinite Surface has an optional argument (about triangle orientation) that is not implemented here.
            _insertInSetDict(surfaceMarkers, marker, ID)

        # WRITE VOLUMES:
        for ID, [outerLoop, holes, ID, marker, isStructured] in self.geometry.volumes.items():
            # Surface loops for the volume boundary and holes (if any):
            self._writeSurfaceLoop(outerLoop, ID+1)
            holeIDs = []
            for hole, i in zip(holes, range(len(holes))):
                # ID-number for the hole surface loop
                holeID = 10000 * (ID+1) + 10 * i + 7
                self._writeSurfaceLoop(hole, holeID)
                holeIDs.append(holeID)
            # Write the volume itself:
            # If we have hole we want to include them in the surface.
            holeString = "" if not holeIDs else " , " + _formatList(holeIDs)
            # Like "Plane Surface(2) = {4, 2, 6, 8}
            self.geofile.write(
                "Volume(%i) = {%s%s};\n" % (ID+1, ID+1, holeString))
            # Lastly, we make the volume transfinite if it is a structured volume:
            if isStructured:
                self.geofile.write("Transfinite Volume{%i} = {};\n" % (ID+1))
                # We don't find the corner points of the structured volume like we did with the surfaces. Gmsh can actually
                # find the corners automatically.
            _insertInSetDict(volumeMarkers, marker, ID)

        # MAYBE MAKE QUADS:
        if(self.el_type in self._ElementsWithQuadFaces):  # If we have quads surfaces on the elements
            self.geofile.write("Mesh.RecombineAll = 1;\n")

        # WRITE POINT MARKERS:
        for marker, IDlist in pointMarkers.items():
            if marker != 0:
                self.geofile.write("Physical Point(%i) = {%s};\n" % (
                    marker, _formatList(IDlist, 1)))

        # WRITE CURVE MARKERS:
        for marker, IDlist in curveMarkers.items():
            self.geofile.write("Physical Line(%i) = {%s};\n" % (
                marker, _formatList(IDlist, 1)))

        # WRITE SURFACE MARKERS:
        for marker, IDlist in surfaceMarkers.items():
            self.geofile.write("Physical Surface(%i) = {%s};\n" % (
                marker, _formatList(IDlist, 1)))

        # WRITE SURFACE MARKERS:
        for marker, IDlist in volumeMarkers.items():
            self.geofile.write("Physical Volume(%i) = {%s};\n" % (
                marker, _formatList(IDlist, 1)))

        # If the element type is of an incomplete second order type
        # (i.e it is an 2nd order element without nodes in the middle of the element face),
        # then we need to specify this in the geo-file:
        if self.el_type in self._2dOrderIncompleteElms:
            self.geofile.write("Mesh.SecondOrderIncomplete=1;\n")

    def _writeLineLoop(self, lineIndices, loopID):
        # endPoints is used to keep track of at which points the curves start and end (i.e the direction of the curves)
        endPoints = []
        # lineIndices is a list of curve indices (0-based here, but 1-based later in the method)
        for i in lineIndices:
            curvePoints = self.geometry.curves[i][1]
            endPoints.append([curvePoints[0], curvePoints[-1]])

        # We need the indices to be 1-based rather than 0-based in the next loop. (Some indices will be preceded by a minus-sign)
        lineIndices = _offsetIndices(lineIndices, 1)
        isFirstLine = True
        nbrLinesinLoop = len(lineIndices)
        # In this loop we reverse the direction of some lines in the LineLoop to make them conform to the format that Gmsh expects.
        for k in range(nbrLinesinLoop):
            if isFirstLine and nbrLinesinLoop > 1:
                isFirstLine = False
                # If last point of the first line exists in the endpoints of the second line... Do nothing
                if endPoints[0][1] in endPoints[1]:
                    pass
                # Else if the first point in the first line exists in the endpoints of the second line:
                elif endPoints[0][0] in endPoints[1]:
                    endPoints[0].reverse()
                    lineIndices[0] *= -1  # Reverse the direction of the line
                else:
                    raise Exception(
                        "ERROR: The first curve of line-loop %i does not link up to the subsequent curve" % loopID)
            elif endPoints[k][0] == endPoints[k-1][1]:
                pass
            elif endPoints[k][1] == endPoints[k-1][1]:
                endPoints[k].reverse()
                lineIndices[k] *= -1  # Reverse the direction of the line
            else:
                raise Exception(
                    "ERROR: The %i th curve (starting from 0) of a line-loop %i does not link up with the preceding curve" % (k, loopID))
            if k == nbrLinesinLoop-1 and endPoints[k][1] != endPoints[0][0]:
                # If last line AND the last point of the last curve not equal the first point of the first curve:
                raise Exception(
                    "ERROR: The last curve of a line-loop %i does not join up with the first curve" % loopID)

        # If the model is in 2D we need to make all line loops counter-clockwise so surface normals point in the positive z-direction.
        if not self.geometry.is3D:
            lineIndices = self._makeCounterClockwise(lineIndices)

        self.geofile.write("Line Loop(%i) = {%s};\n" % (loopID, _formatList(
            lineIndices)))  # (lineIndices are alreay 1-based here)

    def _makeCounterClockwise(self, lineIndices):
        '''If the lineIndices describe a line loop that is not counterclockwise,
        this function will return a counterclockwise version of lineIndices
        (i.e. all indices multiplied by -1).
        lineIndices is a list of integers (1-based line indices) that may be negative, but not 0'''
        # Method described at http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        summa = 0.0  # Counter-clockwise if the sum ends up negative.
        for index in lineIndices:
            sign = -1 if index < 0 else 1
            # Make a copy of the line index that is positive and 0-based.
            realIndex = sign*index - 1
            curveType = self.geometry.curves[realIndex][0]
            pointIDs = self.geometry.curves[realIndex][1]
            if curveType in ['Spline', 'BSpline']:
                points = [self.geometry.points[ID][0]
                          for ID in pointIDs]  # [[x,y,z], [x,y,z], ...]
                # Reverse the order of the points if the curve direction is reversed.
                points = points if sign == 1 else points[::-1]
                # For every point along the curve except the last:
                for i in range(len(pointIDs)-1):
                    # (x2-x1)(y2+y1).
                    summa += (points[i+1][0] - points[i][0]) * \
                        (points[i+1][1] + points[i][1])
            elif curveType == 'Circle':
                # We will find a point 'd' on the middle of the circle arc, and use a-d-c as approximation of the arc.
                # 3-by-3 array. The rows are start-center-end points and columns are x,y,z.
                points = np.array([self.geometry.points[ID][0]
                                   for ID in pointIDs])
                # Reverse the order of the points if the curve direction is reversed.
                points = points if sign == 1 else points[::-1]
                a = points[0, :]  # start
                b = points[1, :]  # center
                c = points[2, :]  # end
                r = np.linalg.norm(a-b)  # radius
                d = b + r * (a + 2*b + c) / np.linalg.norm(a + 2*b + c)
                approxArc = np.vstack((a, d, c))
                for i in range(len(approxArc)-1):
                    # (x2-x1)(y2+y1).
                    summa += (approxArc[i+1, 0] - approxArc[i, 0]) * \
                        (approxArc[i+1, 1] + approxArc[i, 1])
            elif curveType == 'Ellipse':
                # We will find a point 'd' near the middle of the circle arc, and use a-d-c as approximation of the arc.
                # The only difference from the circle above, is that the radius at d is approximated as the mean distance between
                # the center and the two other points.
                # 4-by-3 array. The rows are start-center-majAxis-end points and columns are x,y,z.
                points = np.array([self.geometry.points[ID][0]
                                   for ID in pointIDs])
                # skip the major axis point (row 2)
                points = points[[0, 1, 3], :]
                # Reverse the order of the points if the curve direction is reversed.
                points = points if sign == 1 else points[::-1]
                a = points[0, :]  # start
                b = points[1, :]  # center
                c = points[2, :]  # end
                r = (np.linalg.norm(a-b) + np.linalg.norm(c-b)) / \
                    2  # approximate radius
                d = b + r * (a + 2*b + c) / np.linalg.norm(a + 2*b + c)
                approxArc = np.vstack((a, d, c))
                for i in range(len(approxArc)-1):
                    # (x2-x1)(y2+y1).
                    summa += (approxArc[i+1, 0] - approxArc[i, 0]) * \
                        (approxArc[i+1, 1] + approxArc[i, 1])
        # If the sum is positive the loop (closed polygon) is clockwise, so reverse the direction of all curves:
        if summa > 0:
            lineIndices = [-x for x in lineIndices]
        return lineIndices

    def _writeSurfaceLoop(self, outerLoop, ID):
        self.geofile.write("Surface Loop(%i) = {%s};\n" % (
            ID, _formatList(outerLoop, 1)))

    # --- Compatibility properties

    @property
    def elType(self):
        return self.el_type

    @elType.setter
    def elType(self, value):
        self.el_type = value

    @property
    def elSizeFactor(self):
        return self.el_size_factor

    @elSizeFactor.setter
    def elSizeFactor(self, value):
        self.el_size_factor = value

    @property
    def dofsPerNode(self):
        return self.dofs_per_node

    @dofsPerNode.setter
    def dofsPerNode(self, value):
        self.dofs_per_node = value

    @property
    def gmshExecPath(self):
        return self.gmsh_exec_path

    @gmshExecPath.setter
    def gmshExecPath(self, value):
        self.gmsh_exec_path = value

    @property
    def minSize(self):
        return self.min_size

    @minSize.setter
    def minSize(self, value):
        self.min_size = value

    @property
    def maxSize(self):
        return self.max_size

    @maxSize.setter
    def maxSize(self, value):
        self.max_size = value

    @property
    def meshingAlgorithm(self):
        return self.meshing_algorithm

    @meshingAlgorithm.setter
    def meshingAlgorithm(self, value):
        self.meshing_algorithm = value

    @property
    def additionalOptions(self):
        return self.additional_options

    @additionalOptions.setter
    def additionalOptions(self, value):
        self.additional_options = value

    @property
    def meshDir(self):
        return self.mesh_dir

    @meshDir.setter
    def meshDir(self, value):
        self.mesh_dir = value

    @property
    def returnBoundaryElements(self):
        return self.return_boundary_elements

    @returnBoundaryElements.setter
    def returnBoundaryElements(self, value):
        self.return_boundary_elements = value


GmshMesh = GmshMeshGenerator


def trimesh2d(vertices, segments=None, holes=None, maxArea=None, quality=True, dofs_per_node=1, logFilename="tri.log", triangleExecutablePath=None):
    """
    Triangulates an area described by a number vertices (vertices) and a set
    of segments that describes a closed polygon. 

    Parameters:

        vertices            array [nVertices x 2] with vertices describing the geometry.

                            [[v0_x, v0_y],
                             [   ...    ],
                             [vn_x, vn_y]]

        segments            array [nSegments x 3] with segments describing the geometry.

                            [[s0_v0, s0_v1,marker],
                             [        ...        ],
                             [sn_v0, sn_v1,marker]]

        holes               [Not currently used]

        maxArea             Maximum area for triangle. (None)

        quality             If true, triangles are prevented having angles < 30 degrees. (True)

        dofs_per_node         Number of degrees of freedom per node.

        logFilename         Filename for triangle output ("tri.log")

    Returns:

        coords              Node coordinates

                            [[n0_x, n0_y],
                             [   ...    ],
                             [nn_x, nn_y]]

        edof                Element topology

                            [[el0_dof1, ..., el0_dofn],
                             [          ...          ],
                             [eln_dof1, ..., eln_dofn]]

        dofs                Node dofs

                            [[n0_dof1, ..., n0_dofn],
                             [         ...         ],
                             [nn_dof1, ..., nn_dofn]]

        bdofs               Boundary dofs. Dictionary containing lists of dofs for
                            each boundary marker. Dictionary key = marker id.

    """

    # Check for triangle executable

    triangleExecutable = triangleExecutablePath

    if triangleExecutable == None:
        triangleExecutable = ""
        if sys.platform == "win32":
            triangleExecutable = which("triangle.exe")
        else:
            triangleExecutable = which("triangle")
    else:
        if not os.path.exists(triangleExecutable):
            triangleExecutable = None

    if triangleExecutable == None:
        error("Error: Could not find triangle. Please make sure that the \ntriangle executable is available on the search path (PATH).")
        return None, None, None, None

    # Create triangle options

    options = ""

    if maxArea != None:
        options += "-a%f " % maxArea + " "
    if quality:
        options += "-q"

    # Set initial variables

    nSegments = 0
    nHoles = 0
    nAttribs = 0
    nBoundaryMarkers = 1
    nVertices = len(vertices)

    # All files are created as temporary files

    if not os.path.exists("./trimesh.temp"):
        os.mkdir("./trimesh.temp")

    filename = "./trimesh.temp/polyfile.poly"

    if not segments is None:
        nSegments = len(segments)

    if not holes is None:
        nHoles = len(holes)

    # Create a .poly file

    polyFile = open(filename, "w")
    polyFile.write("%d 2 %d \n" % (nVertices, nAttribs))

    i = 0

    for vertex in vertices:
        polyFile.write("%d %g %g\n" % (i, vertex[0], vertex[1]))
        i = i + 1

    polyFile.write("%d %d \n" % (nSegments, nBoundaryMarkers))

    i = 0

    for segment in segments:
        polyFile.write("%d %d %d %d\n" %
                       (i, segment[0], segment[1], segment[2]))
        i = i + 1

    polyFile.write("0\n")

    polyFile.close()

    # Execute triangle

    os.system("%s %s %s > tri.log" % (triangleExecutable, options, filename))

    # Read results from triangle

    strippedName = os.path.splitext(filename)[0]

    nodeFilename = "%s.1.node" % strippedName
    elementFilename = "%s.1.ele" % strippedName
    polyFilename = "%s.1.poly" % strippedName

    # Read vertices

    allVertices = None
    boundaryVertices = {}

    if os.path.exists(nodeFilename):
        nodeFile = open(nodeFilename, "r")
        nodeInfo = list(map(int, nodeFile.readline().split()))

        nNodes = nodeInfo[0]

        allVertices = np.zeros([nNodes, 2], 'd')

        for i in range(nNodes):
            vertexRow = list(map(float, nodeFile.readline().split()))

            boundaryMarker = int(vertexRow[3])

            if not (boundaryMarker in boundaryVertices):
                boundaryVertices[boundaryMarker] = []

            allVertices[i, :] = [vertexRow[1], vertexRow[2]]
            boundaryVertices[boundaryMarker].append(i+1)

        nodeFile.close()

    # Read elements

    elements = []

    if os.path.exists(elementFilename):
        elementFile = open(elementFilename, "r")
        elementInfo = list(map(int, elementFile.readline().split()))

        nElements = elementInfo[0]

        elements = np.zeros([nElements, 3], 'i')

        for i in range(nElements):
            elementRow = list(map(int, elementFile.readline().split()))
            elements[i, :] = [elementRow[1]+1,
                              elementRow[2]+1, elementRow[3]+1]

        elementFile.close()

    # Clean up

    try:
        pass
        # os.remove(filename)
        # os.remove(nodeFilename)
        # os.remove(elementFilename)
        # os.remove(polyFilename)
    except:
        pass

    # Add dofs in edof and bcVerts

    dofs = cfc.createdofs(np.size(allVertices, 0), dofs_per_node)

    if dofs_per_node > 1:
        expandedElements = np.zeros(
            (np.size(elements, 0), 3*dofs_per_node), 'i')
        dofs = cfc.createdofs(np.size(allVertices, 0), dofs_per_node)

        elIdx = 0

        for elementTopo in elements:
            for i in range(3):
                expandedElements[elIdx, i*dofs_per_node:(
                    i*dofs_per_node+dofs_per_node)] = dofs[elementTopo[i]-1, :]
            elIdx += 1

        for bVertIdx in boundaryVertices.keys():
            bVert = boundaryVertices[bVertIdx]
            bVertNew = []
            for i in range(len(bVert)):
                for j in range(dofs_per_node):
                    bVertNew.append(dofs[bVert[i]-1][j])

            boundaryVertices[bVertIdx] = bVertNew

        return allVertices, np.asarray(expandedElements), dofs, boundaryVertices

    return allVertices, elements, dofs, boundaryVertices
