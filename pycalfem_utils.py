#!/bin/env python
# -*- coding: utf-8 -*-

import os, sys

from numpy import *
from pycalfem import *

def readInt(f):
	"""
	Läs in heltal från en rad i en fil.
	Returnera en lista av värden
	"""
	return map(int, f.readline().split())
	
def readFloat(f):
	"""
	Läs in flyttal från en rad i en fil.
	Returnera lista av värden
	"""
	return map(float, f.readline().split())
	
def readSingleInt(f):
	"""
	Läser enstaka integer från fil.
	"""
	return readInt(f)[0] 

def readSingleFloat(f):
	"""
	Läser enstaka float från fil.
	"""
	return readFloat(f)[0]
	
def which(filename):
	if not os.environ.has_key('PATH') or os.environ['PATH'] == '':
		p = os.defpath
	else:
		p = os.environ['PATH']
	
	pathlist = p.split (os.pathsep)
	pathlist.append(".")
	
	for path in pathlist:
		f = os.path.join(path, filename)
		if os.access(f, os.X_OK):
			return f
	return None

	
def trimesh2d(vertices, segments = None, holes = None, maxArea=None, quality=True, dofsPerNode=1, logFilename="tri.log"):
	"""
	Triangulerar ett område beskrivet av ett antal punkter (vertices) samt 
	ett antal linjer som utgör en slute polygon (segments), maxArea anger 
	maximal tillåten storlek på de genererade trianglarna. Rutinen
	anropar triangle för meshgenerering. Temporära filer skapas i 
	användarens temp-katalog.
	"""
	
	# Check for triangle executable
	
	triangleExecutable = which("triangle")
	if triangleExecutable==None:
		print "Error: Could not find triangle. Please make sure that the \ntriangle executable is available on the search path (PATH)."
		return None, None, None, None
	
	# Create triangle options
	
	options = ""
	
	if maxArea!=None:
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
	
	if not os.path.exists("./triangle"):
		os.mkdir("./triangle")
		
	filename = "./triangle/polyfile.poly"
	
	if segments!=None:
		nSegments = len(segments)
	
	if holes!=None:
		nHoles = len(holes)
	
	# Create a .poly file
	
	polyFile = file(filename, "w")
	polyFile.write("%d 2 %d \n" % (nVertices, nAttribs))
	
	i = 0
	
	for vertex in vertices:
		polyFile.write("%d %g %g\n" % (i, vertex[0], vertex[1])) 
		i = i + 1
		
	polyFile.write("%d %d \n" % (nSegments, nBoundaryMarkers))
		
	i = 0
		
	for segment in segments:
		polyFile.write("%d %d %d %d\n" % (i, segment[0], segment[1], segment[2]))
		i = i + 1
		
	polyFile.write("0\n")
	
	polyFile.close()

	# Execute triangle
	
	os.system("%s %s %s > tri.log" % ("triangle", options, filename))
	
	# Read results from triangle
	
	strippedName = os.path.splitext(filename)[0]
	
	nodeFilename = "%s.1.node" % strippedName
	elementFilename = "%s.1.ele" % strippedName
	polyFilename = "%s.1.poly" % strippedName
	
	# Read vertices
	
	allVertices = None
	boundaryVertices = {}
	
	if os.path.exists(nodeFilename):
		nodeFile = file(nodeFilename, "r")
		nodeInfo = map(int, nodeFile.readline().split())
		
		nNodes = nodeInfo[0]
		
		allVertices = zeros([nNodes,2], 'd')
		
		for i in range(nNodes):
			vertexRow = map(float, nodeFile.readline().split())
			
			boundaryMarker = int(vertexRow[3])
			
			if not boundaryVertices.has_key(boundaryMarker):
				boundaryVertices[boundaryMarker] = []
			
			allVertices[i,:] = [vertexRow[1], vertexRow[2]]
			boundaryVertices[boundaryMarker].append(i+1)
			
		nodeFile.close()
		
	# Read elements
			
	elements = []
		
	if os.path.exists(elementFilename):
		elementFile = file(elementFilename, "r")
		elementInfo = map(int, elementFile.readline().split())
		
		nElements = elementInfo[0]
		
		elements = zeros([nElements,3],'i')
		
		for i in range(nElements):
			elementRow = map(int, elementFile.readline().split())
			elements[i,:] = [elementRow[1]+1, elementRow[2]+1, elementRow[3]+1]
			
		elementFile.close()
			
	# Clean up
	
	try:
		os.remove(filename)
		os.remove(nodeFilename)
		os.remove(elementFilename)
		os.remove(polyFilename)
	except:
		pass
	
	# Add dofs in edof and bcVerts
	
	dofs = createdofs(size(allVertices,0),dofsPerNode)
	
	if dofsPerNode>1:
		expandedElements = zeros((size(elements,0),3*dofsPerNode),'i')
		dofs = createdofs(size(allVertices,0),dofsPerNode)
		
		elIdx = 0
		
		for elementTopo in elements:		
			for i in range(3):
				expandedElements[elIdx,i*dofsPerNode:(i*dofsPerNode+dofsPerNode)] = dofs[elementTopo[i]-1,:]
			elIdx += 1
			
		for bVertIdx in boundaryVertices.keys():
			bVert = boundaryVertices[bVertIdx]
			bVertNew = []
			for i in range(len(bVert)):
				for j in range(dofsPerNode):
					bVertNew.append(dofs[bVert[i]-1][j])
					
			boundaryVertices[bVertIdx] = bVertNew
			
		return allVertices, expandedElements, dofs, boundaryVertices
		
	
	return allVertices, elements, dofs, boundaryVertices