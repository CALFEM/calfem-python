#TODO: Consider adding more checks whether entities exist.

class GeoData:
    '''
    Instances of GeoData can hold geometric data and be passed to 
    GmshMesher in pycalfem_Mesh to mesh the geometry.
    '''

    def __init__(self):
        self.points = {}    #dict of [[x, y, z], elSize, marker]
        self.curves = {}    #dict of [curvTypestring, [p1, p2, ... pn], marker, elementsOnCurve, distributionString, distributionVal]
        self.surfaces = {}  #dict of [SurfaceTypeString, [c1, c2 ... cn], [[c1, c2 ... cm], ... [c1, ... ck]], ID, marker, isStructured]. c means curve-ID.
        self.volumes = {}   #dict of [[s1, s2 ..], [[s1,s2...],[s1,s2..],..], ID, marker, isStructured]    
        self.is3D = False   #This is automatically set to True if a 3D point is added.
        self._pointIDspecified = False
        self._nextPointID = 0
        self._curveIDspecified = False
        self._nextcurveID = 0
        self._surfaceIDspecified = False
        self._nextsurfaceID = 0
        self._volumeIDspecified = False
        self._nextvolumeID = 0
    
    
    def removePoint(self, ID):
        '''Removes the point with this ID'''
        self.points.pop(ID)
    
    
    def removeCurve(self, ID):
        '''Removes the curve with this ID'''
        self.curve.pop(ID)
    
    
    def removeSurface(self, ID):
        '''Removes the surface with this ID'''
        self.surfaces.pop(ID)
    
    
    def removeVolume(self, ID):
        '''Removes the volume with this ID'''
        self.volumes.pop(ID)
    
    
    def getPointCoords(self, IDs=None):
        '''
        Returns an N-by-3 list of point coordinates if the parameter is
        a list of IDs. If the parameter is just a single integer then 
        a single coordinate (simple 3-element list) is returned.
        If the parameter is undefined (or None) all point coords will be returned
        '''
        if IDs == None:
            return [p[0] for p in self.points.values()]
        try:
            pointCoords = [self.points[pID][0] for pID in IDs]
        except TypeError: #IDs was not iterable. Probably just a single ID.
            pointCoords = self.points[IDs][0]
        return pointCoords
    
    
    def pointsOnCurves(self, IDs):
        '''
        Returns a list of all geometric points (not nodes) on the curves
        specified in IDs. IDs may be an integer or a list of integers.
        '''
        return self._subentitiesOnEntities(IDs, self.curves, 1)
    
    
    def stuffOnSurfaces(self, IDs):
        '''
        Returns lists of all geometric points and curves on the surfaces
        specified in IDs. IDs may be an integer or a list of integers
        '''
        curveSet = self._subentitiesOnEntities(IDs, self.surfaces, 1) #Curves on the outer edges
        curveSet.update( self._subentityHolesOnEntities(IDs, self.surfaces, 2) ) #Curves on the holes
        pointList = self.pointsOnCurves(curveSet) #Points on the curves of these surfaces.
        return pointList, list(curveSet)
        
    def stuffOnVolumes(self, IDs):
        '''
        Returns lists of all geometric points, curves, and surfaces on the volumes
        specified in IDs. IDs may be an integer or a list of integers
        '''
        surfaceSet = self._subentitiesOnEntities(IDs, self.surfaces, 0)
        surfaceSet.update( self._subentitiesOnEntities(IDs, self.surfaces, 1) )
        pointList, curveList = self.stuffOnSurfaces(surfaceSet)
        return pointList, curveList, list(surfaceSet)
    
    def _subentitiesOnEntities(self, IDs, entityDict, index):
        '''
        Duplicate code. Gets the IDs of the subentities that
        make up an entity, i.e. the points that define a curve or
        the curves that define a surface. Note that only the outer
        subentities of surfaces and volumes can be extracted with
        this function. For holes use _subentityHolesOnEntities().
        '''
        theSet = set()
        try:
            for ID in IDs:
                theSet.update(entityDict[ID][index])
        except TypeError: #IDs is not iterable, so it is probably a single ID
            theSet.update(entityDict[IDs][index])
        return theSet
    
    def _subentityHolesOnEntities(self, IDs, entityDict, index):
        '''Duplicate code. Does the same thing as _subentitiesOnEntities(), but for holes'''
        theSet = set()
        try:
            for ID in IDs:
                for hole in entityDict[ID][index]:
                    theSet.update(hole)
        except TypeError: #IDs is not iterable, so it is probably a single ID
            for hole in entityDict[IDs][index]:
                theSet.update(hole)
        return theSet
        
    
    def addPoint(self, coord, ID=None, marker=0, elSize=1):
        '''
        Adds a point.
        
        Parameters:
        coord     - [x, y] or [x, y, z].
                    List, not array.
        
        ID        - Positive integer ID of this point. If left unspecified the
                    point will be assigned the smallest unused point-ID.
                    It is recommended to specify all point-IDs or none.
                    
        marker    - Marker applied to this point. Default 0.
                    It is not a good idea to apply non-zero markers to points
                    that are control points on B-splines or center points on 
                    circles/ellipses, since this can lead to "loose" nodes
                    that are not part of any elements.
        
        elSize    - The size of elements at this point. Default 1. Use to make
                    a mesh denser or sparser here. Only affects unstructured
                    meshes
        '''
        if len(coord)==3: #A 3D point is inserted.
            self.is3D = True 
        else: #The point is in 2D (we assume)
            coord = coord+[0] #Pad with a 0. (we store points as 3D points for consistency's sake)
            
        if ID==None: #ID is not specified. Get an ID for this point:
            ID = self._getNewPointID()
        else:
            self._pointIDspecified = True
            
        self.points[ID] = [coord, elSize, marker]
    
    
    def addSpline(self, points, ID=None, marker=0, elOnCurve=None, elDistribType=None, elDistribVal=None):
        '''
        Adds a Spline curve
        
        points    - List of indices of control points that make a Spline
                    [p1, p2, ... , pn]
                    
        ID        - Positive integer ID of this curve. If left unspecified the
                    curve will be assigned the smallest unused curve-ID.
                    It is recommended to specify all curve-IDs or none.
                    
        marker    - Integer. Marker applied to this curve. Default 0.
                    
        elOnCurv  - Positive integer. Elements on curve. 
                    The number of element edges that will be distributed
                    along this curve. Only works for structured meshes. 
                    
        elDistribType -
                    String. Either "bump" or "progression". 
                    Determines how the density of elements vary along the curve
                    for structured meshes. Only works for structured meshes.
                    elOnCurv and elDistribVal must be be defined if this param
                    is used.
                                        
        elDistribVal -
                    Float. Determines how severe the element distribution is.
                    Only works for structured meshes. elOnCurv and 
                    elDistribType must be be defined if this param is used.
                    
                        bump:
                    Smaller value means elements are bunched up at the edges
                    of the curve, larger means bunched in the middle.
                    
                        progression:
                    The edge of each element along this curve (from starting
                    point to end) will be larger than the preceding one by 
                    this factor.
                    elDistribVal = 2 meaning for example that each line element 
                    in the series will be twice as long as the preceding one.
                    elDistribVal < 1 makes each element smaller than the 
                    preceeding one.
        '''
        self._addCurve("Spline", points, ID, marker, elOnCurve, elDistribType, elDistribVal)
        
        
    def addBSpline(self, points, ID=None, marker=0, elOnCurve=None,  elDistribType=None, elDistribVal=None):
        '''
        Adds a B-Spline curve
        
        points    - List of indices of control points that make a B-spline
                    [p1, p2, ... , pn]
                    
        ID        - Positive integer ID of this curve. If left unspecified the
                    curve will be assigned the smallest unused curve-ID.
                    It is recommended to specify all curve-IDs or none.
                    
        marker    - Integer. Marker applied to this curve. Default 0.
                    
        elOnCurv  - Positive integer. Elements on curve. 
                    The number of element edges that will be distributed
                    along this curve. Only works for structured meshes. 
                    
        elDistribType -
                    String. Either "bump" or "progression". 
                    Determines how the density of elements vary along the curve
                    for structured meshes. Only works for structured meshes.
                    elOnCurv and elDistribVal must be be defined if this param
                    is used.
                                        
        elDistribVal -
                    Float. Determines how severe the element distribution is.
                    Only works for structured meshes. elOnCurv and 
                    elDistribType must be be defined if this param is used.
                    
                        bump:
                    Smaller value means elements are bunched up at the edges
                    of the curve, larger means bunched in the middle.
                    
                        progression:
                    The edge of each element along this curve (from starting
                    point to end) will be larger than the preceding one by 
                    this factor.
                    elDistribVal = 2 meaning for example that each line element 
                    in the series will be twice as long as the preceding one.
                    elDistribVal < 1 makes each element smaller than the 
                    preceeding one.
        '''
        self._addCurve("BSpline", points, ID, marker, elOnCurve,  elDistribType, elDistribVal)


    def addCircle(self, points, ID=None, marker=0, elOnCurve=None,  elDistribType=None, elDistribVal=None):
        '''
        Adds a Circle arc curve.
        
        points    - list of 3 indices of point that make a circle arc smaller
                    than Pi.
                    [startpoint, centerpoint, endpoint]
        
        ID        - Positive integer ID of this curve. If left unspecified the
                    curve will be assigned the smallest unused curve-ID.
                    It is recommended to specify all curve-IDs or none.
                    
        marker    - Marker applied to this curve. Default 0.
                    
        elOnCurv  - Elements on curve.
                    The number of element edges that will be distributed
                    along this curve. Only works for structured meshes.
                    
        elDistribType -
                    String. Either "bump" or "progression". 
                    Determines how the density of elements vary along the curve
                    for structured meshes. Only works for structured meshes.
                    elOnCurv and elDistribVal must be be defined if this param
                    is used.
                                        
        elDistribVal -
                    Float. Determines how severe the element distribution is.
                    Only works for structured meshes. elOnCurv and 
                    elDistribType must be be defined if this param is used.
                    
                        bump:
                    Smaller value means elements are bunched up at the edges
                    of the curve, larger means bunched in the middle.
                    
                        progression:
                    The edge of each element along this curve (from starting
                    point to end) will be larger than the preceding one by 
                    this factor.
                    elDistribVal = 2 meaning for example that each line element 
                    in the series will be twice as long as the preceding one.
                    elDistribVal < 1 makes each element smaller than the 
                    preceeding one.
        '''
        if len(points) != 3:
            raise IndexError("Circle: points must be a list of 3 positive integers denoting point indices")
        self._addCurve("Circle", points, ID, marker, elOnCurve,  elDistribType, elDistribVal)
        
        
    def addEllipse(self, points, ID=None, marker=0, elOnCurve=None, elDistribType=None, elDistribVal=None):
        '''
        Adds a Ellipse arc curve.
        
        points    - List of 4 indices of point that make a ellipse arc smaller
                    than Pi.
                    [startpoint, centerpoint, mAxisPoint, endpoint]
                    Startpoint is the starting point of the arc.
                    Centerpoint is the point at the center of the ellipse.
                    MAxisPoint is any point on the major axis of the ellipse.
                    Endpoint is the end point of the arc.
        
        ID        - Positive integer ID of this curve. If left unspecified the
                    curve will be assigned the smallest unused curve-ID.
                    It is recommended to specify all curve-IDs or none.
                    
        marker    - Integer. Marker applied to this curve. Default 0.
                    
        elOnCurv  - Positive integer. Elements on curve. 
                    The number of element edges that will be distributed
                    along this curve. Only works for structured meshes. 
                    
        elDistribType -
                    String. Either "bump" or "progression". 
                    Determines how the density of elements vary along the curve
                    for structured meshes. Only works for structured meshes.
                    elOnCurv and elDistribVal must be be defined if this param
                    is used.
                                        
        elDistribVal -
                    Float. Determines how severe the element distribution is.
                    Only works for structured meshes. elOnCurv and 
                    elDistribType must be be defined if this param is used.
                    
                        bump:
                    Smaller value means elements are bunched up at the edges
                    of the curve, larger means bunched in the middle.
                    
                        progression:
                    The edge of each element along this curve (from starting
                    point to end) will be larger than the preceding one by 
                    this factor.
                    elDistribVal = 2 meaning for example that each line element 
                    in the series will be twice as long as the preceding one.
                    elDistribVal < 1 makes each element smaller than the 
                    preceeding one.                      
        '''
        if len(points) != 4:
            raise IndexError("Ellipse: points must be a list of 4 positive integers denoting point indices")
        self._addCurve("Ellipse", points, ID, marker, elOnCurve, elDistribType, elDistribVal)

        
    def _addCurve(self, name, points, ID, marker, elOnCurve, elDistribType, elDistribVal):
        '''Duplicate code goes here!'''
        if ID==None:
            ID = self._getNewCurveID()
        else:
            self._curveIDspecified = True
        
        if elDistribType != None:
            elDistribType = elDistribType.lower().title() #transform into lowercase except the first letter upper case.
            if elDistribType not in ["Bump", "Progression"]:
                raise ValueError("elDistribType must be a string, either \"bump\" or \"progression\". Curve with ID=%i was incorrect." %ID)
            if elDistribVal == None:
                raise ValueError("If elDistribType is defined then elDistribVal must be given a (float) value")
        
        self.curves[ID] = [name, points, marker, elOnCurve, elDistribType, elDistribVal]
    
    def addSurface(self, outerLoop, holes=[], ID=None, marker=0):
        '''
        Adds a plane surface (flat).
        Parameters:
        outerLoop - List of curve IDs that make up the outer boundary of
                    the surface. The curves must lie in the same plane.
        
        holes     - List of lists of curve IDs that make up the inner
                    boundaries of the surface. The curves must lie in the
                    same plane. 
                    
        ID        - Positive integer ID of this surface. If left unspecified
                    the surface will be assigned the smallest unused surface-ID.
                    It is recommended to specify all surface-IDs or none.
                    
        marker    - Integer. Marker applied to this surface. Default 0.
        '''
        #TODO: Possibly check if outerLoop is an actual loop and if the holes are correct.
        self._addSurf("Plane Surface", outerLoop, holes, ID, marker, isStructured=False)
        
        
    def addRuledSurface(self, outerLoop, ID=None, marker=0):
        '''
        Adds a Ruled Surface (bent surface).
        Parameters:
        outerLoop - List of 3 or 4 curve IDs that make up the boundary of
                    the surface.
                    
        ID        - Positive integer ID of this surface. If left unspecified
                    the surface will be assigned the smallest unused surface-ID.
                    It is recommended to specify all surface-IDs or none.
                    
        marker    - Integer. Marker applied to this surface. Default 0.
        '''
        if len(outerLoop) not in [3, 4]:
            raise IndexError("Ruled Surface: outerloop must be a list of 3 or 4 positive integers denoting curve indices")
        self._addSurf("Ruled Surface", outerLoop, [], ID, marker, isStructured=False)
    
    
    def addStructuredSurface(self, outerLoop, ID=None, marker=0):
        '''
        Adds a Structured Surface.
        Parameters:
        outerLoop - List of 4 curve IDs that make up the boundary of
                    the surface. The curves must be structured, i.e. their
                    parameter 'elOnCurv' must be defined.
                    
        ID        - Positive integer ID of this surface. If left unspecified
                    the surface will be assigned the smallest unused surface-ID.
                    It is recommended to specify all surface-IDs or none.
                    
        marker    - Integer. Marker applied to this surface. Default 0.
        '''
        self._checkIfProperStructuredQuadBoundary(outerLoop, ID)
        self._addSurf("Ruled Surface", outerLoop, [], ID, marker, isStructured=True) 
    
    
    def _addSurf(self, name, outerLoop, holes, ID, marker, isStructured):
        '''For duplicate code'''
        #TODO: check if the curves in outerloop actually exist. Maybe print a warning.
        if ID==None:
            ID = self._getNewSurfaceID()
        else:
            self._surfaceIDspecified = True
            
        for hole in holes: #Catch the easy mistake of making holes a list of ints rather than a list of lists of ints.
            try:
                [h for h in hole]
            except TypeError:
                raise TypeError("Hole " + str(hole) + " is not iterable. Parameter holes must be a list of lists of integers")
        
        self.surfaces[ID] = [name, outerLoop, holes, ID, marker, isStructured]
    
    
    def addVolume(self, outerSurfaces, holes=[], ID=None, marker=0):
        '''Adds a Volume
        Parameters:
        outerSurfaces - List of surface IDs that make up the outer boundary of
                        the volume.
        
        holes         - List of lists of surface IDs that make up the inner
                        boundaries of the volume.
                    
        ID            - Positive integer ID of this volume. If left unspecified
                        the volume will be assigned the smallest unused volume-ID.
                        It is recommended to specify all volume-IDs or none.
                    
        marker        - Integer. Marker applied to this volume. Default 0.'''
        self._addVolume(outerSurfaces, holes, ID, marker, isStructured=False)
    
    
    def addStructuredVolume(self, outerSurfaces, ID=None, marker=0):
        '''Adds a Structured Volume
        Parameters:
        outerSurfaces - List of surface IDs that make up the outer boundary of
                        the volume. The surfaces must be Structured Surfaces.
                    
        ID            - Positive integer ID of this volume. If left unspecified
                        the volume will be assigned the smallest unused volume-ID.
                        It is recommended to specify all volume-IDs or none.
                    
        marker        - Integer. Marker applied to this volume. Default 0.'''
        #TODO: Check input. (see if surfaces are structured)
        self._addVolume(outerSurfaces, [], ID, marker, isStructured=True)
    
    
    def _addVolume(self,  outerSurfaces, holes, ID, marker, isStructured):
        '''Duplicate code'''
        #TODO: Check input (outerSurfaces and holes[i] must be closed surfaces)
        if ID==None:
            ID = self._getNewVolumeID()
        else:
            self._volumeIDspecified = True
        self.volumes[ID] = [outerSurfaces, holes, ID, marker, isStructured] 
        
        
    def setPointMarker(self, ID, marker):
        '''Sets the marker of the point with the ID'''
        self.points[ID][2] = marker
    
    
    def setCurveMarker(self, ID, marker):
        '''Sets the marker of the curve with the ID'''
        self.curves[ID][2] = marker
        
        
    def setSurfaceMarker(self, ID, marker):
        '''Sets the marker of the surface with the ID'''
        self.surfaces[ID][4] = marker
        
        
    def setVolumeMarker(self, ID, marker):
        '''Sets the marker of the volume with the ID'''
        self.volumes[ID][3] = marker
    
    
    def _checkIfProperStructuredQuadBoundary(self, outerLoop, ID):
        '''Checks if the four edges of a quad-shaped superelement exist and
        are correct, i.e elOnCurve of opposite curves are equal.'''
        if len(outerLoop) != 4:
            raise IndexError("Structured Surface: outerloop must be a list of 4 positive integers denoting curve indices")
        
        try:
            c0 = self.curves[outerLoop[0]]
            c1 = self.curves[outerLoop[1]]
            c2 = self.curves[outerLoop[2]]
            c3 = self.curves[outerLoop[3]]
        except KeyError:
            raise KeyError("Structured Surface: Attempted construction of StructuredSurface with ID=%s from a curve that does not exist" % ID)
        
        if None in [c0, c1, c2, c3]:
            raise Exception("Attempted to create structured surface from non-structured boundary curves.")
        
        if( c0[-3] != c2[-3] or c1[-3] != c3[-3] ): #Check if the number of elements on opposite curves match.
            raise Exception("Structured Surface: The outerLoop of StructuredSurface %i is not properly " + 
                            "constructed. The reason could be that the number of elements (elOnCurv) on " + 
                            "opposite pairs of curves are different")


    def _getNewPointID(self):
        if not self._pointIDspecified:
            self._nextPointID += 1
            return self._nextPointID - 1
        else:
            return self._smallestFreeKey(self.points)
        
        
    def _getNewCurveID(self):
        if not self._curveIDspecified:
            self._nextcurveID += 1
            return self._nextcurveID - 1
        else:
            return self._smallestFreeKey(self.curves)
        
        
    def _getNewSurfaceID(self):
        if not self._surfaceIDspecified:
            self._nextsurfaceID += 1
            return self._nextsurfaceID - 1
        else:
            return self._smallestFreeKey(self.surfaces)
        
        
    def _getNewVolumeID(self):
        if not self._volumeIDspecified:
            self._nextvolumeID += 1
            return self._nextvolumeID - 1
        else:
            return self._smallestFreeKey(self.volumes)
    
    
    def _smallestFreeKey(self, dictionary):
        '''Finds the smallest unused key in the dict.'''
        sortedkeys = sorted(dictionary)
        for i in range(len(dictionary)):
            if sortedkeys[i] != i:
                return i
