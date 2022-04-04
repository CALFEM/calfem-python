import pickle
import scipy.io
import numpy as np

'''
Handle reading and writing of geometry and generated mesh from the program
'''


def loadGeometry(name):
    with open(name, 'rb') as file:
        test = pickle.load(file)
    return test


def saveGeometry(g, name="Untitled"):
    if not name.endswith(".cfg"):
        name = name + ".cfg"
    with open(name, 'wb') as file:
        pickle.dump(g, file)


def loadMesh(name):
    with open(name, 'rb') as file:
        mesh = pickle.load(file)
    return mesh


def saveMesh(mesh, name="Untitled"):
    if not name.endswith(".cfm"):
        name = name + ".cfm"
    with open(name, 'wb') as file:
        pickle.dump(mesh, file)


def saveArrays(coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict ,name="Untitled"):
    if not name.endswith(".cfma"):
        name = name + ".cfma"
    with open(name, 'wb') as file:
        pickle.dump(coords, file)
        pickle.dump(edof, file)
        pickle.dump(dofs, file)
        #for key in bdofs.items():
        #    print(key, markerDict[key])
        pickle.dump(bdofs, file)
        pickle.dump(elementmarkers, file)
        pickle.dump(boundaryElements, file)
        pickle.dump(markerDict, file)


def loadArrays(name):
    with open(name, 'rb') as file:
        coords = pickle.load(file)
        edof= pickle.load(file)
        dofs = pickle.load(file)
        bdofs = pickle.load(file)
        elementmarkers = pickle.load(file)
        boundaryElements = pickle.load(file)
        markerDict = pickle.load(file)

    return coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict


def saveMatlabArrays(coords, edof, dofs, bdofs, elementmarkers, boundaryElements, markerDict, name="Untitled"):
    if not name.endswith(".mat"):
        name = name + ".mat"
    saveDict = {}
    saveDict["coords"] = coords.astype('double')
    # Convert to CALFEM Edof definition with element number as first index
    new_column = np.arange(1, np.size(edof, 0) + 1)[:, np.newaxis]
    edof = np.append(new_column, edof, axis=1)

    saveDict["edof"] = edof.astype('double')
    saveDict["dofs"] = dofs.astype('double')
   # bdofs = {str(k): v for k, v in bdofs.items()} # MATLAB struct needs keys as strings
    #print(markerDict)
    newBdof = {}
    for name, index in bdofs.items():
        print(name, index)
        if index == 0:
            newBdof["None"] = 0
        else:
            newBdof[markerDict[index]] = name

    saveDict["bdofs"] = newBdof
    elementmarkers = np.asarray(elementmarkers)
    elementmarkers = elementmarkers + 1  # To avoid problems with one indexing in MATLAB
    saveDict["elementmarkers"] = elementmarkers
    scipy.io.savemat(name, saveDict)

