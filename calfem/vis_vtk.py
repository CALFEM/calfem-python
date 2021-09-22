# -*- coding: utf-8 -*-
"""
CALFEM Visualisation module (matplotlib)

Contains all the functions implementing visualisation routines.
"""

import vtk

def draw_mesh(coords, edof, el_type):

    colors = vtk.vtkNamedColors()
    points = vtk.vtkPoints()

    for i, coord in enumerate(coords):
        points.InsertPoint(i, coord)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(edof.shape[0])

    for dofs in edof:
        if (el_type == 4):
            ugrid.InsertNextCell(vtk.VTK_TETRA, 4, dofs-1)
        elif (el_type == 5):
            ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, dofs-1)
        else:
            print("Unsupported element type.")

    ugrid.SetPoints(points)

    ugridMapper = vtk.vtkDataSetMapper()
    ugridMapper.SetInputData(ugrid)

    ugridActor = vtk.vtkActor()
    ugridActor.SetMapper(ugridMapper)
    ugridActor.GetProperty().SetColor(colors.GetColor3d('Peacock'))
    ugridActor.GetProperty().EdgeVisibilityOn()

    renderer = vtk.vtkRenderer()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    renderer.AddActor(ugridActor)
    renderer.SetBackground(colors.GetColor3d('Beige'))

    renderer.ResetCamera()
    renderer.GetActiveCamera().Elevation(60.0)
    renderer.GetActiveCamera().Azimuth(30.0)
    renderer.GetActiveCamera().Dolly(1.0)

    renWin.SetSize(640, 480)
    renWin.SetWindowName('UGrid')

    # Interact with the data.
    renWin.Render()

    iren.Start()


