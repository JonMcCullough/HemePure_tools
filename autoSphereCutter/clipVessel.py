#!/usr/bin/env python

import os.path
import sys
from vtk import vtkStructuredPointsReader
# noinspection PyUnresolvedReferences
# import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
# from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonDataModel import (
    # vtkPlane,
    vtkSphere,
    # vtkPolyData
)
from vtkmodules.vtkFiltersCore import (
    vtkClipPolyData,
)
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkIOGeometry import (
    vtkBYUReader,
    vtkOBJReader,
    vtkSTLReader,
    vtkSTLWriter
)
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkRenderingCore import (
    vtkDataSetMapper
)
import json
import sys
# import csv
# import numpy as np

def get_program_parameters():
    import argparse
    description = 'Clip polydata using a plane.'
    epilogue = '''
    This is an example using vtkClipPolyData to clip input polydata, if provided, or a sphere otherwise.
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue)
    parser.add_argument('filename', nargs='?', default=None, help='Optional input filename e.g cow.g.')
    args = parser.parse_args()
    return args.filename

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 clipVessel.py input.stl input.json setRadius output.stl")
        sys.exit(1)  # Exit with a non-zero value to indicate an error

    filePath = sys.argv[1]
    print(filePath)

    if filePath and os.path.isfile(filePath):
        polyData = ReadPolyData(filePath)
        if not polyData:
            polyData = GetSpherePD()
    else:
        polyData = GetSpherePD()
        
    print("loading input stl finished!")
    f = open(sys.argv[2])
    data = json.load(f)

    # print(data['markups'][0]["controlPoints"])
    # Iterating through the json
    # list
    # set single sphere
    sphere = vtkSphere()

    numberOfSpheres = len(data['markups'][0]["controlPoints"])
    # spheres = list()
    for i in range(0, numberOfSpheres):
      # Create a sphere
      print("Cutting caps #: ", i)
      sphere.SetRadius(float(sys.argv[3]))
      sphere.SetCenter(data['markups'][0]["controlPoints"][i]["position"])

      clipper = vtkClipPolyData()
      clipper.SetInputData(polyData)
      clipper.SetClipFunction(sphere)
      clipper.SetValue(0)
      clipper.Update()
      polyData = clipper.GetOutput()
    
    clipMapper = vtkDataSetMapper()
    clipMapper.SetInputData(polyData)
    

    # use stl writer to write clipped cube into stl file
    # change from unstructured grid to polydata
    # stlFilter = vtkDataSetSurfaceFilter()
    # stlFilter.SetInputConnection(clipper.GetOutputPort(1))
    # stlFilter.Update()

    stlWriter = vtkSTLWriter()
    stlWriter.SetFileName(sys.argv[4])
    stlWriter.SetInputConnection(clipper.GetOutputPort())
    stlWriter.Write()

    # Below is only for the visualization purpose
    # Define colors
    # colors = vtkNamedColors()
    # backgroundColor = colors.GetColor3d('steel_blue')
    # boundaryColor = colors.GetColor3d('Banana')
    # clipColor = colors.GetColor3d('Tomato')
    # clipActor = vtkActor()
    # clipActor.SetMapper(clipMapper)
    # clipActor.GetProperty().SetDiffuseColor(clipColor)
    # clipActor.GetProperty().SetInterpolationToFlat()
    # clipActor.GetProperty().EdgeVisibilityOn()

    # # Now extract feature edges
    # boundaryEdges = vtkFeatureEdges()
    # boundaryEdges.SetInputData(polyData)
    # boundaryEdges.BoundaryEdgesOn()
    # boundaryEdges.FeatureEdgesOff()
    # boundaryEdges.NonManifoldEdgesOff()
    # boundaryEdges.ManifoldEdgesOff()

    # boundaryStrips = vtkStripper()
    # boundaryStrips.SetInputConnection(boundaryEdges.GetOutputPort())
    # boundaryStrips.Update()

    # # Change the polylines into polygons
    # boundaryPoly = vtkPolyData()
    # boundaryPoly.SetPoints(boundaryStrips.GetOutput().GetPoints())
    # boundaryPoly.SetPolys(boundaryStrips.GetOutput().GetLines())

    # boundaryMapper = vtkPolyDataMapper()
    # boundaryMapper.SetInputData(boundaryPoly)

    # boundaryActor = vtkActor()
    # boundaryActor.SetMapper(boundaryMapper)
    # boundaryActor.GetProperty().SetDiffuseColor(boundaryColor)

    # # create renderer render window, and interactor
    # renderer = vtkRenderer()
    # renderWindow = vtkRenderWindow()
    # renderWindow.AddRenderer(renderer)
    # interactor = vtkRenderWindowInteractor()
    # interactor.SetRenderWindow(renderWindow)

    # # set background color and size
    # renderer.SetBackground(backgroundColor)
    # renderWindow.SetSize(640, 480)

    # # add our actor to the renderer
    # renderer.AddActor(clipActor)
    # renderer.AddActor(boundaryActor)

    # # Generate an interesting view
    # renderer.ResetCamera()
    # renderer.GetActiveCamera().Azimuth(30)
    # renderer.GetActiveCamera().Elevation(30)
    # renderer.GetActiveCamera().Dolly(1.2)
    # renderer.ResetCameraClippingRange()

    # renderWindow.Render()
    # renderWindow.SetWindowName('CapClip')
    # renderWindow.Render()

    # interactor.Start()


def ReadPolyData(file_name):
    import os
    path, extension = os.path.splitext(file_name)
    extension = extension.lower()
    if extension == '.ply':
        reader = vtkPLYReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtp':
        reader = vtkXMLpoly_dataReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.obj':
        reader = vtkOBJReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.stl':
        reader = vtkSTLReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtk':
        reader = vtkStructuredPointsReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.g':
        reader = vtkBYUReader()
        reader.SetGeometryFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    else:
        # Return a None if the extension is unknown.
        poly_data = None
    return poly_data


def GetSpherePD():
    '''
    :return: The PolyData representation of a sphere.
    '''
    source = vtkSphereSource()
    source.SetThetaResolution(20)
    source.SetPhiResolution(11)
    source.Update()
    return source.GetOutput()


if __name__ == '__main__':
    main()
