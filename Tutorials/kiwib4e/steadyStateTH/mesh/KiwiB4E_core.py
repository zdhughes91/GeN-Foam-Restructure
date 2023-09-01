#!/usr/bin/env python

# Author: Thomas Guilbaud, 8/10/2022

###
### This file is generated automatically by SALOME v9.9.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'D:/Projects/Nerva/KiwiB4E/coreWedge/GeNFoam/rootCase/mesh')
sys.path.insert(0, r'D:/Projects/Nerva/KiwiB4E')

from dimensions import *

####################################################
##       Begin of NoteBook variables section      ##
####################################################
# notebook.set("centralUnloadedFuelElementRadius", 4.7879)
# notebook.set("fuelElementPitch", 19.1516)
# notebook.set("halfFuelElementPitch", "fuelElementPitch/2")
# notebook.set("fuelElementRadius", "halfFuelElementPitch/1.7320508076")
####################################################
##        End of NoteBook variables section       ##
####################################################


#------------------------------------------------------------------------------*
#                                   GEOMETRY
#------------------------------------------------------------------------------*

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import numpy as np

# Initialize the geometry builder
geompy = geomBuilder.New()

# Create origin for vizualisation
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

def polygon(r: float, n: int):
    """
    Create regular polygon.
        r: radius
        n: number of edges
    return a Salome Face
    """
    vertexList, edgeList = [], []

    # Create vertices
    for theta in np.linspace(0, 2*np.pi, n+1)[:-1]:
        vertex = geompy.MakeVertex(r * np.cos(theta), r * np.sin(theta), 0)
        vertexList.append(vertex)

    # Create edges
    for p1, p2 in zip(vertexList[:-1], vertexList[1:]):
        edge = geompy.MakeEdge(p1, p2)
        edgeList.append(edge)
    edge = geompy.MakeEdge(vertexList[-1], vertexList[0])
    edgeList.append(edge)

    # Create face
    face = geompy.MakeFaceWires(edgeList, 1)
    return(face)

# Create fuel element face
hexagonalFace = polygon(fuelElementPitch/np.sqrt(3), 6)


#------------------------------------------------------------------------------*
#                                     MESH
#------------------------------------------------------------------------------*

# Import Salome Mesh builder
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

# Mesh resolution
ns = 1  # Number of segments per hexagonal side
nz = 25 # Number of segments along the Z-axis

# Initialize the mesh and meshing parameters
hexagonalElement_mesh = smesh.Mesh(hexagonalFace, 'hexagonalElement_mesh')
#Regular_1D = hexagonalElement_mesh.Segment()
#Number_of_Segments_1 = Regular_1D.NumberOfSegments(ns)
#Quadrangle_2D = hexagonalElement_mesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)
#Quadrangle_Parameters_1 = Quadrangle_2D.QuadrangleParameters(
#    smeshBuilder.QUAD_QUADRANGLE_PREF, -1, [], []
#)


NETGEN_1D_2D = hexagonalElement_mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Simple_Parameters = NETGEN_1D_2D.Parameters(smeshBuilder.SIMPLE)
status = hexagonalElement_mesh.AddHypothesis(NETGEN_2D_Simple_Parameters)

# Create mesh face from geometry
base = hexagonalElement_mesh.GroupOnGeom(hexagonalFace, 'base', SMESH.FACE)
# edge = hexagonalElement_mesh.GroupOnGeom(hexagonalFace, 'edge', SMESH.NODE)

# Compute the mesh
isDone = hexagonalElement_mesh.Compute()

# Recover the face
[ base ] = hexagonalElement_mesh.GetGroups()

# Extruction of the meshed face
[ base_extruded, base_top ] = hexagonalElement_mesh.ExtrusionSweepObjects(
    [hexagonalElement_mesh], [hexagonalElement_mesh], [hexagonalElement_mesh],
    [ 0, 0, coreHeight/nz ],    # Long per step
    nz,                         # Nb of segments along the extruction
    1,
    [  ],
    0,
    [ 0, 0, 0 ],                # Base point
    [  ],
    0
)

# Create baffle by generating a global face from the extruded volume
nbAdded, hexagonalElement_meshTemp, tempBaffle = hexagonalElement_mesh.MakeBoundaryElements(
    SMESH.BND_2DFROM3D, 'tempBaffle', 'hexagonalElement_mesh', 1, [ base_extruded ]
)
# And removing the inlet and outlet faces (base, base_top)
baffle = hexagonalElement_meshTemp.GetMesh().CutListOfGroups(
    [ tempBaffle ], [ base, base_top ], 'baffle'
)
# Then remove
hexagonalElement_meshTemp.RemoveGroup( tempBaffle )

## some objects were removed
# aStudyBuilder = salome.myStudy.NewBuilder()
# SO = salome.myStudy.FindObjectIOR(salome.myStudy.ConvertObjectToIOR(tempBaffle))
# if SO: aStudyBuilder.RemoveObjectWithChildren(SO)


base.SetName('outletCentralUnloaded')
base_top.SetName('inletCentralUnloaded')
base_extruded.SetName('centralUnloadedFuelElement')

# Create fuel assembly meshes by translation
fuelElementList = [hexagonalElement_mesh, hexagonalElement_meshTemp]
for theta in np.linspace(0, 2*np.pi, 7)[:-1]:
    # Duplication
    fuelElement_mesh = hexagonalElement_mesh.TranslateObjectMakeMesh(
        hexagonalElement_mesh,  # Base mesh
        [
            fuelElementPitch * np.cos(theta + np.pi/6),
            fuelElementPitch * np.sin(theta + np.pi/6),
            0
        ],                      # Translation vector
        1,                      # Create groups
        'fuelElement'           # Name of the mesh
    )
    fuelElementList.append(fuelElement_mesh)

    # Rename the group volume to be merge into fuelElement
    liste = fuelElement_mesh.GetGroups()
    for e in liste:
        print(e.GetName(), e.GetType())
        if (e.GetType() == SMESH.VOLUME):
            e.SetName("fuelElement")
        elif (e.GetType() == SMESH.FACE and "outlet" in e.GetName()):
            e.SetName("outletFuelElement")
        elif (e.GetType() == SMESH.FACE and "inlet" in e.GetName()):
            e.SetName("inletFuelElement")


# Merge the meshes to form a unique mesh
fuelAssembly = smesh.Concatenate(
    fuelElementList,
    1,
    1,
    1e-05,  # Tolerance
    False   # Create sub groups from fuel elements, if False, keep the previous
            # groups names
)
smesh.SetName(fuelAssembly, 'fuelAssembly')




## Set names of Mesh objects
# smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')
# smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
# smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
# smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
