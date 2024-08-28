# ShapeRecogn
A tool leveraging the MEDCoupling library for recognizing canonical shapes in 3D mesh files.

## Table Of Contents

1. [Quickstart](#quickstart)
2. [Installation](#installation)
3. [Description](#description)

## Quickstart

### With output file

```python
>> import ShapeRecogn as sr
>>
>> shape_recogn_builder = sr.ShapeRecognMeshBuilder("resources/ShapeRecognCone.med")
>> shape_recogn = shape_recogn_builder.recognize()
>> #shape_recogn.save("ShapeRecognCone_areas.med")
```

### Without output file

```python
>> import ShapeRecogn as sr
>>
>> shape_recogn_builder = sr.ShapeRecognMeshBuilder("resources/ShapeRecognCone.med")
>> shape_recogn = shape_recogn_builder.recognize()
>> radius_field = shape_recogn.getRadius()
>> radius_field
MEDCouplingFieldDouble C++ instance at 0x55f3418c6700. Name : "Radius (Area)".
Nature of field : NoNature.
P1 spatial discretization.

Mesh info : MEDCouplingUMesh C++ instance at 0x55f34171a2d0. Name : "Mesh_1". Mesh dimension : 2. Space dimension : 3.

Array info : DataArrayDouble C++ instance at 0x55f3418b40e0. Number of tuples : 957. Number of components : 1.
[9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, 9.8948383131651862, ... ]
```

## Installation

### Prerequisites

**LAPACKE Library**: ShapeRecogn requires the [LAPACKE](https://www.netlib.org/lapack/lapacke.html) library for numerical computations. Ensure it is installed on your system.

### CMAKE configuration

Run CMake to configure the build. Make sure to enable ShapeRecogn by setting the `-DMEDCOUPLING_ENABLE_SHAPERECOGN=ON` option.

## Description

ShapeRecogn is a tool that leverages the MEDCoupling library to recognize
canonical shapes from 3D meshes using triangle cells.

The tool is based on the thesis work of Roseline Bénière *[Reconstruction d’un modèle B-Rep à partir d’un maillage 3D](https://theses.fr/2012MON20034)*, and it recognizes five canonical shapes:
 - Plane (0)
 - Sphere (1)
 - Cylinder (2)
 - Cone (3)
 - Torus (4)
 - Unknown (5) (When the algorithm failed the recognition)

The recognized shapes are divided into areas within the mesh.
The tool also generates an output file with detailed fields describing the areas and their properties. This makes it easier to analyze or further process the mesh data.
 - Area Id: To distinguish the areas
 - Primitive Type (Area): One of the canonical shape with the id describe above
 - Normal (Area):
    * Normal of a plane area
 - Minor Radius (Area)
    * Minor radius of a torus area
 - Radius (Area)
    * Radius of a sphere area
    * Radius of a cylinder area
    * Radius of the base of a cone area
    * Major radius of a torus area
 - Angle (Area)
    * Angle between the central axis and the slant of a cone area
 - Center (Area)
    * Center of a sphere area
    * Center of a torus area
 - Axis (Area)
    * Central axis of a cylinder area
    * Central axis of a cone area
 - Apex (Area)
    * Apex of a cone area

Some intermediate computation values concerning the nodes are also available:
 - K1 (Node) and K2 (Node): the *[principal curvatures](https://en.wikipedia.org/wiki/Principal_curvature)*
 - Primitive Type (Node) : One of the canonical shape with the id describe above. The primitive type is deduced usint K1 and K2 values.
 - Normal (Node): Normal of the nodes using neighbor nodes

Each field can be retrieved as a `MEDCouplingDoubleField` using the following methods of the `ShapeRecognMesh` class:
 - `getAreaId()`
 - `getAreaPrimitiveType()`
 - `getAreaNormal()`
 - `getMinorRadius()`
 - `getRadius()`
 - `getAngle()`
 - `getCenter()`
 - `getAxis()`
 - `getApex()`
 - `getNodeK1()`
 - `getNodeK2()`
 - `getNodePrimitiveType()`
 - `getNodeNormal()`
