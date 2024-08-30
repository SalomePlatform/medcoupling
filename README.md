<a id="readme-top"></a>

# MEDCoupling

MEDCoupling is a powerful library to manipulate meshes and fields. 

<!-- TABLE OF CONTENTS -->
<details>
    <summary>Table of Contents</summary>
    <ol>
        <li>
            <a href="#about-medcoupling">About MEDCoupling</a>
        </li>
        <li>
            <a href="#getting-stated">Getting started</a>
        </li>
        <ul>
            <li><a href="#installation">Installation</a></li>
            <li><a href="#usage">Usage</a></li>
        </ul>
        <li>
            <a href="#roadmap">Roadmap</a>
        </li>
        <li>
            <a href="#building-from-source">Building from source</a>
        </li>
    </ol>
</details>

<!-- CONTENT -->

## About MEDCoupling

MEDCoupling is part of the Salome project. It is a powerfull mesh and field
library able to compute mesh intersections, volumes, ...
The library is co-developed between the CEA and EDF.

## Getting started

### Installation

You can find MEDCoupling binaries for multiple platforms on the
[Salome website](https://salome-platform.org). It is available under Windows
and a few Linux distributions.

### Usage

MEDCoupling can be used to:
- create structured/unstructured meshes by hand (no geometrical auto-tetra meshing)
- intersecting meshes, merging nodes, extracting part of a mesh, ...
- creating custom fields on meshes (by constant, expression, measure, ...)
- load and write mesh in the `medfile` file format (`.med` extension)
- coupling codes by transfering fields from one mesh to another
- post-processing fields and mesh informations

The MEDCoupling format is quite versatile and allows to manage mesh with different kind of cells (tetra, hexa, higher order, ...) and fields on cells, on nodes, or on gauss points.

## Roadmap

The MEDCoupling library is under a major refactoring for the v10. It will allow to:
- Make usage easier
    - Better documentation
        - [ ] Build a new unified documentation with Sphinx and its amazing elastic
          search integration
        - [ ] Take advantage of `breathe` to integrate `MEDCoupling` API to the doc
        - [ ] Replace the `.rst` files of the tutorial doc with executable notebooks
    - Better compilation process
        - [ ] Usage of modern cmake and clearer targets names
        - [ ] Allow python native compilation
        - [ ] Adding `spack` and `conan` recipes
    - Empowering medcoupling core mesh abililties
        - [ ] Separate the usefull `MEDFileXXX` objects from the `medfile` dependency
        - [ ] Make `medfile` one of the mesh backends
        - [ ] Add a new file format backend (namely CGNS)
- Make contributions easier
    - Modular architecture
        - [ ] Clarify what is part of medcoupling core data structure and what is not
          with the introduction of `medcoupling_tools`
        - [ ] Separate the remapper to make it into a tool
    - Modern and standard C++
        - [ ] Modernize the whole repository by adding standard tooling (clang-tidy,
          clang-format, pre-commit, ...)
    - Public CI/CD
        - [ ] Taking advantage of GitHub by adding github workflows

## Building from source

For now the recommended way to compile is to use Salome homemade package
manager, `sat`. It will install medcoupling dependencies such as HDF5 and
medfile depending on the choosen configuration (native or not).
