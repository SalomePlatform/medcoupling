<a id="readme-top"></a>

# MEDCoupling

**MEDCoupling** is a robust and versatile C++/Python library for advanced mesh
and field manipulation. It is part of the [SALOME platform](https://salome-platform.org)
and co-developed by **CEA** and **EDF**.

Whether you're performing mesh intersections, transferring fields between
meshes, or analyzing simulation results, MEDCoupling offers the building blocks
you need for scientific computing workflows.

---

## Table of Contents

- [About MEDCoupling](#about-medcoupling)
- [Getting Started](#getting-started)
  - [Installation](#installation)
  - [Usage Highlights](#usage-highlights)
  - [Resources](#resources)
- [Developer Guide](#developer-guide)
  - [Building from Source](#building-from-source)
  - [Contributing](#contributing)
  - [Roadmap](#roadmap)

---

## About MEDCoupling

**MEDCoupling** is a powerful library designed for the manipulation and
analysis of structured and unstructured meshes, as well as associated field
data. It supports:

- Mesh intersection, node merging, extraction, and more
- Field creation and transformation
- Data transfer between non-conforming meshes
- I/O in the `.med` (MEDFile) format

The library handles complex geometries with hybrid cell types (tetrahedra,
hexahedra, higher-order elements, etc.), and supports fields defined on nodes,
cells, or Gauss points.

---

## Getting Started

### Installation

Precompiled binaries for multiple platforms (Windows, Linux) are available on
the [SALOME platform website](https://salome-platform.org).

For more flexibility, you can also build MEDCoupling from source—see the
[Developer Guide](#building-from-source) section below.

### Usage Highlights

MEDCoupling provides a rich set of tools for building, analyzing, and
processing meshes and fields in scientific computing workflows:

#### 🔨 Mesh Creation & Editing
- Create structured and unstructured meshes programmatically.
- Build hybrid meshes with mixed element types (e.g., triangles and quads).
- Modify connectivity or extract sub-meshes.

#### 🔁 Mesh Operations
- Merge or aggregate meshes, with optional node deduplication.
- Intersect meshes to compute shared domains or regions (in 2d).
- Compute cell barycenter on any cell type.

#### 🌡️ Field Definition & Transformation
- Define scalar/vector fields on nodes, cells, or integration points.
- Compute measure fields (cell volumes) on any cell type.
- Use constants, expressions, or interpolated data.
- Transfer fields between non-matching meshes (remap).
- Combine or transform fields using arithmetic or geometric logic.

#### 💾 Input/Output & Interoperability
- Read/write `.med` files (MEDFile format), compatible with SALOME, Code_Saturne, etc.
- Support for time-dependent data and multi-step simulations.
- Export data for post-processing in ParaView or SALOME.
- Used a main `ICoCo` mesh and field transfer format for code coupling.

#### 🐍 Python Integration
- Access full functionality via Python bindings.
- Interact easily with NumPy arrays.
- Script mesh operations and field processing in notebooks or automation pipelines.


### Resources

- 📘 **Tutorials**: [MEDCoupling Tutorials on GitHub](https://github.com/SalomePlatform/MEDCoupling_tutos)
- 📚 **Documentation**: [Official MEDCoupling User & API Docs](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html)
  > *Note: Despite the title, this is the main user and API reference documentation.*

---

## Developer Guide

### Building from Source

The recommended way to build MEDCoupling is using SALOME’s official package
manager, **SAT**, which handles all dependencies (e.g., HDF5, MEDFile) and
supports both native and custom configurations.

---

### Contributing

We welcome contributions!

To ensure consistency, please install and use
[`pre-commit`](https://pre-commit.com):

```bash
pip install pre-commit
cd $PATH_TO_MEDCOUPLING/
pre-commit install
```

This will automatically run code formatters and linters on modified files
before each commit. If hooks fail, there are two possibilities:

- It was a formatting issue — the formatter applied changes. You need to `git
  add` the file(s) and retry the commit.
- It was a linter error — the linter will tell you what to fix. Apply the
  corrections, then `git add` and commit again.

If you're in a hurry or blocked by linter issues, you can:

- Skip a specific linter hook:
  ```bash
  SKIP=linter-name git commit -m "wip: commit bypassing linter"
  ```

- Skip all hooks (**not recommended**, as formatting won't be applied):
  ```bash
  git commit -n -m "wip: commit skipping all pre-commit hooks"
  ```

---

### Roadmap

We're actively improving MEDCoupling. Upcoming features and goals include:

#### 🔧 Usability Improvements
- Unified and searchable documentation (via Sphinx + Breathe)
- Executable notebooks replacing static `.rst` tutorials
- Modernized CMake targets
- Native Python builds
- `spack` and `conan` package support

#### 🛠️ Developer Experience
- Clearer modular architecture (`medcoupling_tools`)
- Modern C++ tooling (`clang-tidy`, `clang-format`, etc.)
- GitHub-based CI/CD workflows

---

> For more information about SALOME and the broader ecosystem, visit [salome-platform.org](https://salome-platform.org).
> We look forward to your feedback and contributions!
