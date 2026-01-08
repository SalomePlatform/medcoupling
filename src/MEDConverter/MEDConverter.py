# Copyright (C) 2023-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

"""
This package defines the *engine* of the *MEDConverter* plugin.
"""

from .abaqus import MEDConverterAbaqus
from .ansys import MEDConverterAnsys
from .aster import MEDConverterAster
from .systus import MEDConverterSystus
from .tetgen import MEDConverterTetgen
from .zset import MEDConverterZset


# ABAQUS mesh
def LoadINPFileInMEDFileUMeshInstance(inputMailFilePath: str, verbose=False):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterAbaqus.convert_abaqus_to_med(inputMailFilePath, verbose)


def ConvertFromINPToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadINPFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath


# ANSYS mesh
def LoadCDBFileInMEDFileUMeshInstance(inputMailFilePath: str, verbose=False):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterAnsys.convert_ansys_to_med(inputMailFilePath, verbose)


def ConvertFromCDBToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadCDBFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath


# ASTER mesh
def LoadMailFileInMEDFileUMeshInstance(inputMailFilePath: str, verbose=False):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterAster.convert_aster_to_med(inputMailFilePath, verbose)


def ConvertFromMailToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadMailFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath


# SYSTUS mesh
def LoadASCFileInMEDFileUMeshInstance(
    inputMailFilePath: str, skip_types=[], verbose=False
):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterSystus.convert_systus_to_med(
        inputMailFilePath, skip_types, verbose
    )


def ConvertFromASCToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadASCFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath


# TETGEN mesh
def LoadTetgenFileInMEDFileUMeshInstance(inputMailFilePath: str, verbose=False):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterTetgen.convert_tetgen_to_med(inputMailFilePath, verbose)


def ConvertFromTetgenToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadTetgenFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath


# ZSET mesh
def LoadGeofFileInMEDFileUMeshInstance(inputMailFilePath: str, verbose=False):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Build MEDCouplingUMesh from data structures
    return MEDConverterZset.convert_zset_to_med(inputMailFilePath, verbose)


def ConvertFromGeofToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadGeofFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath
