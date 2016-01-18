#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
# Author : Adrien BRUNETON (CEA/DEN)

import fileinput
import re, os
import argparse

from optparse import OptionParser
import os

DEFAULT_EXT = ['.hxx', '.cxx', '.txx', '.py', '.i', '.dox', '.rst', '.h', '.hh', '.hpp', '.c', '.cpp', ]

## The API changes:
REPLACEMENTS = [("RevIntegral",  "IntensiveConservation"), 
                ("ConservativeVolumic", "IntensiveMaximum"),
                ("IntegralGlobConstraint", "ExtensiveConservation"),
                ("Integral", "ExtensiveMaximum"),
                ("MEDCouplingAutoRefCountObjectPtr", "MCAuto"),
                ("deepCpy", "deepCopy") 
                ]   

__myName = os.path.abspath(__file__)

def convert(root_dir, ext, quiet): 
  if not os.path.isdir(root_dir):
    raise ValueError("%s is not a valid directory!" % root_dir)
  
  for root, _, fNames in os.walk(root_dir, followlinks=False):
    for fName in fNames:
      fileName = os.path.join(root, fName)
      # Skip this script!
      if fileName == __myName:
        if not quiet:
          print "!!! Skipping script %s !!!" % __myName
        continue
      ok = False
      for e in ext:
        if fileName[-len(e):] == e:
          ok = True
          break
      if not ok: continue # skip file
      if not quiet:  print "Handling %s ..." % fileName
      for line in fileinput.input(fileName, inplace=1, backup='.bak'):
        for before, after in REPLACEMENTS:
          line = re.sub("(\W|^)(%s)(\W|$)" % before, r"\1%s\3" % after, line.rstrip())
        print(line)  # print in file

parser = OptionParser()
parser.set_usage("Automatic code renaming, corresponding to the API changes between version 7x and 8x of MEDCoupling.\n"
                 "Original files are backup adding '.bak' to the original name.\n"
                 "On Linux file systems, links are not followed!\n\n"
                 "   %prog [options] root_directory")
ext_disp = ",".join(DEFAULT_EXT)
parser.add_option("-e", "--extension", dest="extension", default="", metavar="EXTENSIONS",
                  help="Comma separated list of handled extensions (default is %s)" % ext_disp)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="Don't print status messages to stdout")
(opts, args) = parser.parse_args()

if len(args) != 1:
  parser.print_usage()
  exit(1)

if opts.extension != "":
  ext = [s.strip() for s in opts.extension.split(",")]
else:
  ext = DEFAULT_EXT
  
dirName = os.path.abspath(args[0])
convert(dirName, ext, not opts.verbose)
  
