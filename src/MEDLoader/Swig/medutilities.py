# -*- coding: iso-8859-1 -*-
# --
# Copyright (C) 2009-2016  CEA/DEN, EDF R&D
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
# Author : Erwan ADAM (CEA), Anthony GEAY (CEA)
# --

from MEDLoader import *

def my_remove(f):
    from os import remove
    try:
        remove(f)
    except OSError:
        pass
    return

def convert(file_in, driver_in, driver_out, format=1, file_out=None):
    #
    print(file_in)
    #
    if file_out is None:
        file_out = file_in
        if driver_out == "GIBI":
            file_out += ".sauv"
        elif driver_out == "MED":
            file_out += ".med"
        else:
            msg = "Driver out %s is unknown"%(driver_out)
            raise NotImplementedError(msg)
        pass
    print(file_out)
    #
    if driver_in == "GIBI":
        sr = SauvReader.New(file_in)
        mfd= sr.loadInMEDFileDS()
        pass
    elif driver_in == "MED":
        mfd = MEDFileData(file_in)
        pass
    else:
        raise NotImplementedError("Driver in %s is unknown"%(driver_in))
    #
    my_remove(file_out)
    #
    if driver_out == "GIBI":
        sw=SauvWriter.New()
        sw.setMEDFileDS(mfd,0);#0 ?
        sw.write(file_out)
        #
        mesh = mfd.getMeshes()[0]
        mesh_dim = mesh.getSpaceDimension()
        if mesh_dim >= 3:
            from sys import platform
            if platform in ["win32"]:
                f = open(file_out)
                content = f.read()
                f.close()
                content = content.replace("IFOUR  -1", "IFOUR   2")
                content = content.replace("IFOMOD  -1", "IFOMOD   2")
                f = open(file_out, "w")
                f.write(content)
                f.close()
            else:
                cmd  = "sed"
                cmd += ' -e "s/IFOUR  -1/IFOUR   2/g"'
                cmd += ' -e "s/IFOMOD  -1/IFOMOD   2/g"'
                # cmd += ' -e "s/IECHO   1/IECHO   0/g"'
                cmd += ' %s > .dummy'%(file_out)
                cmd += ' && '
                cmd += ' mv -f .dummy %s'%(file_out)
                from os import system
                system(cmd)
                pass
            pass
        #
        if format == 0:
            from castemlauncher import CastemLauncher
            dgibi_stream  = "\n"
            dgibi_stream += "OPTI REST FORMAT '%s' ;\n"%(file_out)
            dgibi_stream += "REST FORMAT;\n"
            file_out = file_out.replace('__format__', '')
            dgibi_stream += "OPTI SAUV '%s' ;\n"%(file_out)
            dgibi_stream += "SAUV ;\n"
            cl = CastemLauncher(dgibi_stream)
            cl.addTmpFiles(file_out+'__format__', "UTILNOTI", "UTILPROC")
            cl.run()
            pass
        return
    elif driver_out == "MED":
        mfd.write(file_out,2)
        return
    else:
        raise NotImplementedError("Driver in %s is unknown"%(driver_in))

def sauv2med(*argv):
    argv = list(argv)
    for arg in argv:
        convert(arg, "GIBI", "MED")
        pass
    return

def med2sauv(*argv):
    argv = list(argv)
    format = 1
    for arg in argv[:]:
        if arg.find('--format') == 0:
            argv.remove(arg)
            try:
                value = arg.split("=")[1]
            except IndexError:
                usage(1)
                pass
            try:
                value = int(value)
            except ValueError:
                usage(1)
                pass
            format = value
            pass
        pass
    for arg in argv:
        convert(arg, "MED", "GIBI", format)
        pass
    return
