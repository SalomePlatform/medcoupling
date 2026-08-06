# coding=utf-8

# Copyright 2019 EDF R&D
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, you may download a copy of license
# from https://www.gnu.org/licenses/gpl-3.0.
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

"""This package defines the logger of the MEDConverter plugin."""

import logging


class MEDConverterLogger:
    def __init__(self, level=logging.INFO):
        internal_logger = logging.getLogger("med_convert")
        internal_logger.setLevel(level)

        # Avoid installing a new handler each time this module is reloaded.
        if not internal_logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter(" %(message)s"))
            internal_logger.addHandler(handler)

        self._log = internal_logger

    def setLevel(self, level):
        """Set the level of the logger."""
        self._log.setLevel(level)

    def debug(self, msg, *args, **kwargs):
        self._log.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._log.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._log.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._log.error(msg, *args, **kwargs)

    def exception(self, msg, *args, **kwargs):
        self._log.exception(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._log.critical(msg, *args, **kwargs)


logger = MEDConverterLogger()
