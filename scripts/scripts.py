#!/usr/bin/env python

# DOUG - Domain decomposition On Unstructured Grids
# Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
# Department of Mathematics, University of Bath
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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or contact the authors (University of Tartu, Faculty of Computer Science, Chair
# of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
# mailto:info(at)dougdevel.org)

"""Common classes used by all DOUG scripts"""

class ScriptException (Exception):
    """Script running error that is understood by script, ie failure running some programs.

    Some day later, additional information can be added to this class that makes possible logs,
    error files, input files or other information that can help understand error source."""

    def __init__(self, descr, causeException=None):
        Exception.__init__(self, descr)
        self.cause = causeException
        self.files = []

    def addFile(self, fname, description=None):
        self.files.append((fname, description))

    def __str__(self):
        str = Exception.__str__(self)
        if self.cause:
            return "\n".join([str, self.cause.__str__()])
        else:
            return str
                   
