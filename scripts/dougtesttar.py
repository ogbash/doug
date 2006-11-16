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

from scripts import ScriptException
import tarfile
import unittest
import logging
import os

LOG = logging.getLogger("dougtesttar")

class DougTarTestResult(unittest.TestResult):
	"""Tars DOUG test results into one file."""
	
	def __init__(self, tarfname):
		unittest.TestResult.__init__(self)
		self.tarFileName = tarfname
		self.tarFile = tarfile.open(self.tarFileName, "w")

	def startTest(self, test):
		unittest.TestResult.startTest(self, test)
		d = self._getTestDirectory(test)

	def stopTest(self, test):
		unittest.TestResult.stopTest(self, test)

	def addError(self, test, err):
		unittest.TestResult.addError(self, test, err)
		errtp, errval, errtb = err
		d = self._getTestDirectory(test)
		if issubclass(errtp, ScriptException):
			for fname, descr in errval.files:
				self._addFile(fname, descr, test)

	def addFailure(self, test, err):
		unittest.TestResult.addFailure(self, test, err)

	def addSuccess(self, test):
		unittest.TestResult.addSuccess(self, test)

	def stop(self):
		unittest.TestResult.stop(self)

	def close(self):
		self.tarFile.close()

	def _getTestDirectory(self, test):
		import time
		dirName = "/".join([test.testname,
				    test.executable,
				    "s%dm%dl%dnp%d" % (test.solver, test.method, test.levels, test.nproc),
				    ""])
		try:
			testDir = self.tarFile.getmember(dirName)
		except KeyError:
			testDir = tarfile.TarInfo(dirName)
			testDir.type = tarfile.DIRTYPE
			testDir.mode = 0755
			testDir.mtime = time.time()
			self.tarFile.addfile(testDir)
		return testDir

	def _addFile(self, fname, descr, test):
		LOG.debug("Adding file %s to tar" % fname)
		d = self._getTestDirectory(test)
		arcname = "/".join([d.name, os.path.basename(fname)])
		self.tarFile.add(fname, arcname)
