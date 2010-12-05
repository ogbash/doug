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
import pickle
from StringIO import StringIO
import time
from config import DOUGConfigParser

LOG = logging.getLogger("dougtesttar")

class DougTarTestResult(unittest.TestResult):
	"""Tars DOUG test results into one file."""
	
	def __init__(self, tarfname, conf):
		unittest.TestResult.__init__(self)
		self.tarFileName = tarfname
		self.conf = DOUGConfigParser()
		self.conf.addConfig(conf)
		self.conf.add_section("testrun")
		self.testCount = 0
		self.tarFile = tarfile.open(self.tarFileName, "w")
		self.testDir = None

	def _addConf(self):
		conf = self.conf
		content = StringIO()
		conf.write(content)
		content.seek(0)
		self._addFileGlobal('global.conf', "Configuration for all tests", fd=content)

	def _addTestConf(self):
		conf = self.testConf
		content = StringIO()
		conf.write(content)
		content.seek(0)
		self._addFile('test.conf', "Configuration for the test", fd=content)

	def startTest(self, test):
		unittest.TestResult.startTest(self, test)
		test.acquire()
		self.testConf = DOUGConfigParser()
		self.testConf.add_section("test")
		self.testConf.set("test", "starttime", str(time.time()))
		self.testDir = self._newTarDirectory()

	def stopTest(self, test):
		unittest.TestResult.stopTest(self, test)
		self.testConf.set("test", "stoptime", str(time.time()))
		self.testConf.addConfig(test.dougExecution.config)
		self._addTestConf()
		test.free()
		self.testDir = None

	def addError(self, test, err):
		unittest.TestResult.addError(self, test, err)
		errtp, errval, errtb = err
		self._addFiles(test.files)
		self._addResult(test, 'error')
		self._addException(test, errval)

	def addFailure(self, test, err):
		unittest.TestResult.addFailure(self, test, err)
		self._addFiles(test.files)
		self._addResult(test, 'failure')
		self._addException(test, errval)

	def addSuccess(self, test):
		unittest.TestResult.addSuccess(self, test)
		self._addFiles(test.files)
		self._addResult(test, 'success')

	def _addResult(self, test, status):
		rfilepath = os.path.join(test.dougExecution.workdir,'result.conf')
		test.resultConfig.set('doug-result', 'status', status)
		f = open(rfilepath, 'w')
		test.resultConfig.write(f)
		f.close()
		self._addFile(rfilepath, "Results of the test execution.")

	def _addException(self, test, err):
		rfilepath = os.path.join(test.dougExecution.workdir,'exception.pickle')
		f = open(rfilepath, 'w')
		pickle.dump(err, f)
		f.close()
		self._addFile(rfilepath, "Exception")

	def startTestRun(self, test):
		# python <2.7 does not have it
		#unittest.TestResult.startTestRun(self)
		self.conf.set("testrun", "starttime", str(time.time()))

	def stopTestRun(self, test):
		# python <2.7 does not have it
		#unittest.TestResult.stopTestRun(self)
		self.conf.set("testrun", "stoptime", str(time.time()))
		self._addConf()

	def close(self):
		self.tarFile.close()

	def _newTarDirectory(self):
		import time
		dirName = "%04d" % self.testCount
		self.testCount += 1
		try:
			testDir = self.tarFile.getmember(dirName)
		except KeyError:
			testDir = tarfile.TarInfo(dirName)
			testDir.type = tarfile.DIRTYPE
			testDir.mode = 0755
			testDir.mtime = time.time()
			self.tarFile.addfile(testDir)
		return testDir

	def _addFiles(self, files):
		"Add files from the test directory, where each file is its name and description."
		d = self.testDir
		for fname, descr in files:
			self._addFile(fname, descr)

	def _addFile(self, fname, descr, fd=None):
		LOG.debug("Adding file %s to tar test" % fname)
		d = self.testDir
		arcname = os.path.join(d.name, os.path.basename(fname))
		self._tarAddFile(arcname, fname, fd)

	def _addFileGlobal(self, fname, descr, fd=None):
		LOG.debug("Adding file %s to tar" % fname)
		arcname = os.path.basename(fname)
		self._tarAddFile(arcname, fname, fd)

	def _tarAddFile(self, arcname, fname=None, fd=None):
		if fd is None:
			self.tarFile.add(fname, arcname)
		else:
			tarinfo = tarfile.TarInfo(name=arcname)
			tarinfo.mtime = int(time.time())
			tarinfo.size = len(fd.getvalue())
			self.tarFile.addfile(tarinfo, fd)

