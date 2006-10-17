
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
				    "s%dm%dnp%d" % (test.solver, test.method, test.nproc),
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
