#! /usr/bin/env python

from scripts import ScriptException
import logging
import os
import sys
import re
import unittest
import popen2

__all__ = ['TestCase', 'getDefaultConfig']

_defaultConfig="""
[dougtest]
# directory where doug_main and doug_aggr are located
dougbindir:

# MPI
mpiboot: lamboot
mpirun: mpirun

doug-outfilename: dougtest.out
doug-errfilename: dougtest.err
"""

LOG = logging.getLogger("dougtest")

def getDefaultConfig():
	return _defaultConfig

class ControlFile:	
	re_assignment = re.compile("(\S+)\s+(.*)")

	def __init__(self, filename):
		self.name = filename
		self.options = {}
		self._parse()

	def _parse(self):
		f = open(self.name, 'r')
		try:
			for line in f:
				match = self.re_assignment.match(line)
				if match:
					self.options[match.group(1)] = match.group(2)
		finally:
			f.close()
			
	def save(self, filename):
		f = open(filename, "w")
		try:
			for key in self.options:
				f.write("%s %s\n" % (key, self.options[key]))
		finally:
			f.close()

class TestCase (unittest.TestCase):
	def __init__(self, testname, datadir, ctrlfname, solutionfname, conf, solver, method, nproc):
		unittest.TestCase.__init__(self)
		self.testname = testname
		self.datadir = datadir
		self.ctrlfname = ctrlfname
		self.solutionfname = solutionfname
		self.conf = conf
		self.solver = solver
		self.method = method
		self.nproc = nproc
	
	def setUp(self):
		LOG.debug("Preparing testing environment")
		tmpdir = os.tempnam(None, 'doug-')
		os.mkdir(tmpdir)
		LOG.debug("Temporary directory %s created" % tmpdir)
		self.tmpdir = tmpdir
		try:
			# copy control file
			self.controlFile = ControlFile(self.ctrlfname)
			self.controlFile.options['plotting'] = '0'
			self.controlFile.options['solver'] = str(self.solver)
			self.controlFile.options['method'] = str(self.method)
			self.testctrlfname = os.path.join(self.tmpdir,
							  os.path.basename(self.ctrlfname))
			self.controlFile.save(self.testctrlfname)

			# make hard link to all other files
			files = os.listdir(self.datadir)
			for fname in files:
				fullfname = os.path.join(self.datadir, fname)

				if fname == os.path.basename(self.testctrlfname):
					continue
				if os.path.isdir(fullfname):
					continue
				
				LOG.debug("Creating hard link for %s" % fullfname)
				os.link(fullfname,
					os.path.join(self.tmpdir, fname))

			# run lamboot (temporary solution for developing)
			mpibootname = self.conf.get("dougtest", "mpiboot")
			mpiboot = popen2.Popen3("%s > %s.out 2> %s.err" % (mpibootname, mpibootname, mpibootname))
			res = mpiboot.wait()
			if res:
				raise ScriptException("Error running %s (%d)" % (mpiboot, res))
		except:
			self._clean()
			raise
		
			
	def tearDown(self):
		LOG.debug("Cleaning after test")
		self._clean()

	def _clean(self):
		if hasattr(self, 'tmpdir') and self.tmpdir:
			os.system('rm -rf %s' % self.tmpdir)
			LOG.debug("Temporary directory %s deleted" % self.tmpdir)

	def runTest(self):
		LOG.debug("Running test")
		LOG.info("solver=%d, method=%d, nproc=%d" % (self.solver, self.method, self.nproc))
		mpirun = self.conf.get("dougtest", "mpirun")
		dougbindir = os.path.abspath(self.conf.get("dougtest", "dougbindir"))
		main = os.path.join(dougbindir, "doug_main")
		errfname = os.path.join(self.tmpdir, self.conf.get("dougtest", "doug-errfilename"))
		outfname = os.path.join(self.tmpdir, self.conf.get("dougtest", "doug-outfilename"))
		
		curdir = os.getcwd()
		
		LOG.debug("Changing directory to %s" % self.tmpdir)
		os.chdir(self.tmpdir)
		try:
			LOG.debug("Running %s -np 1 %s -f %s" % (mpirun, main, self.testctrlfname))
			doug = popen2.Popen3('%s -np %d %s -f %s > %s 2> %s'%
					    (mpirun, self.nproc, main, self.testctrlfname,
					     outfname, errfname))
            
			value = doug.wait()
			LOG.debug("Finished %s with code %d" % (mpirun, value))
			if value != 0:
				se = ScriptException("Error occured while running doug (value=%d), "
						     "inspect output files (%s, %s) for error description." %
						     (value, outfname, errfname))
				se.addFile(outfname, "%s standard output" % mpirun)
				se.addFile(errfname, "%s standard error" % mpirun)
				raise se
			
			#res = os.spawnlp(os.P_WAIT, mpirun, mpirun, "-np", str(self.nproc),
			#		 main, "-f", self.testctrlfname)
		finally:
			LOG.debug("Changing directory to %s" % curdir)
			os.chdir(curdir)

	def run(self, result=None):
		if result is None: result = self.defaultTestResult()
		result.startTest(self)
		testMethod = getattr(self, self.__testMethodName)
		try:
			try:
				self.setUp()
			except KeyboardInterrupt:
				raise
			except:
				result.addError(self, self.__exc_info())
				return

			try:
				testMethod()
				result.addSuccess(self)
			except self.failureException:
				result.addFailure(self, self.__exc_info())
			except KeyboardInterrupt:
				raise
			except:
				result.addError(self, self.__exc_info())

			try:
				self.tearDown()
			except KeyboardInterrupt:
				raise
			except:
				result.addError(self, self.__exc_info())
		finally:
			result.stopTest(self)

	def __str__(self):
		return "%s: solver=%d, method=%d, processors=%d" % \
		       (self.testname, self.solver, self.method, self.nproc)

class CombinedTestResult(unittest.TestResult):
	"""Holder for multiple test result objects.

	All test result objects are called in order they are passed constructor."""
	
	def __init__(self, testResults):
		unittest.TestResult.__init__(self)
		self.testResults = testResults

	def startTest(self, test):
		unittest.TestResult.startTest(self, test)
		for testResult in self.testResults:
			testResult.startTest(test)

	def stopTest(self, test):
		unittest.TestResult.stopTest(self, test)
		for testResult in self.testResults:
			testResult.stopTest(test)

	def addError(self, test, err):
		unittest.TestResult.addError(self, test, err)
		for testResult in self.testResults:
			testResult.addError(test, err)

	def addFailure(self, test, err):
		unittest.TestResult.addFailure(self, test, err)
		for testResult in self.testResults:
			testResult.addFailure(test, err)

	def addSuccess(self, test):
		unittest.TestResult.addSuccess(self, test)
		for testResult in self.testResults:
			testResult.addSuccess(test)

	def stop(self):
		unittest.TestResult.stop(self)
		for testResult in self.testResults:
			testResult.stop()

class TestRunner:
	"""A test runner class that uses console to print results or tars result
	files into specified archive.
	"""
	def __init__(self, testResults = None):
		stream = sys.stderr
		self.stream = unittest._WritelnDecorator(stream)
		self.descriptions = True
		self.verbosity = 1
		if not testResults:
			testResults = [unittest._TextTestResult(self.stream, self.descriptions, self.verbosity)]
		self.testResults = testResults

	def run(self, test):
		"Run the given test case or test suite."
		import time
		result = CombinedTestResult(self.testResults)
		
		startTime = time.time()
		test(result)			
		stopTime = time.time()
		
		timeTaken = stopTime - startTime

		if self.verbosity:
			for testResult in result.testResults:
				if isinstance(testResult, unittest._TextTestResult):
					self._printErrors(testResult, timeTaken)

		return result

	def _printErrors(self, result, timeTaken):
		"Ugly hack to support TextTestResult."
		result.printErrors()
		self.stream.writeln(result.separator2)
		run = result.testsRun
		self.stream.writeln("Ran %d test%s in %.3fs" %
				    (run, run != 1 and "s" or "", timeTaken))
		self.stream.writeln()
		if not result.wasSuccessful():
			self.stream.write("FAILED (")
			failed, errored = map(len, (result.failures, result.errors))
			if failed:
				self.stream.write("failures=%d" % failed)
			if errored:
				if failed: self.stream.write(", ")
				self.stream.write("errors=%d" % errored)
			self.stream.writeln(")")
		else:
			self.stream.writeln("OK")
