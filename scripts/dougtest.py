#! /usr/bin/env python

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
import logging
import os
import sys
import re
import unittest
import popen2

_defaultConfig="""
[dougtest]
# all info-* attributes are not necessary, but may be used by test results
# name of the server tests are run on
info-server:
# DOUG svn version
info-svn:
# fortran compiler
info-fc:
# MPI version
info-mpi:

"""

__all__ = ['TestCase', 'getDefaultConfig']


LOG = logging.getLogger("dougtest")

def getDefaultConfig():
	return _defaultConfig


class TestFailure (ScriptException):
	pass

class TestCase (unittest.TestCase):
	failureException=TestFailure

	def __init__(self, dougExecution):
		unittest.TestCase.__init__(self, '_test')
		self.dougExecution = dougExecution

	def setUp(self):
		self.dougExecution.setUp()

	def tearDown(self):
		self.dougExecution.tearDown()


	def _test(self):
		self.dougExecution.run()

	def __run(self, result=None):
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

        def _assertSolution(self):
            from array import array
            solCorrect = array('d')
            sol = array('d')
	    # we need 4 byte integer, ugly hack for 64bit platforms
            intarr = array('l')
	    if intarr.itemsize != 4:
		    intarr = array('i')

            solfilename = os.path.abspath(self.solutionfname)
            f=open(solfilename, "rb")
            try:
                try:
                    intarr.fromfile(f, 1)
                    smarker = intarr[-1]
                    sol.fromfile(f, smarker/sol.itemsize)
                    intarr.fromfile(f, 1)
                    emarker = intarr[-1]
                except Exception, e:
                    se = ScriptException("Error reading solution, investigate '%s' file."
                                         % os.path.basename(solfilename), e)
		    se.addFile(solfilename, "solution file")
                    raise se
            finally:
                f.close()
            
            f=open(self.csolutionfname, "rb")
            try:
                intarr.fromfile(f, 1)
                smarker = intarr[-1]
                solCorrect.fromfile(f, smarker/solCorrect.itemsize)
                intarr.fromfile(f, 1)
                emarker = intarr[-1]
            finally:
                f.close()

	    self.assertAlmostEqual(sol, solCorrect, 1E-6)

	def assertAlmostEqual(self, v1, v2, acceptedTolerance):
		if len(v1) != len(v2):
			raise self.failureException("Sizes of solutions are different: %d, %d"
						    % (len(v1), len(v2)) )
		tolerance = 0.
		for i in xrange(0,len(v1)):
			tolerance += pow(v1[i]-v2[i], 2)
		tolerance = pow(tolerance, 0.5)

		if tolerance > acceptedTolerance:
			raise self.failureException("Solution is wrong with difference: %f"
						    % (tolerance) )

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
