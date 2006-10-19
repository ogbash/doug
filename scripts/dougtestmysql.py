
from scripts import ScriptException
import unittest
import logging
import os
import MySQLdb
import datetime
from datetime import datetime

LOG = logging.getLogger("dougtestmysql")

class TestStatus:
    NONE=None
    RUNNING=1
    SUCCESS=2
    FAILURE=3
    ERROR=4

class DougMySQLTestResult(unittest.TestResult):
	"""Saves DOUG test results into MySQL DB."""
	
	def __init__(self, host, user, password, database):
		unittest.TestResult.__init__(self)
                self.host = host
                self.user = user
                self.database = database
                self.connection = MySQLdb.connect(host, user, password, database)
                self.cursor = self.connection.cursor()
                self._createTestRun()

	def startTest(self, test):
            unittest.TestResult.startTest(self, test)
            self.cursor.execute("insert into testresults (name, testrun_ID, starttime, "
                                " method, solver, nproc) values (%s, %s, %s, %s, %s, %s)",
                                (test.testname, self.ID, datetime.today(), test.method, test.solver, test.nproc))
            test.ID = self.cursor.lastrowid
            self.cursor.execute("update testresults set status=%s where id=%s",
                                (TestStatus.RUNNING, test.ID))

	def stopTest(self, test):
            unittest.TestResult.stopTest(self, test)
            self.cursor.execute("update testresults set endtime=%s where id=%s",
                                (datetime.today(), test.ID))

	def addError(self, test, err):
		unittest.TestResult.addError(self, test, err)
                self.cursor.execute("update testresults set status=%s where id=%s",
                                    (TestStatus.ERROR, test.ID))
                self._storeError(test, err)

	def addFailure(self, test, err):
		unittest.TestResult.addFailure(self, test, err)
                self.cursor.execute("update testresults set status=%s where id=%s",
                                    (TestStatus.FAILURE, test.ID))
                self._storeError(test, err)

	def addSuccess(self, test):
		unittest.TestResult.addSuccess(self, test)
                self.cursor.execute("update testresults set status=%s where id=%s",
                                    (TestStatus.SUCCESS, test.ID))

	def stop(self):
		unittest.TestResult.stop(self)

	def close(self):
            try:
                self._finishTestRun()
            finally:
                try:
                    self.cursor.close()
                except Exception, e:
                    LOG.warning("Could not close sursor: %s" % e)

                self.connection.close()

        def _createTestRun(self):
            self.cursor.execute("insert into testruns (servername, starttime)"
                                " values (%s, %s)", ('kheiron', datetime.today()))
            self.ID = self.cursor.lastrowid

        def _finishTestRun(self):
            self.cursor.execute("update testruns set endtime=%s"
                                " where id=%s", (datetime.today(), self.ID))

        def _storeError(self, test, err):
		errtp, errval, errtb = err
                self.cursor.execute("update testresults set errortext=%s where id=%s",
                                    (str(errval), test.ID))
                
		if issubclass(errtp, ScriptException):
			for fname, descr in errval.files:
				self._addFile(fname, test)

        def _addFile(self, name, test):
            self.cursor.execute("insert into testresultfiles (testresult_ID, name)"
                                " values (%s, %s)", (test.ID, os.path.basename(name)))
            fileID = self.cursor.lastrowid

            fd = open(name, "r")
            try:
                fcontent = fd.read()
            except IOError, e:
                LOG.error("Could not open output file %s: %s" % (name, e))
            else:
                self.cursor.execute("update testresultfiles set content=%s"
                                    " where id=%s", (MySQLdb.Binary(fcontent), fileID))
