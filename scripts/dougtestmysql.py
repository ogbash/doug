
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
	
	def __init__(self, host, user, password, database, conf):
		unittest.TestResult.__init__(self)
                self.host = host
                self.user = user
                self.database = database
                self.conf = conf
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
            revision = self.conf.get("dougtest","info-svn")
            revision = revision or None
            compiler = self.conf.get("dougtest","info-fc")
            mpi = self.conf.get("dougtest","info-mpi")
            
            self.cursor.execute("insert into testruns (servername, svnrevision, fcompiler, mpi, starttime)"
                                " values (%s, %s, %s, %s, %s)",
                                ('kheiron', revision, compiler, mpi, datetime.today()))
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
