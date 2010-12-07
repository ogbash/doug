
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
    STRS = ['NONE','RUNNING','SUCCESS','FAILURE','ERROR']

class MysqlTestStore:
    "Methods to store test to Mysql database."

    def __init__(self, host, user, password, database):
        self.host = host
        self.user = user
        self.password = password
        self.database = database
        self.connection = MySQLdb.connect(host, user, password, database)
        self.cursor = self.connection.cursor()

    def insertTestRun(self, conf, time):
        servername = conf.get("dougtest", "info-server")
        revision = conf.get("dougtest","info-svn")
        revision = revision or None
        gitversion = conf.get("dougtest","info-git")
        gitversion = gitversion or None
        compiler = conf.get("dougtest","info-fc")
        mpi = conf.get("dougtest","info-mpi")
            
        self.cursor.execute("insert into testruns (servername, svnrevision, gitversion, fcompiler, mpi, starttime)"
                            " values (%s, %s, %s, %s, %s, %s)",
                            (servername, revision, gitversion, compiler, mpi, time))
        self.ID = self.cursor.lastrowid

    def updateTestRunStop(self, time):
        self.cursor.execute("update testruns set endtime=%s"
                            " where id=%s", (time, self.ID))

    def updateTestRunError(self, errdesc):
        self.cursor.execute("update testruns set errortext=%s where id=%s",
                            (errdesc, self.ID))


    def insertTest(self, testname, time, conf):
        "time - test start time"
        self.cursor.execute("insert into testresults (name, testrun_ID, starttime, "
                            " executable, inputtype, method, solver, nproc, levels, overlap, smoothers)"
                            " values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                            (testname, self.ID, time,
                             conf.get("doug", "executable"),
                             conf.get("doug-controls", "input_type"),
                             conf.get("doug-controls", "method"),
                             conf.get("doug-controls", "solver"),
                             conf.get("doug", "nproc"),
                             conf.get("doug-controls", "levels"),
                             conf.get("doug-controls", "overlap"),
                             conf.get("doug-controls", "smoothers")))
        ID = self.cursor.lastrowid
        self.cursor.execute("update testresults set status=%s where id=%s",
                            (TestStatus.RUNNING, ID))
        return ID

    def updateTestStop(self, ID, time):
        self.cursor.execute("update testresults set endtime=%s where id=%s",
                            (time, ID))

    def updateTestStatus(self, ID, status):
        self.cursor.execute("update testresults set status=%s where id=%s",
                            (status, ID))

    def updateTestError(self, ID, errdesc):
        self.cursor.execute("update testresults set errortext=%s where id=%s",
                            (errdesc, ID))

    def updateTestProfileInfo(self, ID, conf):
        c = conf
        sqlr = ("update testresults set iterations=%d, iterations_time=%f, preconditioner_time=%f," + 
                " fine_aggrs=%d, coarse_aggrs=%d where id=%s") % (
            c.getint('doug-profile', 'pcg-iterations'),
            c.getfloat('doug-profile', 'iterations-time'),
            c.getfloat('doug-profile', 'preconditioner-time'),
            c.has_option('doug-profile', 'fine-aggregates') and c.getint('doug-profile', 'fine-aggregates') or -1,
            c.has_option('doug-profile', 'coarse-aggregates') and c.getint('doug-profile', 'coarse-aggregates') or -1,
            ID)
        self.cursor.execute(sqlr)

    def insertTestFileContent(self, ID, filename, content):
        self.cursor.execute("insert into testresultfiles (testresult_ID, name)"
                            " values (%s, %s)", (ID, filename))
        fileID = self.cursor.lastrowid
        self.cursor.execute("update testresultfiles set content=%s"
                            " where id=%s", (MySQLdb.Binary(content), fileID))

                        
class DougMySQLTestResult(unittest.TestResult):
	"""Saves DOUG test results into MySQL DB."""
	
	def __init__(self, host, user, password, database, conf):
		unittest.TestResult.__init__(self)
                self.store = MysqlTestStore(host,user,password,database)
                self.conf = conf
                self.store.insertTestRun(self.conf, datetime.today())

	def startTest(self, test):
            test.acquire()
            unittest.TestResult.startTest(self, test)
            conf = test.dougExecution.config
            test.ID = self.store.insertTest(test.testname, datetime.today(), conf)

	def stopTest(self, test):
            test.free()            
            unittest.TestResult.stopTest(self, test)
            self.store.updateTestStop(test.ID, datetime.today())

	def addError(self, test, err):
            if test!=None:
                # test error
                unittest.TestResult.addError(self, test, err)
                self.store.updateTestStatus(test.ID, TestStatus.ERROR)
                self._storeError(test, err)
            else:
                # generic error. hack :(, must be removed?
		errtp, errval, errtb = err
                self.store.updateTestRunError(str(errval))

	def addFailure(self, test, err):
		unittest.TestResult.addFailure(self, test, err)
                self.store.updateTestStatus(test.ID, TestStatus.FAILURE)
                self._storeError(test, err)

	def addSuccess(self, test):
		unittest.TestResult.addSuccess(self, test)
                self.store.updateTestStatus(test.ID, TestStatus.SUCCESS)

		if hasattr(test, "files"):
			for fname, descr in test.files:
				self._addFile(fname, test)

                # update profiling values
                if test.resultConfig.has_section('doug-profile'):
                    conf = test.resultConfig
                    self.store.updateTestProfileInfo(test.ID, conf)

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

        def _finishTestRun(self):
            self.store.updateTestRunStop(datetime.today())

        def _storeError(self, test, err):
            errtp, errval, errtb = err
            self.store.updateTestError(test.ID, str(errval))
                
            if hasattr(test, "files"):
                for fname, descr in test.files:
                    self._addFile(fname, test)

        def _addFile(self, name, test):
            fd = open(name, "r")
            try:
                fcontent = fd.read()
            except IOError, e:
                LOG.error("Could not open output file %s: %s" % (name, e))
            else:
                self.store.insertTestFileContent(test.ID, os.path.basename(name), fcontent)
            fd.close()
