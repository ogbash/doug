#! /usr/bin/env python

import MySQLdb
import sys
import getopt
import logging
logging.basicConfig()
LOG = logging.getLogger('testcompare')

from doug import testmysql
from doug.config import DOUGConfigParser

def usage():
    sys.stderr.write("Usage: %s\n"
                     "\t [--conf=<filename>]...\n"
                     % (sys.argv[0],))


try:
    opts, extra = getopt.getopt(sys.argv[1:], "", ["conf="])
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(extra)>0:
    usage()
    sys.exit(1)

confFileNames = []
for opt in opts:
    if opt[0]=="--conf":
        confFileNames.append(opt[1])

# config and tar file
conf = DOUGConfigParser()
conf.read(confFileNames)

host = conf.get("testscript", "mysql-host")
user = conf.get("testscript", "mysql-user")
password = conf.get("testscript", "mysql-password")
database = conf.get("testscript", "mysql-database")

class TestResult:
    def __init__(self, **kvargs):
        self.__dict__.update(kvargs)

    def __str__(self):
        return "TestResult(ID=%d)"%self.ID

    def getKeyHeader(self):
        return ('name','inputtype','solver','method','levels','overlap','smoothers','nproc')        

    def getKey(self):
        return tuple(map(lambda attr: getattr(self, attr), self.getKeyHeader()))

def readTestResults():
    cursor = connection.cursor()

    # get all test runs
    cursor.execute("select id from testruns")
    ids = list(cursor)
    ids = map(lambda x: x[0], ids)

    req="SELECT tru.*, count(tr.ID) ntestresults, count(tr.status=3 || NULL) ntestfails, count(tr.status=4 || NULL) ntesterrors" \
         + " FROM testruns tru LEFT JOIN testresults tr" \
         + " ON tru.ID=tr.testrun_ID WHERE tru.ID IN (" + ",".join(map(str,ids)) + ")" \
         + " GROUP BY tru.ID" \
         + " ORDER BY tru.starttime DESC"

    cursor.execute(req)
    testruns = list(cursor)

    # ask for 2 test runs
    for i,r in enumerate(testruns):
        print "%d: "%i, map(str,r)
    print "select two test runs (space sep): 0..%d" % (len(testruns)-1)
    runIds = map(int,map(str.strip, sys.stdin.readline().split(" ")))

    dcursor = connection.cursor(MySQLdb.cursors.DictCursor)

    # get test results
    testResultss = []
    for run_id in runIds:
        req = "SELECT * FROM testresults WHERE testrun_ID=%d" % testruns[run_id][0]
        dcursor.execute(req)
        testResultss.append(map(lambda d: TestResult(**d), list(dcursor)))

    return testResultss

def compareTestResults(testResultss):
    # compare test results
    testDict = {}
    # classify by key
    LOG.info("Scanning 1st test run results")
    for testResult in testResultss[0]:
        key = testResult.getKey()
        if testDict.has_key(key):
            LOG.error("Cannot have the same test in one test run: %s, %s", testResult, testDict[key])
            raise Exception()
        testDict[key] = testResult

    def __compare(tr1, tr2):
        if tr1.status!=tr2.status:
            print "Error for %d,%d: unequal statuses %s,%s for\n%s\n%s" % (
                tr1.ID, tr2.ID,
                testmysql.TestStatus.STRS[tr1.status],
                testmysql.TestStatus.STRS[tr2.status],
                "\t".join(map(str,tr1.getKeyHeader())),
                "\t".join(map(str,tr1.getKey())))
        elif tr1.iterations!=tr2.iterations:
            print "Error for %d,%d: unequal iterations %s,%s for\n%s\n%s" % (
                tr1.ID, tr2.ID,
                tr1.iterations, tr2.iterations,
                "\t".join(map(str,tr1.getKeyHeader())),
                "\t".join(map(str,tr1.getKey())))            

    LOG.info("Comparing with 2nd test run results")
    for testResult in testResultss[1]:
        key = testResult.getKey()
        if testDict.has_key(key):
            __compare(testResult, testDict[key])
        else:
            print "Error: test is missing for %s" % str(map(str,key))

def login():
    global connection
    connection = MySQLdb.connect(host, user, password, database)

login()

if len(sys.argv)>2 and sys.argv[1]=='--getctl':
    cursor = connection.cursor()

    testResultId = int(sys.argv[2])
    req = "SELECT ID,name FROM testresultfiles WHERE testresult_ID=%d" % testResultId
    cursor.execute(req)
    files = list(cursor)
    for i,f in enumerate(files):
        print "%d: %s" % (i,f)
    print "Select file:"
    iFile = int(sys.stdin.readline())
    req = "SELECT content FROM testresultfiles WHERE ID=%d" % files[iFile][0]
    cursor.execute(req)
    content = list(cursor)[0][0]
    print "File content:"
    print content
else:
    testResultss = readTestResults()
    compareTestResults(testResultss)
