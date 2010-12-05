#! /usr/bin/env python

import sys
import getopt
import tarfile
import re
import pickle
from doug import testtar, testmysql
from doug.config import DOUGConfigParser
from datetime import datetime

def usage():
    sys.stderr.write("Usage: %s\n"
                     "\t [--conf=<filename>]...\n"
                     "\t <tar tests file>\n"
                     % (sys.argv[0],))


try:
    opts, extra = getopt.getopt(sys.argv[1:], "", ["conf="])
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(extra)>1 or len(extra)<1:
    usage()
    sys.exit(1)

confFileNames = []
for opt in opts:
    if opt[0]=="--conf":
        confFileNames.append(opt[1])

# config and tar file
conf = DOUGConfigParser()
conf.read(confFileNames)

tarFilename = extra[0]
tarFile = tarfile.TarFile(tarFilename)

host = conf.get("testscript", "mysql-host")
user = conf.get("testscript", "mysql-user")
password = conf.get("testscript", "mysql-password")
database = conf.get("testscript", "mysql-database")

store = testmysql.MysqlTestStore(host, user, password, database)

def createTestRun():
    gconfFile = tarFile.extractfile(testtar.GLOBAL_CONF)
    gconfContent = gconfFile.read()
    config = DOUGConfigParser()
    config.addConfigContents(gconfContent)
    
    startTime = float(config.get("testrun", "starttime"))
    stopTime = float(config.get("testrun", "stoptime"))
    store.insertTestRun(config, datetime.fromtimestamp(startTime))
    store.updateTestRunStop(datetime.fromtimestamp(stopTime))

def createTest(testdir, members):
    del members[testtar.TEST_CONF]
    testConfigFile = tarFile.extractfile("/".join([testdir.name, testtar.TEST_CONF]))
    testConfigContent = testConfigFile.read()
    del members[testtar.RESULT_CONF]
    testResultFile = tarFile.extractfile("/".join([testdir.name, testtar.RESULT_CONF]))
    testResultContent = testResultFile.read()

    config = DOUGConfigParser()
    config.addConfigContents(testConfigContent)
    config.addConfigContents(testResultContent)
    
    # store test times and status
    testname = config.get("test", "name")
    startTime = float(config.get("test", "starttime"))
    stopTime = float(config.get("test", "stoptime"))

    ID = store.insertTest(testname, datetime.fromtimestamp(startTime), config)
    store.updateTestStop(ID, datetime.fromtimestamp(stopTime))
    status = config.get("doug-result", "status").upper()
    status = testmysql.TestStatus.STRS.index(status)
    store.updateTestStatus(ID, status)
    
    # store error
    if members.has_key(testtar.EXCEPTION_FILE):
        del members[testtar.EXCEPTION_FILE]
        exceptionFile = tarFile.extractfile("/".join([testdir.name, testtar.EXCEPTION_FILE]))
        exc = pickle.load(exceptionFile)
        store.updateTestError(ID, str(exc))

    # store profile info
    if config.has_section('doug-profile'):
        store.updateTestProfileInfo(ID, config)

    # left files should be stored to database
    for name in members.keys():
        file = tarFile.extractfile("/".join([testdir.name, name]))
        fileContent = file.read()
        store.insertTestFileContent(ID, name, fileContent)

# main
createTestRun()

# collect tests info
testRE = '([0-9]{4})/(.*)'
members = tarFile.getmembers()
dirFiles = {}
for member in members:
    m = re.match(testRE, member.name)
    if not member.isdir() and m:
        dirName = m.group(1)
        fileName = m.group(2)
        if not dirFiles.has_key(dirName):
            dirFiles[dirName] = {}
        dirFiles[dirName][fileName] = member

# insert tests
testDirRE = '([0-9]{4})$'
dirs = filter(lambda x: x.isdir() and re.match(testDirRE, x.name), members)
for testdir in dirs:
    createTest(testdir, dirFiles[testdir.name])

