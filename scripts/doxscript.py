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
# or contact the author (University of Tartu, Faculty of Computer Science, Chair
# of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
# mailto:info(at)dougdevel.org)

from scripts import ScriptException
import svnscripts
import autotools
import logging
import getopt
import sys
from ConfigParser import SafeConfigParser
from StringIO import StringIO

defaultConfig = """
[doxscript]
doxygen: /usr/bin/doxygen
docdir: doc
docdirhtml: %(docdir)s/html
doxygen-outfilename: doxygen.out
doxygen-errfilename: doxygen.err

[apache]
docdir: /var/www/docs
"""

logconfFileNames = []
confFileNames = []

def usage():
    sys.stderr.write("Usage:\n"
                     "\t%s --conf=<filename> --logconf=<filename>\n"
                     "\t%s -g (to generate default config to stdout)\n"
                     % (sys.argv[0], sys.argv[0]))

try:
    opts, extra = getopt.getopt(sys.argv[1:], "g", ["logconf=", "conf="])
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(extra)>0:
    usage()
    sys.exit(1)

for opt in opts:
    if opt[0]=="--logconf":
        logconfFileNames.append(opt[1])
    if opt[0]=="--conf":
        confFileNames.append(opt[1])
    if opt[0]=="-g":
        sys.stdout.write(defaultConfig)
        sys.stdout.write(svnscripts.defaultConfig)
        sys.stdout.write(autotools.defaultConfig)
        sys.exit(0)

import logging.config
for fname in logconfFileNames:
    logging.config.fileConfig(fname)
else:
    logging.basicConfig()
    
LOG = logging.getLogger('doxscript')

try:
    import time
    LOG.info("Running doxscript at %s" % time.asctime())

    # svn
    svndoug = svnscripts.run(confFileNames)

    # autotools
    autotools.run(confFileNames)

    # doxygen 
    LOG.info("Preparing to run doxygen")
    import popen2
    import os

    conf = SafeConfigParser()
    conf.readfp(StringIO(defaultConfig))
    conf.read(confFileNames)

    cmd = conf.get('doxscript', 'doxygen')
    docdir = conf.get('doxscript', 'docdir')
    docdirhtml = conf.get('doxscript', 'docdirhtml')
    outfname = conf.get('doxscript', 'doxygen-outfilename')
    errfname = conf.get('doxscript', 'doxygen-errfilename')

    curdir = os.getcwd()
    os.chdir(docdir)
    LOG.debug('Changed directory to %s', docdir)

    try:
        LOG.info("Running doxygen")
        dx = popen2.Popen3("%s > %s 2> %s" % (cmd, outfname, errfname))
        dx.wait()

        dx_value = dx.poll()
        if dx_value != 0:
            raise ScriptException("Error occured while running doxygen (value=%d), "
                                  "inspect output files (%s, %s) for error description." %
                                  (dx_value, outfname, errfname))

    finally:
        os.chdir(curdir)
        LOG.debug('Changed directory to %s', curdir)

    # copying to apache dir
    apache_docdir = conf.get('apache', 'docdir')

    LOG.info("Copying documentation to apache")
    dx = popen2.Popen3("cp -R %s/* %s >> %s 2>> %s" % (docdirhtml, apache_docdir, outfname, errfname))
    dx.wait()

    dx_value = dx.poll()
    if dx_value != 0:
        raise ScriptException("Error occured while copying documentation (value=%d), "
                  "inspect output files (%s, %s) for error description." %
                  (dx_value, outfname, errfname))

    LOG.info("Ended doxscript at %s" % time.asctime())

except ScriptException, e:
    LOG.critical("Error while running script: %s" % e)
    sys.exit(1)
except Exception, e:
    LOG.critical("Unknown exception occured: %s" % e)
    raise
