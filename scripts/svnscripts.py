#!/usr/bin/env python

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

import popen2
import os
import time
import logging
from scripts import ScriptException

defaultConfig="""
# This are default (ie reference) values for svnscripts helper package.
# See python ConfigParser package to parse this file.

[svnscripts]
# repository URL like svn://...
repository:

# this is place where modules will be checkout and looked for
workingrootdir: .

# comma separated pairs of module name and workingdir values to update in run(),
# eg: doug/trunk doug_trunk, more/module, another/module/name dir_to_checkout_to
# All non-whitespace from comma to the first whitespace is moduel name,
# after until comma is working directory name. You can skip working dir.
modules:

# time to sleep between svn completion polls (in seconds)
svn-waittime: 3

svn-outfilename: svn.out
svn-errfilename: svn.err
"""
from StringIO import StringIO
defaultConfigFile = StringIO(defaultConfig)

if __name__=='__main__':
    logging.basicConfig()

LOG = logging.getLogger('svnscripts')

def getRevision(dir):
    stream = os.popen('svn info %s | sed "/^Revision:/{s/Revision: \\(.*\\)/\\1/; q;} ;d"' % dir)
    s = stream.read()
    revision = int(s)
    return revision

def run(configFileNames, defaultDict={}):
    LOG.info("Running SVN scripts")

    from ConfigParser import SafeConfigParser
    conf = SafeConfigParser(defaultDict)
    conf.readfp(defaultConfigFile)
    conf.read(configFileNames)

    modulesStr = conf.get('svnscripts', 'modules')
    modules = modulesStr.split(",")

    for moduleStr in modules:
        module = moduleStr.strip().split(" ", 1)
        moduleName = module[0]
        workingDir = moduleName
        if len(module)>1:
            workingDir = module[1]

        work = WorkingDirectory(conf, workingDir, moduleName)
        work.checkoutOrUpdate()
        

class WorkingDirectory:
    """Represents one working SVN directory, whereas one don't have to exist if
    checkout is planned.

    - workingrootdir is used as reference location, unless workingDir is absolute
    - repository is used as location for modules"""

    def __init__(self, conf, workingDir, moduleName = None):
        """If workingDir does not exist then moduleName is obligatory and checkout()
        must be called next."""
        self.conf = conf
        self.workingDir = workingDir
        self.moduleName = moduleName

    def checkout(self):
        LOG.info('Running checkout')
        
        curdir = os.getcwd()
        waittime = self.conf.getint('svnscripts', 'svn-waittime')
        workingRootDir = self.conf.get('svnscripts', 'workingrootdir')
        CWD = self.conf.get('svnscripts', 'cwd')
        outfname = os.path.join(CWD, self.conf.get('svnscripts', 'svn-outfilename'))
        errfname = os.path.join(CWD, self.conf.get('svnscripts', 'svn-errfilename'))
        repoloc = self.conf.get('svnscripts', 'repository')+'/'+self.moduleName
        
        os.chdir(workingRootDir)
        LOG.debug('Changed directory to %s', workingRootDir)
        try:
            LOG.info("Checking out %s to directory %s" % (repoloc, self.workingDir))
            svn = popen2.Popen3('svn checkout %s %s > %s 2> %s'%
                                (repoloc, self.workingDir or "",
                                 outfname, errfname))
            while svn.poll() == -1:
                LOG.debug("Waiting other %d seconds for svn to complete" % waittime)
                time.sleep(waittime)
            
            value = svn.poll()
            if value != 0:
                raise ScriptException("Error occured while running svn (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (value, outfname, errfname))
            
        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)            

    def update(self):
        LOG.info('Running update')
        
        curdir = os.getcwd()
        CWD = self.conf.get('svnscripts', 'cwd')
        waittime = self.conf.getint('svnscripts', 'svn-waittime')
        workingRootDir = os.path.join(CWD, self.conf.get('svnscripts', 'workingrootdir'))
        outfname = os.path.join(CWD,self.conf.get('svnscripts', 'svn-outfilename'))
        errfname = os.path.join(CWD,self.conf.get('svnscripts', 'svn-errfilename'))
        
        os.chdir(workingRootDir)
        LOG.debug('Changed directory to %s', workingRootDir)
        try:
            LOG.info("Updating %s" % (self.workingDir))
            svn = popen2.Popen3('svn update %s > %s 2> %s'%
                                (self.workingDir,
                                 outfname, errfname))
            while svn.poll() == -1:
                LOG.debug("Waiting other %d seconds for svn to complete" % waittime)
                time.sleep(waittime)

            value = svn.poll()
            if value != 0:
                raise ScriptException("Error occured while running svn (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (value, outfname, errfname))
            
        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)

    def checkoutOrUpdate(self):
        "Run update if directory exists, otherwise checkout."
        workingRootDir = self.conf.get('svnscripts', 'workingrootdir')
        d = os.path.join(workingRootDir, self.workingDir)
        if os.path.isdir(d):
            self.update()
        else:
            self.checkout()


if __name__=='__main__':
    logging.getLogger().setLevel(logging.DEBUG)

    from ConfigParser import SafeConfigParser
    conf = SafeConfigParser()
    conf.readfp(defaultConfigFile)
    conf.set('svnscripts', 'repository', 'svn://kheiron.at.mt.ut.ee')
    conf.set('svnscripts', 'workingrootdir', 'test')
    #os.system('rm -r test')
    #os.mkdir('test')
    
    svndoug = WorkingDirectory(conf, 'doug_trunk', 'doug/trunk')

    svndoug.checkoutOrUpdate()
    
    #svndoug.checkout()
    #svndoug.update()

    conf.set('svnscripts', 'modules', 'doug/trunk doug, doug_examples')

    run(conf)
