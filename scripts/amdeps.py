#!/usr/bin/env python

import sys
import re
import os.path
        
class DevNull:
    def __getattr__(self, name):
        return self.null

    def null(self, *args):
        pass

warn = lambda x: sys.stderr.write('WARN: %s\n' % x)
info = lambda x: sys.stderr.write('INFO: %s\n' % x)
debug = lambda x: sys.stderr.write('DEBUG: %s\n' % x)

debug = DevNull()

# Targets
class Target:
    """Class representing one target (file) in Makefile. It's name is
    stored in target and source (if any) in source. On creation it
    must read the source file and calculate it's dependency and
    provision logical names."""
    
    def __init__(self, targetName, sourceName, targetDir = '', sourceDir = ''):
        """targetDir - directory of Makefile, needed becuase current
        working directory may be different and target may not find
        it's source."""
        self.targetName = targetName
        self.sourceName = sourceName
        self.targetDir = targetDir
        self.sourceDir = sourceDir
    
    def getDeps(self):
        """Method return logical names that this target depends on
        (these may be file names, eg my.h, or fortran module names)."""
        raise NotImplementedError("Must be implemented in derived class")

    def getProvs(self):        
        """Method return logical names that this target provides
        (these may be file names, eg my.h, or fortran module names)."""
        raise NotImplementedError("Must be implemented in derived class")

def createTarget(targetName, sourceName, targetDir='', sourceDir=''):
    "Tries to match sourceName extension and target class."

    from os.path import splitext
    
    if splitext(sourceName)[1] in [".f90",".F90",".f95",".F95"]:
        return F9xTarget(targetName, sourceName, targetDir, sourceDir)

    warn("Unresolved target class for %s" % sourceName)

    return None

class F9xTarget(Target):
    re_hasmod = re.compile("^\s*(!?)\s*module\s+(\S*)")
    re_usesmod = re.compile("^\s*(!?)\s*use\s+([^,\s]*),*")

    def __init__(self, targetName, sourceName, targetDir='', sourceDir=''):
        # TODO oleg:all add targetName, sourceName resolving code if None

        Target.__init__(self, targetName, sourceName, targetDir, sourceDir)
        self._scanProvsAndDeps()

    def getProvs(self):
        return self.modprovs

    def getDeps(self):
        return self.moddeps

    def _scanProvsAndDeps(self):
        "Scan Fortran 90/95 source file for used and contained modules."
        f = open(os.path.join(self.sourceDir, self.sourceName), 'r')
        modhas = []
        moduses = []
        for line in f:
            line = line.lower()
            
            hasmod = F9xTarget.re_hasmod.match(line)
            
            if hasmod:
                if hasmod.group(1) != '!':
                    modhas.append(hasmod.group(2))
            
            usesmod = F9xTarget.re_usesmod.match(line)
            if usesmod:
                if (usesmod.group(1) != '!') and (usesmod.group(2) not in moduses):
                    moduses.append(usesmod.group(2))

        self.modprovs = modhas
        self.moddeps = moduses

    def __str__(self):
        return "Fortran target: %s, %s\n\tmodules provides = %s\n\tmodules depends = %s\n" % (self.targetName, self.sourceName, self.modprovs, self.moddeps)

# Dependency collectors
class Collector:
    """Derived classes must take some kind of info and generate
    dependencies."""
    
    def getDependencies(self):

        
        "Dictionary of file dependencies (or any other?)."
        raise NotImplementedError("Method purpose and exact meaning is revised")

def am_canonicalize(name):    
    "Canonicalize name like Automake does it."
    import string
    # all chars except a-zA-Z_@ are turned into underscore
    l = list(name)
    l_letters = list(string.letters)
    l_digits = list(string.digits)
    
    for i,c in enumerate(l):
        if not c in l_letters and \
           not c in l_digits and \
           not c in ['_','@']:
            l[i] = '_'

    return "".join(l)


re_useVar = re.compile(r'\$[({]?[_@A-Za-z]+[)}]?')
def am_subst(str, vars={}):
    "Replace any variables with values (which may be stupid $(EXEEXT))"
    l_resStr = []

    i = 0
    match = re_useVar.search(str, i)
    while match:
        l_resStr.append(str[i:match.start()])
        if vars.has_key(match.group(0)):
            l_resStr.append(vars[match.group(0)])
        i = match.end()
        match = re_useVar.search(str, i)

    l_resStr.append(str[i:])
    return "".join(l_resStr)
    

class AutomakeCollector(Collector):
    """Takes Automake generated Makefile and AM BIN target name. Then
    scans Makefile for sources and objects and tries to collect all
    dependencies."""

    re_prgs = re.compile(r"^([_@a-zA-Z]+)_PROGRAMS =")
    re_ltlibs = re.compile(r"^([_@a-zA-Z]+)_LTLIBRARIES =")
    re_obj = re.compile(r"^([_@a-zA-Z]+)_OBJECTS =")
    re_src = re.compile(r"^([_@a-zA-Z]+)_SOURCES =")
    re_var = re.compile(r"^([_@a-zA-Z]+) =(?P<val>.*)\n?$")
    re_varCont = re.compile(r"^\s+(?P<val>.*)\n?$")
    re_target = re.compile(r"(\S+)\s*:\s*(\S+)\s*$")

    def __init__(self, fileName="Makefile", binTargets=None, sourceDir=''):
        (self.fileDir, self.fileName) = os.path.split(fileName)
        self.binTargets = binTargets
        self.sourceDir = sourceDir
        self.targets = {}
        
        self._scanMakefile()

    def _scanMakefile(self):
        info("Scanning %s" % os.path.join(self.fileDir, self.fileName))
        
        fd = open(os.path.join(self.fileDir, self.fileName), 'r')

        try:
            # if binTargets is null then scan Makefile and find all
            if self.binTargets == None:
                self.binTargets = []
                fd.seek(0)
                iter_ = iter(fd)
                for line in iter_:
                    # match <target>_PROGRAMS
                    match = self.re_prgs.match(line)
                    if match:
                        var=self._getVariable(line, fd)
                        prgs = var[1].split()
                        self.binTargets.extend(prgs)
                    # match <target>_LTLIBRARIES
                    match = self.re_ltlibs.match(line)
                    if match:
                        var=self._getVariable(line, fd)
                        ltlibs = var[1].split()
                        self.binTargets.extend(ltlibs)
            
            # in first scan get sources
            self.srcs = []
            self.objs = []

            
            binTargets_vr = map(am_subst, self.binTargets) # variables resolved
            info("Bin targets are %s" % str(binTargets_vr))
            binTargets_cn = map(am_canonicalize, binTargets_vr) # canonicalized
            
            fd.seek(0)
            iter_ = iter(fd)
            for line in iter_:
                # match <target>_SOURCES
                match = self.re_src.match(line)
                if match and match.group(1) in binTargets_cn:
                    var=self._getVariable(line, fd)
                    self.srcs.extend(var[1].split())
                # match <target>_OBJECTS
                match = self.re_obj.match(line)
                if match and match.group(1) in binTargets_cn:
                    var=self._getVariable(line, fd)
                    self.objs.extend(var[1].split())

            # in second scan get targets for sources
            # with automake sources may lie in another directory
            srcDir = os.path.join(self.sourceDir, self.fileDir)
            fd.seek(0)
            iter_ = iter(fd)
            for line in iter_:
                match = self.re_target.match(line)
                if match:
                    t,s = match.groups()
                    if s in self.srcs:
                        self.targets[t] = createTarget(t,s, targetDir=self.fileDir, sourceDir=srcDir)

        finally:
            fd.close()

    def _getVariable(self, line, iter_):
        "Scan, possibly multiline, variable definition in Makefile."

        match = self.re_var.match(line)
        (name, matchStr) = match.groups()
        value = ""

        while match:
            if line[match.end()-2]=='\\':
                value = " ".join((value, matchStr[:-1].strip()))
                line = iter_.next()
                match = self.re_varCont.match(line)
                if not match:
                    raise RuntimeError("Error reading variable")
                matchStr = match.group(1)
            else:
                value = " ".join((value, matchStr.strip()))
                match = None
        
        return (name, value)

    def __str__(self):
        import StringIO
        io = StringIO.StringIO()
        io.write("Makefile collector: %s, %s\n" % (self.fileName,self.binTarget))
        for tkey in self.targets.iterkeys():
            io.write('%s = %s\n' % ( tkey, (self.targets[tkey] or '')))

        s = io.getvalue()
        io.close()
        return s

# General dependency resolver
class DependencyResolver:
    """Takes logical dependencies from collector, calculates file
    dependencies and writes to file."""

    def __init__(self, collector, collectors = []):
        self.collector = collector
        self.collectors = collectors
        self.deps = {}

        self._calculateDependencies()

    def _calculateDependencies(self):
        "Take collector and fill deps dictionary."
        # first, what files provide what logical names
        # this scan is made for all collectors
        info("Resolving dependencies")
        
        logprovs = {}
        collectors = [self.collector]
        collectors.extend(self.collectors)
        for collector in collectors:
            for t in collector.targets.itervalues():
                if not t: # t may be None, if appropriate Target class not found
                    continue
                for prov in t.getProvs():
                    logprovs[prov] = t

        # now iterate all targets and construct their deps
        # this operation is made only for main Makefile
        for t in self.collector.targets.itervalues():
            if not t:
                continue
            deps = []
            for dep in t.getDeps():
                if logprovs.has_key(dep):
                    prov = logprovs[dep]
                    deps.append(os.path.join(prov.targetDir, prov.targetName))
                else:
                    warn("Unresolved dependency %s for %s" %
                         (dep, t.targetName))
            self.deps[t.targetName] = deps

    def save(self, fileName):
        info("Creating %s" % fileName)
        f = open(fileName, 'w')
        try:
            for tName in self.deps.keys():
                if len(self.deps[tName]) == 0:
                    continue
                f.write('%s:' % tName)
                for dep in self.deps[tName]:
                    f.write(' %s' % dep)
                f.write('\n')    
        except:
            f.close()

if __name__=="__main__":
    import getopt
    #os.chdir("../src/main") # temporary hack for testing only

    amc_names = []

    (opts, args) = getopt.getopt(sys.argv[1:], 'i:', ['srcdir='])
    debug(opts)
    debug(args)

    srcDir = ''
    for op, value in opts:
        if op=='-i':
            amc_names.append(value)
        elif op=='--srcdir':
            srcDir = value
    
    amc = AutomakeCollector(sourceDir=srcDir)
    amcs = []
    for amc_name in amc_names:
        amcs.append(AutomakeCollector(amc_name, sourceDir=srcDir))
    
    resolver = DependencyResolver(amc, amcs)
    resolver.save('Make.deps')

    
