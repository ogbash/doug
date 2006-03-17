import re, os, sys

re_hasmod = re.compile("^\s*(!?)\s*module\s+(\S*)")
re_usesmod = re.compile("^\s*(!?)\s*use\s+([^,\s]*),*")

def setext(filename, newext):
    ind = filename.rfind(".")
    if (ind != -1):
        return filename[:ind] + newext

def getext(filename):
    ind = filename.rfind(".")
    if (ind != -1):
        return filename[ind:]
    else:
        return ""

def make_rel_path(isrc, idest):
    src = os.path.normpath(isrc)
    dest = os.path.normpath(idest)
    bsrc = src + "/"
    commpref = os.path.commonprefix([bsrc, dest])
    commlen = len(commpref)
    nwsrc = bsrc[commlen:]
    nwdest = dest[commlen:]
    level = nwsrc.count("/")
    prefix = "../"*level
#    print src, dest, nwsrc, nwdest, prefix, os.path.normpath(prefix + nwdest)
    return os.path.normpath(prefix + nwdest)
        
                
class depgen:
    def getprovdep(self, fname):
        f = open(fname, 'r')
        lhas = []
        luses = []
        for line in f.readlines():
            line = line.lower()
            
            hasmod = re_hasmod.match(line)
            
            if hasmod:
                if hasmod.group(1) != '!':
                    lhas.append(hasmod.group(2))
            
            usesmod = re_usesmod.match(line)
            if usesmod:
                if (usesmod.group(1) != '!') and (usesmod.group(2) not in luses):
                    luses.append(usesmod.group(2))
#        print fname, 'has', lhas, 'uses', luses
        return (lhas, luses)

    provs = {}
    deps = {}
    
    def parsefiles(self, dirn):    
        print 'Parsing files: ',
        for fname in os.listdir(dirn):
            if not fname.endswith('.f90'): continue
            lhas, luses = self.getprovdep(fname)
            for i in lhas:
                self.provs[i] = fname
            for i in luses:
                self.deps[fname] = luses
        print 'done'
    
    fdeps = {}
    fdeps_done = 0    
    
    def makedeps(self, save=0):
        print 'Building dependencies: ',
        if save: depfile = open(save, 'w')
        for fname in self.deps.keys():
#            print fname
            unresolved = []
            if save: depfile.write(setext(fname, '.o') + ':')
            plname = setext(fname, '')
            if not plname in self.fdeps: self.fdeps[plname] = []
            for d in self.deps[fname]:
                if d in self.provs:
                    if save: depfile.write(' ' + setext(self.provs[d], '.o'))
                    self.fdeps[plname].append(setext(self.provs[d], ''))
                else:              
                    unresolved.append(d)
            depfile.write('\n')
            if unresolved and self.verbose:
                print 'Unresolved deps for %s:' % fname,
                for n in unresolved: print '"%s"' % n,
                print 
        self.fdeps_done = 1
        print 'done'
    
    def alldeps(self, infname):
        #print self.fdeps
        if not self.fdeps_done:
            self.makedeps()
        fname = os.path.normpath(infname)
        fdp = [fname]
        idp = self.fdeps[fname]
        while(idp):
            x = idp.pop(0)
            if (x not in fdp):
                fdp.append(x)
                if x in self.fdeps:
                    for j in self.fdeps[x]:
                        idp.append(j)

        dirn = os.path.dirname(fname)
        
        relfdp = []
        
        for p in fdp:
            relfdp.append(make_rel_path(dirn, p))
                                            
        return relfdp
        
    def getsrcs(self, fname):
        r = self.alldeps(fname)
        s = ''    
        for i in r:
            s += ' ' + i + '.o'
        return s
        
    def parse_recursive(self, topdir):
        if self.verbose: print "Parsing files:"
        for (dirpath, dirnames, filenames) in os.walk(topdir, topdown=True):
            # Remove unwanted dirs so walk() wont enter them later
            for dirn in dirnames:  
#                if os.path.basename(dirn) == ".svn":
                if dirn[:1] == ".":
                    dirnames.remove(dirn)
            for filen in filenames:
                if getext(filen).lower() == ".f90": 
                    fullname = os.path.join(dirpath, filen)
                    if self.verbose: print fullname,
                    lhas, luses = self.getprovdep(fullname)
                    if lhas and self.verbose: print "provides:",
                    for i in lhas:
                        self.provs[i] = fullname
                        if self.verbose: print "'%s'" % i,
                    if luses and self.verbose: print "requires:",
                    self.deps[fullname] = luses
                    if self.verbose:
                        for i in luses:
                            print "'%s'" % i,
                    if self.verbose: print
#            print dirpath
    
    def write_depends(self, topdir, depfilename):        
        for (dirpath, dirnames, filenames) in os.walk(topdir, topdown=True):
            # Assume the dir is irrelevant if it contains no makefile
            if "Makefile" not in filenames: continue  
            self.makedeps_adv(dirpath, depfilename)
            
#            print dirpath, level


    def makedeps_adv(self, dirname, depfilename):
        fulldepfilename = os.path.join(dirname, depfilename)
        print "Creating", os.path.normpath(fulldepfilename)
        depfile = open(fulldepfilename, 'w')
        for fname in self.deps.keys():
            if os.path.dirname(fname) != dirname: continue 
            if self.verbose > 1: print fname, ":", 
            unresolved = []
            depfile.write(setext(os.path.basename(fname), '.o') + ':')
            
            plname = os.path.normpath(setext(fname, ''))
            if not plname in self.fdeps: self.fdeps[plname] = []
            for d in self.deps[fname]:
                if d in self.provs:
                    if self.verbose > 1: print make_rel_path(dirname, self.provs[d]),
                    depfile.write(' ' + make_rel_path(dirname, setext(self.provs[d], '.o')))
                    self.fdeps[plname].append(os.path.normpath(setext(self.provs[d], '')))
                else:              
                    unresolved.append(d)
            depfile.write('\n')
            if self.verbose > 1: print
            self.fdeps_done = 1
#             if unresolved:
#                 print 'Unresolved deps for %s:' % fname,
#                 for n in unresolved: print '"%s"' % n,
#                 print 
#         print 'done'
        
#        depfile = open(, 'w')
 #       for fname in self.deps.keys():
#            print fname        pass
        
        
    def __init__(self, dirn, be_recursive = False, verbose = 0):
        self.verbose = verbose
        if be_recursive:
            self.parse_recursive(dirn)
        else:
            self.parsefiles(dirn)
        
