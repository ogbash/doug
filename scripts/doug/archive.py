import copy
import os

from doug.config import DOUGConfigParser

import logging
LOG = logging.getLogger(__name__)

class Archive(object):
    DIRTY = 'DIRTY'
    SAVED = 'SAVED'
    ARCHIVED = 'ARCHIVED'

    FILETYPES = ['', 'Grid', 'Matrix', 'Vector', 'Vector/RHS', 'Vector/Solution',
                 'Aggregates/Fine', 'Aggregates/Coarse', 'Text', 'Text/Output', 'Text/Error']

    def _setState(self, newstate): self._state=newstate
    state = property(lambda self: self._state, fset=_setState)
    
    def __init__(self, name=None, directoryName=None, archiveType='problem'):
        self.name = name
        self.filetypes = {} #: file name -> type
        
        if directoryName==None:
            directoryName = os.path.abspath(name)
            
        if os.path.isdir(directoryName):
            self._state = Archive.SAVED
            self.directoryName = directoryName
            self.load()
        else:
            os.mkdir(directoryName)
            self._state = Archive.DIRTY
            self.directoryName = directoryName

            self.info = DOUGConfigParser()
            self.info.add_section('general')
            self.info.set('general', 'archive-type', archiveType)
            self.info.add_section('files')

    def save(self):
        if self.state == Archive.DIRTY:
            LOG.info("Saving '%s' to %s" % (self.name, self.directoryName))
            config = copy.deepcopy(self.info)
            
            config.set('general', 'name', self.name)
            configFile = file(os.path.join(self.directoryName, '.info'), 'w')

            files=map(lambda f: "%s:%s"%f, self.filetypes.items())
            config.set('files', 'types', ", ".join(files))
            
            config.write(configFile)
            configFile.close()

    def load(self):
        LOG.info("Loading from %s" % self.directoryName)
        config = DOUGConfigParser()
        config.read([os.path.join(self.directoryName, '.info')])
        self.info = config

        if config.has_section('general'):
            self.name = config.get('general', 'name')
        else:
            config.add_section('general')

        if config.has_section('files'):
            types = config.get('files', 'types', default="")
            types = types.split(",")
            types = filter(bool, types) # filter out empty strings
            types = map(lambda s: map(str.strip, s.split(':')), types)
            
            for filename, filetype in types:
                self.setFileType(filename, filetype)
        else:
            config.add_section('files')

    def close(self):
        self.save()

    def setFileType(self, filename, filetype):
        if filetype==None:
            if self.filetypes.has_key(filename):
                del self.filetypes[filename]
        else:
            self.filetypes[filename] = filetype

        if filetype=='Text/Profile':
            self._readProfileFile(filename)
        
        self.state = Archive.DIRTY

    def _readProfileFile(self, filename):
        f = open(os.path.join(self.directoryName, filename))
        try:
            for line in f:
                line = line.strip()
                proc, name, value = line.split(':')
                name = name.replace(' ', '-')
                self.info.set('doug-profile', name, value, addsection=True)
        finally:
            f.close()

    def getFileType(self, filename):
        return self.filetypes.get(filename, None)

    def getFiles(self, filetype):
        files = map(lambda g: g[0],
                    filter(lambda x: x[1]==filetype,
                           self.filetypes.items()))
        return files

    def __str__(self):
        return "Archive(directory=%s)" % self.directoryName

