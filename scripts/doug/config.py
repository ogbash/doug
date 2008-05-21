
from ConfigParser import SafeConfigParser
import os.path
import re
import StringIO
import logging
logging.basicConfig(level=logging.DEBUG/2)
LOG = logging.getLogger(__name__)

class DOUGConfigParser(SafeConfigParser):
    def __init__(self, *args, **kargs):
        self.__name = ''
        if 'name' in kargs:
            self.__name = kargs.get('name')
            del kargs['name']
        self.description=ConfigDesc("")
        SafeConfigParser.__init__(self, *args, **kargs)

    def getpath(self, section, option):
        value = self.get(section, option)
        if self.has_option(section, 'wd'):
            basepath = self.get(section, 'wd')
            return os.path.join(basepath, value)
        else:
            return value

    def items(self, section, *args, **kargs):
        nodefaults = kargs.get('nodefaults', None)
        if nodefaults!=None:
            del kargs['nodefaults']
        lst = SafeConfigParser.items(self, section, *args, **kargs)
        if nodefaults:
            nlst = []
            for k,v in lst:
                if self._sections[section].has_key(k):
                    nlst.append((k,v))
            return nlst

        return lst

    def addConfig(self, conf):
        LOG.debug("Adding %s to %s" % (conf, self))
        
        def copyOptions(sectionName, section):
            for option in section.keys():
                value = section[option]
                LOG.log(logging.DEBUG/2, "Setting %s in %s to %s", option, sectionName, value)
                self.set(sectionName, option, value)

        # copy defaults
        LOG.debug("Copying default options")
        for option, value in conf._defaults.items():
            LOG.log(logging.DEBUG/2, "Setting %s in DEFAULT to %s", option, value)
            self._defaults[option]=value

        # copy options
        for sectionName, section in conf._sections.items():
            if not self.has_section(sectionName):
                self.add_section(sectionName)
            LOG.debug("Copying options of section '%s'", sectionName)
            copyOptions(sectionName, section)

    def addConfigContents(self, contents):
        self.readfp(StringIO.StringIO(contents))
        LOG.debug("Adding description to %s" % (self,))
        self.description.addDescription(contents)

    def __str__(self):
        return "Config<%s>" % (self.__name)

class ConfigDesc:
    DOCS_RE = re.compile('^#\s*(?:@(\w+))?\s*:(.*)')
    VALS_RE = re.compile('^\s*([^:=\s]+)\s*[:=].*')
    SECT_RE = re.compile('^\s*\[\s*([^\] \t\n]+)\s*\]')
    
    def __init__(self, lines):
        self.docs = {}
        self.__readDesc(lines)

    def addDescription(self, contents):
        self.__readDesc(contents.split('\n'))

    def __readDesc(self, lines):
        doc = {}
        section = None
        
        for line in lines:
            match = ConfigDesc.SECT_RE.match(line)
            if match:
                section, = match.groups()
                continue
            
            match = ConfigDesc.DOCS_RE.match(line)
            if match:
                key, value = match.groups()
                key = key or 'description'
                doc[key] = value.strip()
                continue
            
            match = ConfigDesc.VALS_RE.match(line)
            if match:
                key, = match.groups()
                doc['__section__'] = section
                LOG.log(logging.DEBUG/3, "Found description %s: %s" % (key, doc))
                self.docs[key] = doc
                doc = {}
                
    def __str__(self):
        s=StringIO.StringIO()
        s.write("ConfigDesc{")
        kvs=", ".join(map(
            lambda item: "%s=%s" % item,
            self.docs.items()))
        s.write(kvs)
        s.write("}")
        return s.getvalue()

from Tkinter import *
import Pmw
class ConfigPanel:
    def __init__(self, parent, config, sectionNames=None, readonly=False):
        self.frame=Toplevel(parent)
        self.config=config
        self.readonly=readonly

        self._showConfig(self.frame, sectionNames)
        
        b = Button(self.frame, text="OK", command=self.ok)
        b.pack(pady=5)
        self.frame.grab_set()
        parent.wait_window(self.frame)

    def _showConfig(self, parent, sectionNames):
        notebook = Pmw.NoteBook(parent)
        balloon=Pmw.Balloon(parent)
        
        # prepare docs
        sections = {}
        for key, docs in self.config.description.docs.items():
            section = docs['__section__']
            if not sections.has_key(section):
                sections[section] = {}
            sections[section][key]=docs

        # create widgets
        for sectionName in sectionNames or self.config.sections():
            sectionFrame=notebook.add(sectionName)
            #sectionFrame = LabelFrame(parent, text=section)
            sectionFrame.grid_columnconfigure(0,weight=0)
            sectionFrame.grid_columnconfigure(1,weight=1)

            values=self.config.items(sectionName, nodefaults=True)
            values.sort()
            row=0
            for key, value in values:
                docs=sections.get(sectionName, {}).get(key, {})
                
                label=Label(sectionFrame, text=key, anchor=N)
                label.grid(row=row, column=0, sticky=W)
                field=Entry(sectionFrame)
                field.insert(END, value)
                if self.readonly:
                    field.config(state='readonly')
                field.grid(row=row, column=1, sticky=EW)

                balloon.bind(label, docs.get('description', ""))
                
                row=row+1

        notebook.pack(fill=BOTH, expand=True)

    def ok(self):
        self.frame.destroy()
