# UI
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import Pmw

import doug.config

import logging
LOG = logging.getLogger(__name__)

class ConfigPanel(tkSimpleDialog.Dialog):
    
    def __init__(self, parent, config, sectionNames=None, readonly=False, *args, **kargs):
        self.config=config
        self.readonly=readonly
        self.sectionNames = sectionNames
        self.fields = {}
        self.done=False
        
        tkSimpleDialog.Dialog.__init__(self, parent, *args, **kargs)
        

    def body(self, parent):
        parent.pack_configure(fill=BOTH, expand=True)
        notebook = Pmw.NoteBook(parent)
        balloon=Pmw.Balloon(parent)
        sectionNames=self.sectionNames
        
        # prepare docs
        sections = {}
        subsections = {}
        for key, docs in doug.config.configDesc.docs.items():
            section = docs['__section__']
            if not sections.has_key(section):
                sections[section] = {}
                subsections[section] = set()
            sections[section][key]=docs

            if docs.has_key('subsection'):
                subsections[section].add(docs['subsection'])
                
        # create widgets
        if not sectionNames:
            sectionNames = ['DEFAULT']
            sectionNames.extend(self.config.sections())
        # for all sections, create a tab
        for sectionName in sectionNames:
            sectionFrame=notebook.add(sectionName)

            subsectionFrames={}
            if subsections.has_key(sectionName) and subsections[sectionName]:
                subsectionNotebook=Pmw.NoteBook(sectionFrame)
                subsectionNotebook.pack(expand=True, fill=BOTH)

                # create frames, options go to these frames
                sectionFrame=subsectionNotebook.add('COMMON') # here go options without subsection
                sectionFrame.grid_columnconfigure(0,weight=0)
                sectionFrame.grid_columnconfigure(1,weight=1)

                for subsection in subsections[sectionName]:
                    subsectionFrame=subsectionNotebook.add(subsection)
                    subsectionFrame.grid_columnconfigure(0,weight=0)
                    subsectionFrame.grid_columnconfigure(1,weight=1)                    
                    subsectionFrames[subsection]=subsectionFrame
            else:
                # there are no subsections, do not create internal notebook
                sectionFrame.grid_columnconfigure(0,weight=0)
                sectionFrame.grid_columnconfigure(1,weight=1)

            if self.config.has_section(sectionName) or sectionName=='DEFAULT':
                values=self.config.items(sectionName, raw=True, nodefaults=True)
            else:
                values=[]
            values.sort()
            row=0
            for key, value in values:
                docs=sections.get(sectionName, {}).get(key, {})
                
                # choose frame to add to
                if docs.has_key('subsection'):
                    frame=subsectionFrames[docs['subsection']]
                else:
                    frame=sectionFrame

                # create label
                label=Label(frame, text=key, anchor=N)
                label.grid(row=row, column=0, sticky=W)

                # create field
                if docs.get('type', None)=='directory':
                    field=ConfigPanel.DirectoryField(self, frame, value)
                elif docs.get('type', None)=='list':
                    field=ConfigPanel.OptionField(self, frame, value, docs)
                else:
                    field=ConfigPanel.StringField(self, frame, value)
                self.fields[(sectionName, key)] = field
                field.grid(row=row, column=1, sticky=EW)

                balloon.bind(label, docs.get('description', ""))
                
                row=row+1

        notebook.pack(fill=BOTH, expand=True)

    class StringField(Entry):
        def __init__(self, panel, parent, value):
            Entry.__init__(self,parent)

            self.panel = panel
            self.insert(END, value)
            if panel.readonly:
                self.config(state='readonly')

    class OptionField(Frame):
        def __init__(self, panel, parent, value, docs):
            Frame.__init__(self, parent)
        
            var = StringVar()
            params=docs['type-params']
            params=params.split(',')
            params=map(lambda v: tuple(map(str.strip, v.split(':'))),
                       params)
            self.var = var
            self.params=params

            self.om=OptionMenu(self, var, *map(self._tostring, params))
            self.set(value)
            if panel.readonly:
                self.om.config(state='disabled')
            
            self.om.pack(side=LEFT)

        def set(self, value):
            values=map(lambda e: e[0], self.params)
            index=values.index(str(value))
            self.var.set(index!=-1 and self._tostring(self.params[index]) or '')

        def get(self):
            v=self.var.get()
            return self._fromstring(v)
            
        def _tostring(self, v):
            if len(v)>1:
                return "%s (%s)" % v
            else:
                return v[0]

        def _fromstring(self, s):
            if '(' in s:
                return s[0:s.index('(')].strip()
            else:
                return s

    class DirectoryField(Frame):
        def __init__(self, panel, parent, value):
            Frame.__init__(self, parent)
            self.grid_columnconfigure(0, weight=1)
        
            self.entry = Entry(self)
            self.entry.insert(0, value)
            self.entry.grid(row=0, column=0, sticky=EW)

            if panel.readonly:
                self.entry.config(state='readonly')
                return

            def _openDirectory():
                d = tkFileDialog.askdirectory()
                if d:
                    self.entry.delete(0, END)
                    self.entry.insert(0, d)
                    
            button = Button(self, text='select')
            button.grid(row=0, column=1)
            button['command'] = _openDirectory

        def get(self):
            return self.entry.get()

    def apply(self):
        if self.readonly:
            return

        try:
            for (sectionName, key), field in self.fields.items():
                self.config.set(sectionName, key, field.get())
            self.done=True
        except Exception, e:
            LOG.error(e, exc_info=e)
        
