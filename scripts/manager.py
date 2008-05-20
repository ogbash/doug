#! /usr/bin/env python

# import packages from modules directory
import sys, os
pathname = os.path.dirname(sys.argv[0])
sys.path.append(os.path.join(os.path.abspath(pathname), 'modules'))

from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox
import Pmw
import random
import os
import shutil
import os.path
from ConfigParser import SafeConfigParser
from StringIO import StringIO
import copy
import gridplot
import doug.execution
from doug.config import ConfigDesc, DOUGConfigParser, ConfigPanel
import weakref

import logging
logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
logging.addLevelName(5, 'FINE')
logging.addLevelName(3, 'FINER')
logging.addLevelName(2, 'FINEST')
logging.addLevelName(1, 'TRACE')

LOG = logging.getLogger('app')

_configContent = """
[DEFAULT]
#: DOUG binary directory
#@type: directory
doug-bindir: /usr/bin

"""

#_dougConfigContent = """
#[doug]
#"""

class ArchivePanel:
    def __init__(self,parent,globalConfig, app):
        self.archive = None
        self.app = weakref.proxy(app)
        self.globalConfig = globalConfig
        #self.dougConfig = SafeConfigParser()
        #self.dougConfig.readfp(StringIO(_dougConfigContent))
        
        self.frame = LabelFrame(parent, text="Archive")
        infoFrame = Frame(self.frame)
        self.infoFrame = infoFrame
        row = 0
        
        label = Label(infoFrame, text="name")
        label.grid(row=row)
        self.nameField = Entry(infoFrame, state='readonly', width=40)
        self.nameField.grid(row=row, column=1,sticky=W)
        row=row+1

        label = Label(infoFrame, text="directory")
        label.grid(row=row)
        self.directoryField = Entry(infoFrame, state='readonly', width=80)
        self.directoryField.grid(row=row,column=1,sticky=W)
        row=row+1

        label = Label(infoFrame, text="state")
        label.grid(row=row)
        self.stateField = Label(infoFrame)
        self.stateField.grid(row=row, column=1,sticky=W)
        row=row+1
        
        infoFrame.pack(side=TOP, anchor=W, padx=5, pady=5, ipady=5)

        buttonFrame = Frame(self.frame)
        buttonFrame.pack(side=TOP, fill=X)
        self.buttonFrame = buttonFrame

        filesFrame = LabelFrame(self.frame, text="Files")
        filesFrame.pack(side=TOP, expand=True, fill=BOTH)

        filesButtonPanel = Frame(filesFrame)
        filesButtonPanel.pack(side=BOTTOM, fill=X)
        self.updateButton = Button(filesButtonPanel, text="Update", command=self._update, fg="green")
        self.updateButton.pack(side=LEFT)
        self.resetButton = Button(filesButtonPanel, text="Reset", fg="red")
        self.resetButton.pack(side=LEFT)
        
        self.filesPanel = Frame(filesFrame)
        self.filesPanel.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        self.filesPanel.grid_columnconfigure(2, weight=1) # extend last (empty) column to align others to the right
        
    def setArchive(self, archive):
        LOG.debug("Showing %s" % archive)
        self.archive = archive
        self._setEntryField(self.nameField, archive.name)
        self._setEntryField(self.directoryField, archive.directoryName)
        self.stateField['text'] = archive.state

        self._gridClean(self.filesPanel)
        files = os.listdir(archive.directoryName)
        for i,fileName in enumerate(files):
            l = Label(self.filesPanel, text=fileName)
            l.grid(row=i,sticky=W)
            var = StringVar()
            lb = OptionMenu(self.filesPanel, var, *Archive.FILETYPES)
            lb.var = var
            filetype = archive.getFileType(fileName)
            var.set(filetype or '')
            lb.grid(row=i,column=1,sticky=E)

    def _setEntryField(self, field, value):
        state = field['state']
        field.config(state=NORMAL)
        field.delete(0, END)
        field.insert(END, value)
        field.config(state=state)

    def _setTextField(self, field, value):
        state = field['state']
        field.config(state=NORMAL)
        field.delete('0.0', END)
        field.insert(END, value)
        field.config(state=state)

    def _gridClean(self, grid):
        slaves = grid.grid_slaves()
        for slave in slaves:
            slave.grid_forget()
            slave.destroy()

    def _update(self):
        archive = self.archive
        if archive==None:
            return
        
        ncols,nrows = self.filesPanel.size()
        for i in xrange(nrows):
            slaves = self.filesPanel.grid_slaves(row=i)
            slaves.reverse() # because it starts from the last
            filename = slaves[0]['text']
            filetype = slaves[1]['text']
            archive.setFileType(filename, filetype or None)


class ProblemArchivePanel(ArchivePanel):
    def __init__(self, *args, **kargs):
        ArchivePanel.__init__(self, *args, **kargs)

        buttonFrame = self.buttonFrame
        self.showGridButton = Button(buttonFrame, text="Show grid", command=self._showGrid)
        self.showGridButton.pack(side=RIGHT)
        self.showGridButton = Button(buttonFrame, text="Show RHS", command=self._showRHS)
        self.showGridButton.pack(side=RIGHT)
        self.solveButton = Button(buttonFrame, text="Solve", command=self._solve, fg='blue')
        self.solveButton.pack(side=RIGHT)        
    
    def _solve(self):
        if self.archive==None:
            return

        from doug.execution import DOUGExecution
        import time
        try:
            config = DOUGConfigParser(name='DOUG execution parameters')
            config.addConfigContents(doug.execution.getDefaultConfigContents())
            config.addConfig(self.globalConfig)            
            
            suffix = time.strftime('_run-%y%m%d%H%M%S')
            config.set('doug', 'cwd', self.archive.directoryName+suffix)
            config.set('DEFAULT', 'doug-datadir', self.archive.directoryName)

            # change configurations
            configPanel = ConfigPanel(self.frame, config)
            
            # run            
            execution = doug.execution.DOUGExecution(config)
            archive = Archive(self.archive.name+suffix,
                              directoryName=self.archive.directoryName+suffix,
                              archiveType='solution')
            archive.info.set('general','problem-name',self.archive.name)
            archive.info.addConfig(execution.config)
            execution.run()
            self.app.addArchive(archive)
        except(Exception), e:
            LOG.error(e, exc_info=e)

    def _showGrid(self):
        if self.archive==None:
            return

        gridFiles = self.archive.getFiles(filetype='Grid')
        if gridFiles:
            gridFile = os.path.join(self.archive.directoryName, gridFiles[0])
            plot = gridplot.Plot(gridFile)
            plot.run(dev='xwin')

    def _showRHS(self):
        if self.archive==None:
            return
        
        gridFiles = self.archive.getFiles(filetype='Grid')
        solutionFiles = self.archive.getFiles(filetype='Vector/RHS')
        if gridFiles and solutionFiles:
            gridFile = os.path.join(self.archive.directoryName, gridFiles[0])
            solutionFile = os.path.join(self.archive.directoryName, solutionFiles[0])
            plot = gridplot.Plot(gridFile, solutionFile=solutionFile)
            plot.run('xwin')


class SolutionArchivePanel(ArchivePanel):
    def __init__(self, *args, **kargs):
        ArchivePanel.__init__(self, *args, **kargs)

        buttonFrame = self.buttonFrame
        self.showGridButton = Button(buttonFrame, text="Show solution", command=self._showSolution)
        self.showGridButton.pack(side=RIGHT)
        self.showGridButton = Button(buttonFrame, text="Show fine aggrs", command=self._showFineAggregates)
        self.showGridButton.pack(side=RIGHT)
        self.showGridButton = Button(buttonFrame, text="Show coarse aggrs", command=self._showCoarseAggregates)
        self.showGridButton.pack(side=RIGHT)
        self.showGridButton = Button(buttonFrame, text="Show info", command=self._showInfo)
        self.showGridButton.pack(side=RIGHT)

        row,column = self.infoFrame.size()
        row = row+1

        label = Label(self.infoFrame, text="problem name")
        label.grid(row=row)
        self.problemNameField = Entry(self.infoFrame, state='readonly', width=40)
        self.problemNameField.grid(row=row, column=1, sticky=W)

    def _showInfo(self):
        configPanel=ConfigPanel(self.frame, self.archive.info,
                                readonly=True, sectionNames=['doug', 'doug-controls'])
        
    def _showSolution(self):
        problemName = self.archive.info.get('general','problem-name')
        problemArchive = self.app.getProblemArchive(problemName)
        if problemArchive==None:
            return
            
        gridFiles = problemArchive.getFiles(filetype='Grid')
        solutionFiles = self.archive.getFiles(filetype='Vector/Solution')
        if gridFiles and solutionFiles:
            gridFile = os.path.join(problemArchive.directoryName, gridFiles[0])
            solutionFile = os.path.join(self.archive.directoryName, solutionFiles[0])
            plot = gridplot.Plot(gridFile, solutionFile=solutionFile)
            plot.run('xwin')
            
    def _showFineAggregates(self):
        self._showAggregates('Aggregates/Fine')

    def _showCoarseAggregates(self):
        self._showAggregates('Aggregates/Coarse')

    def _showAggregates(self, aggrFileType):
        problemName = self.archive.info.get('general','problem-name')
        problemArchive = self.app.getProblemArchive(problemName)
        if problemArchive==None:
            return
            
        gridFiles = problemArchive.getFiles(filetype='Grid')
        aggrFiles = self.archive.getFiles(filetype=aggrFileType)
        if gridFiles and aggrFiles:
            gridFile = os.path.join(problemArchive.directoryName, gridFiles[0])
            aggrFile = os.path.join(self.archive.directoryName, aggrFiles[0])
            plot = gridplot.Plot(gridFile, aggregateFile=aggrFile)
            plot.run('xwin')

    def setArchive(self, archive):
        ArchivePanel.setArchive(self, archive)

        problemName = archive.info.get('general', 'problem-name')
        self._setEntryField(self.problemNameField, problemName)
        

class Archive(object):
    DIRTY = 'DIRTY'
    SAVED = 'SAVED'
    ARCHIVED = 'ARCHIVED'

    FILETYPES = ['', 'Grid', 'Matrix', 'Vector', 'Vector/RHS', 'Vector/Solution',
                 'Aggregates/Fine', 'Aggregates/Coarse']

    def _setState(self, newstate): self._state=newstate
    state = property(lambda self: self._state, fset=_setState)
    
    def __init__(self, name=None, directoryName=None, archiveType='problem'):
        self.name = name
        self.filetypes = {} #: file name -> type
        
        if directoryName==None:
            directoryName = name
            
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

            for filename, filetype in self.filetypes.items():
                config.set('files', filename, filetype)
            
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
            config.add_section('files')

        if config.has_section('files'):
            for filename, filetype in config.items('files'):
                self.filetypes[filename] = filetype
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
            
        self.state = Archive.DIRTY

    def getFileType(self, filename):
        return self.filetypes.get(filename, None)

    def getFiles(self, filetype):
        files = map(lambda g: g[0],
                    filter(lambda x: x[1]==filetype,
                           self.filetypes.items()))
        return files

class OptionDialog(tkSimpleDialog.Dialog):
    def __init__(self, master, configDesc, config, *args, **kargs):
        self.res = False
        self.configDesc = configDesc
        self.config = config
        tkSimpleDialog.Dialog.__init__(self, master, *args, **kargs)
    
    def body(self, master):
        self.fields = {}

        for i, (key, doc) in enumerate(self.configDesc.docs.items()):
            Label(master, text=doc['description'] or key).grid(row=i)
            field = Entry(master)
            self.fields[key] = field
            value = self.config.get(doc['__section__'], key, raw=True)
            field.insert(0, value)
            field.grid(row=i, column=1)

            if doc.get('type', None)=='directory':
                def _openDirectory():
                    d = tkFileDialog.askdirectory()
                    if d:
                        field.delete(0, END)
                        field.insert(0, d)
                    
                button = Button(master, text='select')
                button.grid(row=i, column=2)
                button['command'] = _openDirectory

    def apply(self):
        self.res = True
        for option in self.fields.keys():
            section = self.configDesc.docs[option]['__section__']
            self.config.set(section, option, self.fields[option].get())

class App:
    def __reload(self, ev):
        print "Reloading execution"
        reload(doug.execution)
    
    def __init__(self, root):
        self.root = root
        self.problemsArchives = []
        self.solutionsArchives = []
        self.__loadOptions()

        root.protocol('WM_DELETE_WINDOW', self.__closeApp)
        root.bind('<Control-R>', self.__reload)
        
        bottomFrame = Frame(root)
        bottomFrame.pack(side=BOTTOM, fill=X)

        quitb = Button(bottomFrame, text="Quit", command=self.__closeApp)
        quitb.pack(side=RIGHT)

        menubar = Menu(root)
        root.config(menu=menubar)

        # main frame
        mainFrame=Frame(root)
        mainFrame.pack(expand=True, fill=BOTH)
        problemsFrame = LabelFrame(mainFrame, text='Problems', fg='blue')
        problemsFrame.pack(side=TOP, expand=True, fill=BOTH)
        solutionsFrame = LabelFrame(mainFrame, text='Solutions', fg='blue')
        solutionsFrame.pack(side=TOP, expand=True, fill=BOTH)

        # problems listbox
        listbox = Pmw.ScrolledListBox(problemsFrame, selectioncommand=self.__problemsListboxSelect)
        self.problemsListbox = listbox
        listbox.pack(side=LEFT, fill=Y)
        #listbox.bind('<<ListboxSelect>>', self.__problemsListboxSelect)
        listbox.bind('<<AddFiles>>', self.__addFiles)
        listbox.event_add('<<AddFiles>>', '<2>')

        addb = Button(problemsFrame, text="Add archive", command=self.__addArchive)
        addb.pack(side=TOP)

        # solutions listbox
        listbox = Pmw.ScrolledListBox(solutionsFrame, selectioncommand=self.__solutionsListboxSelect)
        self.solutionsListbox = listbox
        listbox.pack(side=LEFT, fill=Y)
        #listbox.bind('<<ListboxSelect>>', self.__solutionsListboxSelect)
        
        listbox.bind('<<AddFiles>>', self.__addFiles)
        listbox.event_add('<<AddFiles>>', '<2>')

        self.problemPanel = ProblemArchivePanel(problemsFrame, self.config, self)
        self.problemPanel.frame.pack(side=TOP, expand=True, fill=BOTH)

        self.solutionPanel = SolutionArchivePanel(solutionsFrame, self.config, self)
        self.solutionPanel.frame.pack(side=TOP, expand=True, fill=BOTH)

        menu = Menu(menubar, tearoff=False)
        menubar.add_cascade(label='File', menu=menu)

        menu = Menu(menubar, tearoff=False)
        menu.add_command(label='Options', command=self.__editOptions)
        menubar.add_cascade(label='Tools', menu=menu)

        self.__scanDirectory('.', {
            'problem': (self.problemsArchives, self.problemsListbox),
            'solution': (self.solutionsArchives, self.solutionsListbox)})

    def __editOptions(self):
        d = OptionDialog(self.root, self.configDesc, self.config)
        if d.res:
            self.__saveOptions()

    def __scanDirectory(self, dirpath, archiveTypes):
        dirs = os.listdir(dirpath)
        dirs = map(lambda x: os.path.join(dirpath, x), dirs) # create full paths
        dirs = filter(os.path.isdir, dirs) # filter directories
        
        for subdir in dirs:
            if os.path.exists(os.path.join(subdir, '.info')):
                abssubdir = os.path.abspath(subdir)
                archive = Archive(directoryName = abssubdir)
                t = archive.info.get('general', 'archive-type')
                if t in archiveTypes:
                    archiveTypes[t][0].append(archive)
                    archiveTypes[t][1].insert(END, archive.name)

    def __problemsListboxSelect(self, ev=None):
        widget=self.problemsListbox
        index = int(widget.curselection()[0])
        archive = self.problemsArchives[index]
        self.problemPanel.setArchive(archive)

    def __solutionsListboxSelect(self, ev=None):
        widget=self.solutionsListbox
        index = int(widget.curselection()[0])
        archive = self.solutionsArchives[index]
        self.solutionPanel.setArchive(archive)

    def addArchive(self, archive):
        archiveType = archive.info.get('general', 'archive-type')
        if archiveType=='problem':
            self.problemsArchives.append(archive)
            self.problemsListbox.insert(END, archive.name)
        elif archiveType=='solution':
            self.solutionsArchives.append(archive)
            self.solutionsListbox.insert(END, archive.name)
            
    def __addArchive(self):
        name = tkSimpleDialog.askstring("Insert name", "Type new archive name (without .tar.gz)")
        if name:
            archive = Archive(name)
            self.addArchive(archive)

    def __addFiles(self, ev=None):
        sel = ev.widget.curselection()
        if not sel:
            return

        archive = self.problemsArchives[int(sel[0])]

        filenames = tkFileDialog.askopenfilenames(title="Select files to include",
                                                      filetypes=[('All',"*"), ("Text data","*.txt")])
        if filenames:
            #isMove= tkMessageBox.askyesno(title="%d files selected. Copy or move?" % len(filenames),
            #                              message="Move files?",
            #                              default="no")
            for filename in filenames:
                newfilename = os.path.join(archive.directoryName, os.path.basename(filename))
                shutil.copy(filename, newfilename)

    def __closeApp(self):
        for archive in self.problemsArchives:
            archive.close()
        for archive in self.solutionsArchives:
            archive.close()
        LOG.info("Closing application")
        self.root.destroy()

    def __loadOptions(self):
        self.config = SafeConfigParser()
        self.config.readfp(StringIO(_configContent))
        if os.path.exists('.manager.conf'):
            self.config.read('.manager.conf')
            
        self.configDesc = ConfigDesc(_configContent.split('\n'))
        
        if os.path.isdir('.manager.conf'):
            self.config.read('.manager.conf')

    def __saveOptions(self):
        f = open('.manager.conf', 'w')
        try:
            LOG.info('Saving %s' % f)
            self.config.write(f)
        finally:
            f.close()

    def getProblemArchive(self, name):
        for archive in self.problemsArchives:
            if archive.name==name:
                return archive
        return None

root = Tk()
app = App(root)
root.mainloop()
