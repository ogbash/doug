from Tkinter import *
import Pmw

import weakref
import time
import os
import copy

from doug.config import DOUGConfigParser
from doug.configui import ConfigPanel
import doug.execution
import gridplot

import logging
LOG = logging.getLogger('app')

class ArchivePanel:
    def __init__(self,parent,globalConfig, app):
        self.archive = None
        self.app = weakref.proxy(app)
        self.globalConfig = globalConfig
        
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

        panel = Pmw.ScrolledFrame(filesFrame)
        panel.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        self.filesPanel = panel.interior()
        self.filesPanel.grid_columnconfigure(2, weight=1) # extend last (empty) column to align others to the right
        
    def setArchive(self, archive):
        LOG.debug("Showing %s" % archive)
        self.archive = archive
        
        if archive==None:
            self._setEntryField(self.nameField, '')
            self._setEntryField(self.directoryField, '')
            self.stateField['text'] = ''
            self._gridClean(self.filesPanel)
            return
        
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

            view=Button(self.filesPanel, text="view",command=ArchivePanel._ViewFile((self.app,archive,fileName)))
            view.grid(row=i,column=2)

    class _ViewFile:
        def __init__(self, params):
            self.params=params

        def __call__(self):
            app, archive, fileName = self.params
            fullpath = os.path.join(archive.directoryName, fileName)
            viewer=app.config.get('global', 'viewer')
            os.spawnlp(os.P_NOWAIT, viewer, viewer, fullpath)
        

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

        try:
            config = DOUGConfigParser(name='DOUG execution parameters')
            config.addConfigContents(doug.execution.getDefaultConfigContents())
            config.addConfig(self.globalConfig)            
            
            suffix = time.strftime('--run-%y%m%d-%H%M%S')
            config.set('DEFAULT', 'doug-workdir', self.archive.directoryName+suffix)
            config.set('DEFAULT', 'doug-datadir', self.archive.directoryName)

            matrixFiles = self.archive.getFiles('Matrix')
            if matrixFiles:
                matrixFile='%(doug-datadir)s/'+matrixFiles[0]
            else:
                matrixFile=''
            config.set('doug-controls', 'assembled_mtx_file', matrixFile)
                           
            rhsFiles = self.archive.getFiles('Vector/RHS')
            if rhsFiles:
                rhsFile='%(doug-datadir)s/'+rhsFiles[0]
            else:
                rhsFile=''
            config.set('doug-controls', 'assembled_rhs_file', rhsFile)

            # change configurations
            configPanel = ConfigPanel(self.frame, config, title="DOUG configuration",
                                      sectionNames=['DEFAULT', 'doug', 'doug-controls'])
            if not configPanel.done:
                return
            
            # run            
            execution = doug.execution.DOUGExecution(config)
            try:
                archive = Archive(self.archive.name+suffix,
                                  directoryName=self.archive.directoryName+suffix,
                                  archiveType='solution')
                archive.info.set('general','problem-name',self.archive.name)
                archive.info.addConfig(execution.config)
                resultConfig=execution.run()
                archive.info.addConfig(resultConfig)

                if archive.info.has_option('doug-result', 'solutionfile'):
                    archive.setFileType(archive.info.get('doug-result', 'solutionfile'), 'Vector/Solution')
                if archive.info.has_option('doug-result', 'errorfile'):
                    archive.setFileType(archive.info.get('doug-result', 'errorfile'), 'Text/Error')
                if archive.info.has_option('doug-result', 'outputfile'):
                    archive.setFileType(archive.info.get('doug-result', 'outputfile'), 'Text/Output')
                if archive.info.has_option('doug-result', 'fineaggrsfile'):
                    archive.setFileType(archive.info.get('doug-result', 'fineaggrsfile'), 'Aggregates/Fine')
                if archive.info.has_option('doug-result', 'coarseaggrsfile'):
                    archive.setFileType(archive.info.get('doug-result', 'coarseaggrsfile'), 'Aggregates/Coarse')
                    
            finally:
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
            device=self.app.config.get('global', 'plplot-device')
            plot.run(dev=device)

    def _showRHS(self):
        if self.archive==None:
            return
        
        gridFiles = self.archive.getFiles(filetype='Grid')
        solutionFiles = self.archive.getFiles(filetype='Vector/RHS')
        if gridFiles and solutionFiles:
            gridFile = os.path.join(self.archive.directoryName, gridFiles[0])
            solutionFile = os.path.join(self.archive.directoryName, solutionFiles[0])
            plot = gridplot.Plot(gridFile, solutionFile=solutionFile)
            device=self.app.config.get('global', 'plplot-device')
            plot.run(dev=device)


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
                                readonly=True, sectionNames=['DEFAULT', 'doug', 'doug-controls', 'doug-result'],
                                title='Execution info')
        
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
            device=self.app.config.get('global', 'plplot-device')
            plot.run(dev=device)
            
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
            device=self.app.config.get('global', 'plplot-device')
            plot.run(dev=device)

    def setArchive(self, archive):
        ArchivePanel.setArchive(self, archive)

        if archive==None:
            self._setEntryField(self.problemNameField, '')
            return
        
        problemName = archive.info.get('general', 'problem-name')
        self._setEntryField(self.problemNameField, problemName)
        

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
            types = config.get('files', 'types', "").split(",")
            types = filter(bool, types) # filter out empty strings
            types = map(lambda s: map(str.strip, s.split(':')), types)
            
            for filename, filetype in types:
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

    def __str__(self):
        return "Archive(directory=%s)" % self.directoryName

