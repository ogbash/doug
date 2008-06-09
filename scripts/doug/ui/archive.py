from Tkinter import *
import Pmw

import weakref
import time
import os

from doug.config import DOUGConfigParser

from doug.configui import ConfigPanel
from doug.archive import Archive
import doug.execution
import gridplot

import logging
LOG = logging.getLogger(__name__)

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

        panel = Pmw.ScrolledFrame(filesFrame)
        panel.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        self.filesPanel = panel.interior()
        self.filesPanel.grid_columnconfigure(2, weight=1) # extend last (empty) column to align others to the right
        
    def setArchive(self, archive):
        LOG.debug("Showing %s" % archive)
        self.archive = archive
        self.updateFields()

    def updateFields(self):
        archive = self.archive
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
        for i,filename in enumerate(files):
            l = Label(self.filesPanel, text=filename)
            l.grid(row=i,sticky=W)
            
            var = StringVar()
            lb = OptionMenu(self.filesPanel, var, *Archive.FILETYPES)
            lb.var = var
            filetype = archive.getFileType(filename)
            var.set(filetype or '')
            lb.grid(row=i,column=1,sticky=E)

            if filetype:
                fields=self._getFileSpecificFields(self.filesPanel, filename, filetype)
                fields.grid(row=i,column=2)

    def _getFileSpecificFields(self, parent, filename, filetype):
        archive = self.archive
        frame = Frame(parent)

        # show text button
        if filetype.startswith('Text'):
            view=Button(frame,
                        width=16, height=16,
                        image=doug.images.getImage("notepad.gif", self.app.root),
                        command=ArchivePanel._ViewFile((self.app,archive,filename)))
            view.pack()

        # show grid button
        if filetype.startswith("Vector"):
            view=Button(frame,
                        width=16, height=16,
                        image=doug.images.getImage("grid.gif", self.app.root),
                        command=ArchivePanel._ViewVector(self.app,archive,None,filename))
            view.pack()

        if filetype.startswith("Aggregates"):
            view=Button(frame,
                        width=16, height=16,
                        image=doug.images.getImage("grid.gif", self.app.root),
                        command=ArchivePanel._ViewAggregates(self.app,archive,None,filename))
            view.pack()
            
        if filetype.startswith("Matrix/Connections"):
            view=Button(frame,
                        width=16, height=16,
                        image=doug.images.getImage("grid.gif", self.app.root),
                        command=ArchivePanel._ViewConnections(self.app,archive,None,filename))
            view.pack()

        if filetype.startswith("Grid"):
            view=Button(frame,
                        width=16, height=16,
                        image=doug.images.getImage("grid.gif", self.app.root),
                        command=ArchivePanel._ViewGrid(self.app,archive,filename,None))
            view.pack()


        return frame

    class _ViewFile:
        def __init__(self, params):
            self.params=params

        def __call__(self):
            app, archive, fileName = self.params
            fullpath = os.path.join(archive.directoryName, fileName)
            viewer=app.config.get('global', 'viewer')
            os.spawnlp(os.P_NOWAIT, viewer, viewer, fullpath)

    class _ViewGrid:
        def __init__(self, app, archive, gridFileName, otherFileName):
            self.app = app
            self.archive = archive
            self.gridFileName = gridFileName
            self.otherFileName = otherFileName

        def getGridFile(self):
            gridFileName = self.gridFileName
            if gridFileName:
                gridFilePath = os.path.join(self.archive.directoryName, gridFileName)
            else:
                gridFiles = self.archive.getFiles(filetype='Grid')
                if gridFiles:
                    gridFilePath = os.path.join(self.archive.directoryName, gridFiles[0])
                else:
                    problemName = self.archive.info.get('general','problem-name')
                    problemArchive = self.app.getProblemArchive(problemName)
                    if problemArchive==None:
                        return            
                    gridFiles = problemArchive.getFiles(filetype='Grid')
                    gridFilePath = os.path.join(problemArchive.directoryName, gridFiles[0])

            return gridFilePath

        def __call__(self):
            gridFile = self.getGridFile()
            self._plot(gridFile)

        def _plot(self, gridFile, solutionFile=None, aggregatesFile=None, connectionsFile=None):
            config = DOUGConfigParser(name='Plot parameters')
            config.addConfig(self.app.config)
            configPanel = ConfigPanel(self.app.root, config, title="Plot configuration",
                                      sectionNames=['gridplot'])
            
            if not configPanel.done:
                return
            
            LOG.info("Show plot of %s", gridFile)

            gridplot=os.path.join(os.path.dirname(sys.argv[0]), 'gridplot.py')
            args = [gridplot, gridplot]
            args.extend(['--gin', gridFile])
            if solutionFile:
                args.extend(['--sin', solutionFile])
            if aggregatesFile:
                args.extend(['--ain', aggregatesFile])
            if connectionsFile:
                args.extend(['--cin', connectionsFile])

            plargs=[]
            device=config.get('gridplot', 'device')
            plargs.extend(['-dev', device])
            aggr=config.get('gridplot', 'aggr', None)
            if aggr:
                args.extend(['--aggr', aggr])
            args.extend(['--plplot', " ".join(plargs)])

            LOG.debug("Spawn '%s'", args)
            os.spawnvp(os.P_NOWAIT, 'python', args)

    class _ViewConnections(_ViewGrid):
        def __call__(self):
            gridFilePath = self.getGridFile()
            connectionsFilePath = os.path.join(self.archive.directoryName, self.otherFileName)
            self._plot(gridFilePath, connectionsFile=connectionsFilePath)

    class _ViewVector(_ViewGrid):
        def __call__(self):
            gridFilePath = self.getGridFile()
            vectorFilePath = os.path.join(self.archive.directoryName, self.otherFileName)
            self._plot(gridFilePath, solutionFile=vectorFilePath)
            
    class _ViewAggregates(_ViewGrid):
        def __call__(self):
            gridFilePath = self.getGridFile()
            aggrsFilePath = os.path.join(self.archive.directoryName, self.otherFileName)
            self._plot(gridFilePath, aggregatesFile=aggrsFilePath)

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

        self.updateFields()

    def _solve(self, initialConfig, problemArchive):
        try:
            config = DOUGConfigParser(name='DOUG execution parameters')
            config.addConfigContents(doug.execution.getDefaultConfigContents())
            config.addConfig(initialConfig)
            
            suffix = time.strftime('--run-%y%m%d-%H%M%S')
            config.set('DEFAULT', 'doug-workdir', problemArchive.directoryName+suffix)
            config.set('DEFAULT', 'doug-datadir', problemArchive.directoryName)

            matrixFiles = problemArchive.getFiles('Matrix')
            if matrixFiles:
                matrixFile='%(doug-datadir)s/'+matrixFiles[0]
            else:
                matrixFile=''
            config.set('doug-controls', 'assembled_mtx_file', matrixFile)
                           
            rhsFiles = problemArchive.getFiles('Vector/RHS')
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
                archive = Archive(problemArchive.name+suffix,
                                  directoryName=problemArchive.directoryName+suffix,
                                  archiveType='solution')
                archive.info.set('general','problem-name',problemArchive.name)
                archive.info.addConfig(execution.config)

                archive.info.set('doug-execution', 'date', time.strftime('%y-%m-%d %H:%M:%S'), addsection=True)
                
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
                if archive.info.has_option('doug-result', 'profilefile'):
                    archive.setFileType(archive.info.get('doug-result', 'profilefile'), 'Text/Profile')
                    
            finally:
                self.app.addArchive(archive)
        except(Exception), e:
            LOG.error(e, exc_info=e)


class ProblemArchivePanel(ArchivePanel):
    def __init__(self, *args, **kargs):
        ArchivePanel.__init__(self, *args, **kargs)

        buttonFrame = self.buttonFrame
        self.solveButton = Button(buttonFrame, text="Solve", command=self._solveCallback, fg='blue')
        self.solveButton.pack(side=RIGHT)        

    def _solveCallback(self):
        if self.archive==None:
            return
        
        self._solve(self.globalConfig, self.archive)
    

class SolutionArchivePanel(ArchivePanel):
    def __init__(self, *args, **kargs):
        ArchivePanel.__init__(self, *args, **kargs)

        buttonFrame = self.buttonFrame
        self.showGridButton = Button(buttonFrame, text="Show info", command=self._showInfo)
        self.showGridButton.pack(side=RIGHT)

        row,column = self.infoFrame.size()
        row = row+1

        label = Label(self.infoFrame, text="problem name")
        label.grid(row=row)
        self.problemNameField = Entry(self.infoFrame, state='readonly', width=40)
        self.problemNameField.grid(row=row, column=1, sticky=W)

        buttonFrame = self.buttonFrame
        self.solveButton = Button(buttonFrame, text="reSolve", command=self._solveCallback, fg='blue')
        self.solveButton.pack(side=RIGHT)        

    def _showInfo(self):
        configPanel=ConfigPanel(self.frame, self.archive.info,
                                readonly=True, sectionNames=['DEFAULT', 'doug', 'doug-controls', 'doug-result'],
                                title='Execution info')
                    
    def _solveCallback(self):
        if self.archive==None:
            LOG.warn('No archive selected')
            return

        problemName = self.archive.info.get('general', 'problem-name')        
        problemArchive = self.app.getProblemArchive(problemName)
        if problemArchive==None:
            LOG.warn('Problem %s not found', problemName)
            return            

        self._solve(self.archive.info, problemArchive)

    def setArchive(self, archive):
        ArchivePanel.setArchive(self, archive)

        if archive==None:
            self._setEntryField(self.problemNameField, '')
            return
        
        problemName = archive.info.get('general', 'problem-name')

        self._setEntryField(self.problemNameField, problemName)
