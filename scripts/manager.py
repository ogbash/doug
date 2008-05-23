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

import os
import os.path
import shutil

from doug.config import DOUGConfigParser
from doug.configui import ConfigPanel
from doug.archive import Archive, ProblemArchivePanel, SolutionArchivePanel

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

[global]
viewer: gedit

#: output for plplot
# @type: list
# @type-params: xwin, ps, psc
plplot-device: xwin

"""

class App:
    def __reload(self, ev):
        import doug.execution
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
        mainFrame=Pmw.NoteBook(root)
        self.notebook=mainFrame
        mainFrame.pack(expand=True, fill=BOTH)
        
        problemsFrame = LabelFrame(mainFrame.add('problems'), text='Problems', fg='blue')
        problemsFrame.pack(side=TOP, expand=True, fill=BOTH)
        solutionsFrame = LabelFrame(mainFrame.add('solutions'), text='Solutions', fg='blue')
        solutionsFrame.pack(side=TOP, expand=True, fill=BOTH)

        # problems listbox
        listbox = Pmw.ScrolledListBox(problemsFrame,
                                      selectioncommand=self.__problemsListboxSelect)
        self.problemsListbox = listbox
        listbox.pack(side=LEFT, fill=Y)
        listbox.component('listbox').bind('<Delete>', self.__problemsListboxDelete)
        listbox.component('listbox').bind('<2>', self.__addFiles)
        #listbox.event_add('<<AddFiles>>', '<2>')


        addb = Button(problemsFrame, text="Add archive", command=self.__addArchive)
        addb.pack(side=TOP)

        self.problemPanel = ProblemArchivePanel(problemsFrame, self.config, self)
        self.problemPanel.frame.pack(side=TOP, expand=True, fill=BOTH)

        # solutions listbox
        listbox = Pmw.ScrolledListBox(solutionsFrame,
                                      selectioncommand=self.__solutionsListboxSelect)
        self.solutionsListbox = listbox
        listbox.pack(side=LEFT, fill=Y)
        #listbox.bind('<<ListboxSelect>>', self.__solutionsListboxSelect)        
        listbox.component('listbox').bind('<Delete>', self.__solutionsListboxDelete)

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

    def __problemsListboxDelete(self, e):
        listbox=self.problemsListbox
        archives=self.problemsArchives
        self.__listboxDelete(listbox, archives, self.problemPanel)

    def __solutionsListboxDelete(self, e):
        listbox=self.solutionsListbox
        archives=self.solutionsArchives
        self.__listboxDelete(listbox, archives, self.solutionPanel)

    def __listboxDelete(self, listbox, archives, panel):
        names=listbox.getvalue()
        if not names:
            return
        for name in names:
            indices=filter(lambda v: v[1].name==name,
                           enumerate(archives))
            indices.reverse()
            for index, name in indices:
                archive=archives[index]
                del archives[index]
                listbox.delete(index)
                LOG.info('Deleting %s', archive.directoryName)
                shutil.rmtree("%s" % (archive.directoryName))
        panel.setArchive(None)

    def __editOptions(self):
        d = ConfigPanel(self.root, self.config, title='Global options')
        if d.done:
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
        widget.component('listbox').focus_set()
        index = int(widget.curselection()[0])
        archive = self.problemsArchives[index]
        self.problemPanel.setArchive(archive)

    def __solutionsListboxSelect(self, ev=None):
        widget=self.solutionsListbox
        widget.component('listbox').focus_set()
        index = int(widget.curselection()[0])
        archive = self.solutionsArchives[index]
        self.solutionPanel.setArchive(archive)

    def addArchive(self, archive):
        archiveType = archive.info.get('general', 'archive-type')
        if archiveType=='problem':
            self.problemsArchives.append(archive)
            self.problemsListbox.insert(END, archive.name)
            self.notebook.selectpage('problems')
            self.problemsListbox.focus()
        elif archiveType=='solution':
            self.solutionsArchives.append(archive)
            self.notebook.selectpage('solutions')
            self.solutionsListbox.focus()
            self.solutionsListbox.insert(END, archive.name)
            self.solutionsListbox.setvalue(archive.name)
            
    def __addArchive(self):
        name = tkSimpleDialog.askstring("Insert name", "Type new archive name (without .tar.gz)")
        if name:
            archive = Archive(name)
            self.addArchive(archive)

    def __addFiles(self, ev=None):
        sel = self.problemsListbox.curselection()
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
        self.config = DOUGConfigParser()
        self.config.addConfigContents(_configContent)
        if os.path.exists('.manager.conf'):
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
