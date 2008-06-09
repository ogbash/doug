from Tkinter import *
import Pmw
import tkFileDialog

import os

import logging
LOG = logging.getLogger(__name__)

import shutil

class ArchiveListbox(Pmw.ScrolledListBox):
    
    def __init__(self, parent, archiveList=None, archivePanel=None, *args, **kwargs):
        Pmw.ScrolledListBox.__init__(self, parent,
                                     selectioncommand=self.__listboxSelect,
                                     *args, **kwargs),
        self.archiveList = archiveList
        self.archivePanel = archivePanel

        self.component('listbox').bind('<Delete>', self.__listboxDelete)
        self.component('listbox').bind('<2>', self.__addFiles)


    def selectArchive(self, archive):
        self.setvalue(archive.name)
        
    def __listboxSelect(self, ev=None):
        widget=self
        widget.component('listbox').focus_set()
        index = int(widget.curselection()[0])
        archive = self.archiveList[index]
        self.archivePanel.setArchive(archive)

    def __listboxDelete(self, e):
        listbox = self
        archives = self.archiveList
        panel = self.archivePanel
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

    def __addFiles(self, ev=None):
        sel = self.curselection()
        if not sel:
            return

        archive = self.archiveList[int(sel[0])]

        filenames = tkFileDialog.askopenfilenames(title="Select files to include",
                                                      filetypes=[('All',"*"), ("Text data","*.txt")])
        if filenames:
            #isMove= tkMessageBox.askyesno(title="%d files selected. Copy or move?" % len(filenames),
            #                              message="Move files?",
            #                              default="no")
            for filename in filenames:
                newfilename = os.path.join(archive.directoryName, os.path.basename(filename))
                shutil.copy(filename, newfilename)

