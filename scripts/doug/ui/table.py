
from Tkinter import *
import Pmw

class _SelectSolution:
    def __init__(self, app, archive):
        self.app = app
        self.archive = archive

    def __call__(self, ev):
        self.app.solutionsListbox.selectArchive(self.archive)
        self.app.notebook.selectpage('solutions')

class ResultTable:

    def __init__(self, app, master):
        self.app = app
        self.master = master
        self.frame = Pmw.ScrolledFrame(master)
        self.__createTable(master)

    def refresh(self):
        self.__createTable(self.master)
        
    def __createTable(self, master):
        for slave in self.frame.interior().grid_slaves():
            slave.grid_remove()
            slave.destroy()

        col=0
        l=Label(self.frame.interior(), text='name', fg='blue')
        l.grid(row=0) ; col=col+1
        l=Label(self.frame.interior(), text='res', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='lev', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='sol', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='cmeth', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        
        l=Label(self.frame.interior(), text='rad1', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='fine aggr', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='rad2', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='coarse aggr', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='iter', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='time', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='prec. time', fg='blue')
        l.grid(row=0, column=col) ; col=col+1
        l=Label(self.frame.interior(), text='date', fg='blue')
        l.grid(row=0, column=col) ; col=col+1

        for i,archive in enumerate(self.app.solutionsArchives):
            col=0
            l=Label(self.frame.interior(), text=archive.name)
            l.bind("<Double-Button-1>", _SelectSolution(self.app, archive))
            l.grid(row=i+1) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-result', 'returnvalue', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'levels'))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'solver'))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'coarse_method'))
            l.grid(row=i+1, column=col) ; col=col+1
            
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'radius1'))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'fine-aggregates', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'radius2'))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'coarse-aggregates', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'pcg-iterations', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'iterations-time', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'preconditioner-time', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
            l=Label(self.frame.interior(), text=archive.info.get('doug-execution', 'date', default=''))
            l.grid(row=i+1, column=col) ; col=col+1
