
from Tkinter import *
import Pmw

class ResultTable:

    def __init__(self, app, master):
        self.frame = Pmw.ScrolledFrame(master)

        l=Label(self.frame.interior(), text='name', fg='blue')
        l.grid(row=0)
        l=Label(self.frame.interior(), text='res', fg='blue')
        l.grid(row=0, column=1)
        l=Label(self.frame.interior(), text='lev', fg='blue')
        l.grid(row=0, column=2)
        l=Label(self.frame.interior(), text='sol', fg='blue')
        l.grid(row=0, column=3)
        l=Label(self.frame.interior(), text='rad1', fg='blue')
        l.grid(row=0, column=4)
        l=Label(self.frame.interior(), text='rad2', fg='blue')
        l.grid(row=0, column=5)
        l=Label(self.frame.interior(), text='iter', fg='blue')
        l.grid(row=0, column=6)
        l=Label(self.frame.interior(), text='prec. t', fg='blue')
        l.grid(row=0, column=7)

        for i,archive in enumerate(app.solutionsArchives):
            l=Label(self.frame.interior(), text=archive.name)
            l.grid(row=i+1)
            l=Label(self.frame.interior(), text=archive.info.get('doug-result', 'returnvalue'))
            l.grid(row=i+1, column=1)
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'levels'))
            l.grid(row=i+1, column=2)
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'solver'))
            l.grid(row=i+1, column=3)
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'radius1'))
            l.grid(row=i+1, column=4)
            l=Label(self.frame.interior(), text=archive.info.get('doug-controls', 'radius2'))
            l.grid(row=i+1, column=5)
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'pcg-iterations', default=''))
            l.grid(row=i+1, column=6)
            l=Label(self.frame.interior(), text=archive.info.get('doug-profile', 'preconditioner-time', default=''))
            l.grid(row=i+1, column=7)
