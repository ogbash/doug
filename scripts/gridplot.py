
from plplot import *
import sys
import numpy
import math
import getopt

class Plot:

    def __init__(self, gridFile, aggregateFile=None, solutionFile=None, contour=False):
        self.gridFile = gridFile
        self.aggregateFile = aggregateFile
        self.solutionFile = solutionFile
        self.contour = contour
        self.mappings = {} #: num of nodes -> array of map indices
        self.markings = {} #: name -> object specifying marks

    def readGrid(self, filename):
        print "Reading grid from %s" % filename
        f = open(filename, 'r')
        try:
            nvertices, nelements = map(int, f.readline().split())
            self.mappings[nvertices] = numpy.arange(nvertices) # default mapping to itself

            xs = numpy.zeros(nvertices)
            ys = numpy.zeros(nvertices)
            for i in xrange(nvertices):
                n,x,y = f.readline().split()
                n,x,y = int(n), float(x), float(y)
                xs[n] = x
                ys[n] = y

            coefs = numpy.zeros(nelements)
            elems = [None for j in xrange(nelements)]
            for j in xrange(nelements):
                elems[j] = map(int, f.readline().split())
                coefs[j] = float(f.readline())

            # try to read node mapping for dirichlet case
            try:
                num=int(f.readline())
            except (Exception), e:
                pass
            mapping = numpy.zeros(num, dtype=int)
            mapping[:] = -1
            for i in xrange(nvertices):
                n = int(f.readline())
                if n!=-1:
                    mapping[n] = i
            print "Putting mapping for %d" % num
            self.mappings[num] = mapping

        finally:
            f.close()

        return xs,ys,elems,coefs

    def readAggregates(self, filename, nvertices):
        print "Reading aggregates from %s" % filename
        f = open(filename)
        try:
            nmarks, nnodes = map(int, f.readline().split())
            aggrs = numpy.zeros(nvertices, dtype=int)
            aggrs[:] = -1
            mapping = self.mappings[nnodes]
            for i in xrange(nnodes):
                n = int(f.readline())
                aggrs[mapping[i]] = n
            self.markings['aggregates'] = aggrs
        finally:
            f.close()

    def readSolution(self, filename, nvertices):
        print "Reading solution from %s" % filename
        f = open(filename)
        try:
            nnodes = int(f.readline())

            solution = numpy.zeros(nvertices)
            mapping = self.mappings[nnodes]
            for i in xrange(nnodes):
                v = float(f.readline())
                solution[mapping[i]] = v
            self.markings['solution'] = solution
        finally:
            f.close()

    def drawGrid(self, x,y,elems,coefs):
        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y) 
        plenv(xmin, xmax, ymin, ymax, 0, 0)

        maxcoef = max(coefs)
        mincoef = min(coefs)

        aggrs = self.markings.get('aggregates', None)
        solution = self.markings.get('solution', None)
        if solution!=None:
            maxsol = max(solution)
            minsol = min(solution)
            pllab("x", "y", "Solution: %s to %s" % (minsol, maxsol))
        else:
            pllab("x", "y", "Coefficient: %s to %s" % (mincoef, maxcoef))

        i = [0.,1]
        r = numpy.array((0.,.5))
        g = numpy.array((.5,0.))
        b = numpy.array((0.,0.))

        #plscmap1n(32)
        plscmap1l(1, i, r, g, b, [0,0])

        #lx = numpy.zeros(4)
        #ly = numpy.zeros(4)
        for i in xrange(len(elems)):
            elem = elems[i]
            lx=x[elem]
            ly=y[elem]

            if solution!=None:
                average = sum(solution[elem])/len(elem)
                if maxsol-minsol==0:
                    plcol1(0.0)
                else:
                    plcol1((average-minsol)/(maxsol-minsol))
            else:
                if maxcoef-mincoef == 0:
                    plcol1(0.0)
                else:
                    plcol1((coefs[i]-mincoef)/(maxcoef-mincoef))
            plpsty(0)
            plfill(lx,ly)

            n = len(lx)
            for j in xrange(n):
                if aggrs!=None and aggrs[elem[j]]==aggrs[elem[(j+1)%n]]:
                    plcol0(aggrs[elem[j]]%14+1) # except 0 (black) and 15 (white)
                elif self.contour:
                    plwid(3)
                    plcol0(1)
                else:
                    plwid(1)
                    plcol0(0)
                pljoin(lx[j],ly[j],lx[(j+1)%n],ly[(j+1)%n])

            # fill
            if aggrs!=None:
                m = aggrs[elem[0]]
                anyoutside = (aggrs[elem]==-1).any()
                if not anyoutside and \
                       (aggrs[elem]==m).all(): # all nodes in the same aggr
                    plpsty(4)
                    plfill(lx,ly)

        # draw aggregate nodes
        if solution!=None:
            for i in xrange(len(x)):
                plcol1((solution[i]-minsol)/(maxsol-minsol+0.001))
                plssym(0,.2)
                plpoin(x[i:i+1],y[i:i+1],4)

    ## main program

    def run(self, dev=None):
        if dev:
            plsdev(dev)
        plinit()

        xs,ys,elems,coefs = self.readGrid(self.gridFile)
        if self.aggregateFile:
            self.readAggregates(self.aggregateFile, len(xs))
        if self.solutionFile:
            self.readSolution(self.solutionFile, len(xs))

        self.drawGrid(xs,ys,elems,coefs)

        plend()

def main():
    options = ""
    loptions = [
        'fortran',
        'plplot=',
        'gin=',
        'ain=',
        'sin='
        ]

    plargs = [sys.argv[0]]
    GRIDIN = 'grid.txt'
    AGGRSIN = None
    SOLUTIONIN = None

    opts = getopt.gnu_getopt(sys.argv, options, loptions)
    for key, value in opts[0]:
        if key=='--plplot':
            plargs.extend(value.split())
        elif key=='--gin':
            GRIDIN=value
        elif key=='--ain':
            AGGRSIN=value
        elif key=='--sin':
            SOLUTIONIN=value

    plparseopts(plargs, PL_PARSE_FULL)

    plot = Plot(GRIDIN, AGGRSIN, SOLUTIONIN)
    plot.run()
    
if __name__=="__main__":
    main()
