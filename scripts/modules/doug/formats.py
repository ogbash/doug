
class Fortran:
    def readVectorFromFile(self, f):
        from array import array
        vec = array('d')
        # we need 4 byte integer, ugly hack for 64bit platforms
        intarr = array('l')
        if intarr.itemsize != 4:
            intarr = array('i')
        
        intarr.fromfile(f, 1) # read start marker (fortran)
        smarker = intarr[-1]
        vec.fromfile(f, smarker/vec.itemsize) # read values
        intarr.fromfile(f, 1) # read end marker
        emarker = intarr[-1]

        if smarker != emarker:
            raise Exception("Invalid fortran file, start (%d) and end markers (%d) are not equal." %
                            (smarker, emarker))

        return vec
