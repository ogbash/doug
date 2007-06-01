
def diffVectors(v1, v2):
    """Takes two vectors of values and returns difference"""

    if len(v1) != len(v2):
            raise self.failureException("Sizes of vectors are different: %d, %d"
                                        % (len(v1), len(v2)) )
    dif = 0.
    for i in xrange(0,len(v1)):
            dif += pow(v1[i]-v2[i], 2)
    dif = pow(dif, 0.5)

    return dif
