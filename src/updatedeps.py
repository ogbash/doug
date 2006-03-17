#!/usr/bin/env python
import sys, getopt, os, fortdeps

def usage():
    print 'usage: -d -s name'
    
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ds:h", ["help", "all"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    depskey = 0
    for o, a in opts:
        if o == "-d":
            depskey = 1
            d = fortdeps.depgen('.')
            d.makedeps('Make.deps')
        if o == "-s":
            if depskey == 0:
                d = fortdeps.depgen('.')
                d.makedeps('Make.deps')
            sf = open('Make.srcs', 'w')
            sf.write('OBJ_'+ a + ' =' + d.getsrcs(a) + '\n')
        if o in ("-h", "--help"):
            usage()
            sys.exit()
            if o == "--all":
                d = fortdeps.depgen('.', be_recursive=True, verbose=2)
                d.write_depends(".", "Make.deps")
                sys.exit()
                
if __name__ == "__main__":
    main()
    
