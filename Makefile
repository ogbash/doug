# Makefile

include src/Make.def

all: lib main aggr

deps:
	cd src/utils      && $(MAKE) deps
	cd src/datatypes  && $(MAKE) deps
	cd src/solvers    && $(MAKE) deps
	cd src/main       && $(MAKE) deps

lib: deps 
	cd lib/            ; rm -f *.so*
	cd src/utils      && $(MAKE) globalLib
	cd src/datatypes  && $(MAKE) globalLib
	cd src/solvers    && $(MAKE) globalLib
	$(DLL) -o lib/lib$(DOUG_LIB).so.$(VMAJOR) $(DOUG_OBJS) && \
	cd lib/ && ln -s lib$(DOUG_LIB).so.$(VMAJOR) lib$(DOUG_LIB).so

main:
	cd src/main       && $(MAKE) main

aggr:
	cd src/main       && $(MAKE) aggr

libs:
	cd src/utils      && $(MAKE) lib
	cd src/datatypes  && $(MAKE) lib
	cd src/solvers    && $(MAKE) lib

drivers:
	cd src/datatypes  && $(MAKE) all_drivers
	cd src/solvers    && $(MAKE) all_drivers

clean:
	cd src/main        && $(MAKE) clean
	cd src/solvers     && $(MAKE) clean_all #cleans in ./drivers as well
	cd src/datatypes   && $(MAKE) clean_all #cleans in ./drivers as well
	cd src/utils       && $(MAKE) clean
	cd src/; rm -f *~
	rm -f *~

clean_lib:
	cd lib/           ; rm -rf *.a *.so*
