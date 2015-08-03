#------------------------------------------------------------------------------
# Makefile
#
# Project Name: surface reflectance
# Makefile that will invoke subsequent Makefiles in subdirectories.
# Assumes the following Makefiles exist:
#	bin/Makefile
#	src/Makefile
#------------------------------------------------------------------------------

all:
	echo "make all in bin"; \
        (cd bin; $(MAKE) all -f Makefile);

	echo "make all in src"; \
        (cd src; $(MAKE) all -f Makefile);

install:
	echo "make install in bin"; \
        (cd bin; $(MAKE) install -f Makefile);
	echo "make install in src"; \
        (cd src; $(MAKE) install -f Makefile);

clean:
	echo "make clean in bin"; \
        (cd bin; $(MAKE) clean -f Makefile);
	echo "make clean in src"; \
        (cd src; $(MAKE) clean -f Makefile);

