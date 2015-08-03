#------------------------------------------------------------------------------
# Makefile
#
# Project Name: surface reflectance
# Makefile that will invoke subsequent Makefiles in subdirectories.
# Assumes the following Makefiles exist:
#	c_version/Makefile
#------------------------------------------------------------------------------

all:
	echo "make all in c_version"; \
        (cd c_version; $(MAKE) all -f Makefile);

install:
	echo "make install in c_version"; \
        (cd c_version; $(MAKE) install -f Makefile);

clean:
	echo "make clean in c_version"; \
        (cd c_version; $(MAKE) clean -f Makefile);

