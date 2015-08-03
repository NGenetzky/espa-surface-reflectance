#------------------------------------------------------------------------------
# Makefile
#
# Project Name: surface reflectance
# Makefile that will invoke subsequent Makefiles in subdirectories.
# Assumes the following Makefiles exist:
#	ledapsSrc/Makefile
#------------------------------------------------------------------------------

all:
	echo "make all in ledapsSrc"; \
        (cd ledapsSrc; $(MAKE) all -f Makefile);

install:
	echo "make install in ledapsSrc"; \
        (cd ledapsSrc; $(MAKE) install -f Makefile);

clean:
	echo "make clean in ledapsSrc"; \
        (cd ledapsSrc; $(MAKE) clean -f Makefile);

