EXTRA   = -D_BSD_SOURCE -Wall -O2

INCDIR  = -I. -I$(TIFFINC) -I$(GEOTIFF_INC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC)
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

EXLIB	= -L$(GEOTIFF_LIB) -lgeotiff -L$(TIFFLIB) -ltiff \
          -L$(HDFEOS_LIB) -lhdfeos -L$(HDFEOS_GCTPLIB) -lGctp \
          -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1 = lndcsm
OBJ1    = lndcsm.o degdms.o param.o input.o lut.o output.o csm.o space.o \
          names.o  myhdf.o mystring.o error.o tiff.o virbuf.o date.o util.o
INC1    = lndcsm.h keyvalue.h param.h input.h lut.h output.h csm.h names.h \
          date.h myhdf.h mystring.h bool.h const.h error.h tiff.h virbuf.h \
          util.h space.h myproj.h myproj_const.h


all: $(TARGET1)

x: $(TARGET1)

$(OBJ1): $(INC1)

$(TARGET1): $(OBJ1)
	$(CC) $(EXTRA) -o $(TARGET1) $(OBJ1) $(LOADLIB)

clean:
	rm -f *.o $(TARGET1)

install:
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET1) $(PREFIX)/bin

#
# Rules
#
.c.o:
	$(CC) $(EXTRA) $(NCFLAGS) -c $< -o $@
