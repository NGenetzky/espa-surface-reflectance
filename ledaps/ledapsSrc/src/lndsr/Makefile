EXTRA   = -g -D_BSD_SOURCE -Wall -O2

INCDIR  = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC) -I$(XML2INC) -I$(ESPAINC)
NCFLAGS  = $(CFLAGS) $(EXTRA) $(INCDIR)

EXLIB	= -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -l_espa_format_conversion \
  -L$(HDFEOS_LIB) -lhdfeos -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg \
  -L$(XML2LIB) -lxml2 -L$(HDFEOS_GCTPLIB) -lGctp -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

TARGET1	= lndsr
OBJ1    = lndsr.o param.o input.o prwv_input.o lut.o output.o sr.o ar.o \
          date.o mystring.o error.o grib.o read_grib_tools.o myhdf.o \
          CHAND.o CSALBR.o sixs_runs.o clouds.o
INC1    = lndsr.h keyvalue.h param.h input.h prwv_input.h lut.h output.h \
          sr.h ar.h date.h mystring.h bool.h const.h error.h grib.h myhdf.h \
          read_grib_tools.h myproj.h myproj_const.h

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
.f.o:
	gfortran -c $< -o $@
