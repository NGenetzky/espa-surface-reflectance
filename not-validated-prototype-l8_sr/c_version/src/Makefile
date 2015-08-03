#------------------------------------------------------------------------------
# Makefile for L8 SR code
#------------------------------------------------------------------------------

# Set up compile options; add -fopenmp to use the OpenMP support for processing
# L8SR.
# EXTRA = -Wall -O2 -fopenmp
CC = gcc
RM = rm -f
EXTRA = -Wall -O2

# Define the include files
INC = common.h date.h input.h output.h lut_subr.h l8_sr.h
INCDIR = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC) -I$(XML2INC) \
-I$(ESPAINC)
NCFLAGS  = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC = compute_refl.c      \
      date.c              \
      get_args.c          \
      input.c             \
      lut_subr.c          \
      output.c            \
      subaeroret.c        \
      l8_sr.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -l_espa_format_conversion \
  -L$(HDFEOS_LIB) -lhdfeos -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg \
  -L$(XML2LIB) -lxml2 -L$(HDFEOS_GCTPLIB) -lGctp -L$(ZLIBLIB) -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define the executable
EXE	= l8_sr

# Target for the executable
all: $(EXE)

l8_sr: $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

install:$(EXE)
	install -d $(PREFIX)/bin
	install -m 755 $(EXE) $(PREFIX)/bin

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

