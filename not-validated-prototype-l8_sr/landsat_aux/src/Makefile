#-----------------------------------------------------------------------------
# Makefile for combine L8 auxiliary code
#-----------------------------------------------------------------------------

# Set up compile options
CC = gcc
RM = rm -f
EXTRA = -Wall -O2

# Define the include files
INC = combine_l8_aux_data.h
INCDIR = -I. -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_GCTPINC) -I$(XML2INC) \
-I$(ESPAINC)
NCFLAGS  = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC = get_args.c          \
      combine_l8_aux_data.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_common \
  -L$(HDFEOS_LIB) -lhdfeos -L$(HDFLIB) -lmfhdf -ldf -L$(JPEGLIB) -ljpeg \
  -L$(XML2LIB) -lxml2 -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define the executable
EXE	= combine_l8_aux_data

# Target for the executable
all: $(EXE)

combine_l8_aux_data: $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

install:
	install -d $(PREFIX)/bin
	install -m 755 $(EXE) $(PREFIX)/bin
	install -m 755 ../scripts/*.py $(PREFIX)/bin

clean:
	$(RM) *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

