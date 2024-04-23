# Directory with header files
IDIR := ./include
# Directory with source code
SRC := ./src
# Directory where .o files will be stored
ODIR := ./obj
# Directory where the executable will be built
BUILD := ./build


# Information about the compiler and the flags one wants to use
CC := gcc
CFLAGS := -I$(IDIR) -g -Wall -O3 -fno-math-errno -fno-trapping-math -fopenmp
EXECUTABLE := exe


# Dynamically linked libraries used (e.g. -lm, -lgsl, -lopenblas, etc...)
LIBS := -lm


# Header files
DEPS := $(find $(IDIR) -name '*.h' -type f)


# Files with .o extension (i.e. single .c files compiled) 
SRCS := $(shell find $(SRC) -name '*.c' | sed 's|^./src/||')
OBJ := $(SRCS:%=$(ODIR)/%.o)


$(ODIR)/%.c.o: $(SRC)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(BUILD)/$(EXECUTABLE) : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
