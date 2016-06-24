# Declaration of variables
CC = g++

# Directories
OBJDIR = obj
SRCDIR = src
INCDIR = include
LIBDIR = lib

# File names
EXEC = nmr_analysis 
SHARED_LIBRARY = libnmranalysis.so
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
SHARED_OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o, $(filter-out src/main.cpp, $(SOURCES)))


ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs) -lThread -lGui -lMinuit

CC_FLAGS = -I $(INCDIR) -w $(ROOT_CFLAGS)

all: $(EXEC) $(SHARED_LIBRARY)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(ROOT_LIBS) -o $(EXEC) -std=c++11 -Wall

#Shared library
$(SHARED_LIBRARY): $(SHARED_OBJECTS)
	$(CC) $(SHARED_OBJECTS) $(ROOT_LIBS) -o $(LIBDIR)/$(SHARED_LIBRARY) -std=c++11 -Wall -shared

# To obtain object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@  -std=c++11 -Wall

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJDIR)/*.o $(LIBDIR)/$(SHARED_LIBRARY)
