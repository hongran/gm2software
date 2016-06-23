# Declaration of variables
CC = g++

# Directories
OBJDIR = obj
SRCDIR = src
INCDIR = include

# File names
EXEC = nmr_analysis
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs) -lThread -lGui -lMinuit

CC_FLAGS = -I $(INCDIR) -w $(ROOT_CFLAGS)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(ROOT_LIBS) -o $(EXEC) -std=c++11 -Wall

# To obtain object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@  -std=c++11 -Wall

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJDIR)/*.o
