CC=g++

SRC_DIR = ./src
INC_DIR = ./include
OBJ_DIR = ./obj
EXE_DIR = ./bin

ROOT_FLAGS = `root-config --cflags`
ROOT_LIBS = `root-config --libs`

BARCODE_DEPS = $(INC_DIR)/gm2_TBarcode.h

_BARCODE_OBJ = gm2_TBarcode.o gm2_test_Barcode_main.o
BARCODE_OBJ = $(patsubst %,$(OBJ_DIR)/%,$(_BARCODE_OBJ))

_EXE = gm2_test_Barcode
EXE = $(patsubst %,$(EXE_DIR)/%,$(_EXE))
CFLAGS = -I$(INC_DIR) 

all : $(EXE)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx $(BARCODE_DEPS)
	$(CC) -c -o $@ $< $(ROOT_FLAGS) $(CFLAGS) -std=c++11 -Wall

$(EXE): $(BARCODE_OBJ)
	$(CC) -o $@ $^ $(ROOT_LIBS) -std=c++11 -Wall

.PHONY: clean
clean:
	rm $(OBJ_DIR)/*.o $(EXE)
