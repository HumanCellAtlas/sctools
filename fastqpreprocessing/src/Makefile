IDIR =../libStatGen/include

CC = g++ -std=c++17 -O4
CFLAGS = -I$(IDIR) -L../libStatGen

ODIR=obj
LDIR =../libStatGen/

LIBS = -lStatGen -lz -lpthread -lstdc++fs

_DEPS = fastqprocess.h utilities.h inputoptions.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

TARGET1 = fastqprocess
_TARGET1_OBJ = fastqprocess.o 
TARGET1_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET1_OBJ))

TARGET2 = TagSort
_TARGET2_OBJ = tagsort.o 
TARGET2_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET2_OBJ))

all: $(TARGET1) $(TARGET2)


_COMMON_OBJ =  utilities.o inputoptions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_COMMON_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET1): $(OBJ) $(TARGET1_OBJ)  $(_DEPS)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

$(TARGET2): $(OBJ) $(TARGET2_OBJ)  $(_DEPS)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -rf $(TARGET1) $(TARGET2)
