IDIR =../libStatGen/include
IDIR2 =../htslib-1.13
#CC = g++ -std=c++17 -fPIC  -DHTSLIB -Wall -O4 
CC = g++ -std=c++17 -fPIC  -DHTSLIB -DGZSTREAM -Wall -O4 
#CC = g++ -std=c++17 -fPIC  -DLIBSTATGEN
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
_TARGET2_OBJ = tagsort.o htslib_tagsort.o gzstream.o globals.o sort_write.o
TARGET2_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET2_OBJ))


install: $(TARGET1) $(TARGET2)
	cp $(TARGET1) ../bin/
	cp $(TARGET2) ../bin/
	cp ../htslib-1.13/*.so.? ../bin/

all: $(TARGET1) $(TARGET2)

_COMMON_OBJ =  utilities.o inputoptions.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_COMMON_OBJ))

$(ODIR)/%.o: %.cpp $(_DEPS)
	$(CC) -c -o $@ $<  -I$(IDIR) -I. -I$(IDIR2)

$(TARGET1): $(OBJ) $(TARGET1_OBJ) 
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

$(TARGET2): $(OBJ) $(TARGET2_OBJ)  
	$(CC) -Wl,-rpath,/usr/local/bin:fastqpreprocessing/bin:bin:. -o $@ $(OBJ)  $(TARGET2_OBJ)  -I. -L. -lstdc++fs -lz  -L../htslib-1.13 -lhts -lpthread

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  *.o *.so *.a
	rm -rf $(TARGET1) $(TARGET2)
