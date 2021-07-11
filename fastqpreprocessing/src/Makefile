IDIR =../libStatGen/include
IDIR2 =../htslib-1.13
CC = g++ -std=c++17 -fPIC  -DHTSLIB -O4
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
_TARGET2_OBJ = tagsort.o 
TARGET2_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET2_OBJ))

all: $(TARGET1) $(TARGET2)

_COMMON_OBJ =  utilities.o inputoptions.o  libstatgen_tagsort.o sort_write.o htslib_tagsort.o
OBJ = $(patsubst %,$(ODIR)/%,$(_COMMON_OBJ))

$(ODIR)/%.o: %.cpp $(_DEPS)
	$(CC) -c -o $@ $<  -I$(IDIR) -I. -I$(IDIR2)

$(TARGET1): $(OBJ) $(TARGET1_OBJ) 
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

liblibstatgen_tagsort.a: libstatgen_tagsort.cpp
	$(CC)  -c  $< -I../libStatGen/include  -L../libStatGen -lStatGen -lz
	ar ru $@ libstatgen_tagsort.o
	ranlib  liblibstatgen_tagsort.a

liblibstatgen_tagsort.so: libstatgen_tagsort.cpp
	$(CC)  -shared -o $@ $< -I../libStatGen/include  -L../libStatGen -lStatGen -lz

libhtslib_tagsort.so: htslib_tagsort.cpp
	$(CC)  -shared -o $@ $< -I../htslib-1.13/ -L../htslib-1.13 -lhts

$(TARGET2): $(OBJ) $(TARGET2_OBJ) #libhtslib_tagsort.so #liblibstatgen_tagsort.so 
	$(CC) -o $@ $(OBJ)  $(TARGET2_OBJ)  -I. -L. -lstdc++fs  -I../libStatGen/include  -L../libStatGen -lStatGen -lz  -L../htslib-1.13 -lhts

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  *.o *.so *.a
	rm -rf $(TARGET1) $(TARGET2)
