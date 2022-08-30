IDIR =../libStatGen/include
IDIR2 =../htslib-1.13
IDIR3 =../gzstream

CC = g++ -std=c++17 -fPIC  -DHTSLIB -Wall -O4  -Wwrite-strings

CFLAGS = -I$(IDIR) -L../libStatGen  -L../gzstream

ODIR=obj
LDIR =../libStatGen/

LIBS = -lStatGen -lz -lpthread -lstdc++fs -lgzstream

_DEPS = fastqprocess.h utilities.h input_options.h fastq_common.h

TARGET1 = fastqprocess
_TARGET1_OBJ = fastqprocess.o fastq_common.o
TARGET1_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET1_OBJ))

TARGET2 = TagSort
_TARGET2_OBJ = tagsort.o htslib_tagsort.o metricgatherer.o
TARGET2_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET2_OBJ))

TARGET3 = fastq_slideseq
_TARGET3_OBJ = fastq_slideseq.o fastq_common.o
TARGET3_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET3_OBJ))

TARGET4 = fastq_metrics
_TARGET4_OBJ = fastq_metrics.o
TARGET4_OBJ = $(patsubst %,$(ODIR)/%,$(_TARGET4_OBJ))

install: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
	cp $(TARGET1) ../bin/
	cp $(TARGET2) ../bin/
	cp $(TARGET3)  ../bin/
	cp $(TARGET4) ../bin/
	cp ../htslib-1.13/*.so.? ../bin/

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

_COMMON_OBJ =  utilities.o input_options.o
OBJ = $(patsubst %,$(ODIR)/%,$(_COMMON_OBJ))

$(ODIR)/%.o: %.cpp $(_DEPS)
	$(CC) -c -o $@ $<  -I$(IDIR) -I. -I$(IDIR2) -I$(IDIR3)

$(TARGET1): $(OBJ) $(TARGET1_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET2): $(OBJ) $(TARGET2_OBJ)
	$(CC) -Wl,-rpath,/usr/local/bin:fastqpreprocessing/bin:bin:. -o $@ $(OBJ)  $(TARGET2_OBJ)  -I. -L. -lstdc++fs -lz  -L../htslib-1.13 -lhts -lpthread

$(TARGET3): $(OBJ) $(TARGET3_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET4): $(OBJ) $(TARGET4_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  *.o *.so *.a
	rm -rf $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
