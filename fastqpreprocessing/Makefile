IDIR1 = libStatGen/include
IDIR2 = htslib-1.13
IDIR3 = gzstream

CC = g++ -std=c++17 -fPIC  -DHTSLIB -Wall -O4  -Wwrite-strings

CFLAGS = -I$(IDIR1) -LlibStatGen  -Lgzstream

LIBS = -LlibStatGen -lStatGen -lz -lpthread -lstdc++fs -Lgzstream -lgzstream

_DEPS = src/utilities.h src/input_options.h src/fastq_common.h

TARGET1 = bin/fastqprocess
TARGET1_OBJ = obj/fastqprocess.o

TARGET2 = bin/TagSort
TARGET2_OBJ = obj/tagsort.o obj/htslib_tagsort.o obj/metricgatherer.o

TARGET3 = bin/fastq_slideseq
TARGET3_OBJ = obj/fastq_slideseq.o

TARGET4 = bin/fastq_metrics
TARGET4_OBJ = obj/fastq_metrics.o

install: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
	cp htslib-1.13/*.so.? bin/

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

COMMON_OBJ = obj/utilities.o obj/input_options.o obj/fastq_common.o

obj/%.o: src/%.cpp $(_DEPS)
	$(CC) -c -o $@ $<  -I$(IDIR1) -I$(IDIR2) -I$(IDIR3)

$(TARGET1): $(COMMON_OBJ) $(TARGET1_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET2): $(COMMON_OBJ) $(TARGET2_OBJ)
	$(CC) -Wl,-rpath,/usr/local/bin:fastqpreprocessing/bin:bin:. -o $@ $(COMMON_OBJ)  $(TARGET2_OBJ) $(LIBS) -Lhtslib-1.13 -lhts

$(TARGET3): $(COMMON_OBJ) $(TARGET3_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET4): $(COMMON_OBJ) $(TARGET4_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

.PHONY: clean
clean:
	rm -f obj/*.o *~ core $(INCDIR)/*~  *.o *.so *.a
	rm -rf $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
