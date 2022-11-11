IDIR1 = libStatGen/include
IDIR2 = htslib-1.13
IDIR3 = gzstream

CC = g++ -std=c++17 -fPIC  -DHTSLIB -Wall -O4  -Wwrite-strings

CFLAGS = -I$(IDIR1) -LlibStatGen  -Lgzstream

LIBS = -LlibStatGen -lStatGen -lz -lpthread -lstdc++fs -Lgzstream -lgzstream

_DEPS = src/utilities.h src/input_options.h src/fastq_common.h src/mitochondrial_gene_selector.h

TARGET_FQPROCESS = bin/fastqprocess
TARGET_FQPROCESS_OBJ = obj/fastqprocess.o

TARGET_TAGSORT = bin/TagSort
TARGET_TAGSORT_OBJ = obj/tagsort.o obj/htslib_tagsort.o obj/metricgatherer.o obj/mitochondrial_gene_selector.o

TARGET_FQSLIDESEQ = bin/fastq_slideseq
TARGET_FQSLIDESEQ_OBJ = obj/fastq_slideseq.o

TARGET_FQMETRICS = bin/fastq_metrics
TARGET_FQMETRICS_OBJ = obj/fastq_metrics.o

TARGET_SAMPLEFQ = bin/samplefastq
TARGET_SAMPLEFQ_OBJ = obj/samplefastq.o

install: $(TARGET_FQPROCESS) $(TARGET_TAGSORT) $(TARGET_FQSLIDESEQ) $(TARGET_FQMETRICS) $(TARGET_SAMPLEFQ)
	cp htslib-1.13/*.so.? bin/

all: $(TARGET_FQPROCESS) $(TARGET_TAGSORT) $(TARGET_FQSLIDESEQ) $(TARGET_FQMETRICS) $(TARGET_SAMPLEFQ)

COMMON_OBJ = obj/utilities.o obj/input_options.o obj/fastq_common.o

obj/%.o: src/%.cpp $(_DEPS)
	$(CC) -c -o $@ $<  -I$(IDIR1) -I$(IDIR2) -I$(IDIR3)

$(TARGET_FQPROCESS): $(COMMON_OBJ) $(TARGET_FQPROCESS_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET_TAGSORT): $(COMMON_OBJ) $(TARGET_TAGSORT_OBJ)
	$(CC) -Wl,-rpath,/usr/local/bin:fastqpreprocessing/bin:bin:. -o $@ $(COMMON_OBJ)  $(TARGET_TAGSORT_OBJ) $(LIBS) -Lhtslib-1.13 -lhts

$(TARGET_FQSLIDESEQ): $(COMMON_OBJ) $(TARGET_FQSLIDESEQ_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET_FQMETRICS): $(COMMON_OBJ) $(TARGET_FQMETRICS_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

$(TARGET_SAMPLEFQ): $(COMMON_OBJ) $(TARGET_SAMPLEFQ_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS)

.PHONY: clean
clean:
	rm -f obj/*.o *~ core $(INCDIR)/*~  *.o *.so *.a
	rm -rf $(TARGET_FQPROCESS) $(TARGET_TAGSORT) $(TARGET_FQSLIDESEQ) $(TARGET_FQMETRICS) $(TARGET_SAMPLEFQ)
