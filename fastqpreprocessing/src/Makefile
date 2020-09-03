IDIR =../libStatGen/include

CC = g++ -std=c++17 -O4
CFLAGS = -I$(IDIR) -L../libStatGen

ODIR=obj
LDIR =../libStatGen/

TARGET = fastqprocess
LIBS = -lStatGen -lz -lpthread

_DEPS = fastqprocess.h utilities.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = fastqprocess.o  utilities.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)  $(_DEPS)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
