#Please see the tutorial (http://sourceforge.net/projects/dna-bison/files/bison_tutorial.tar.gz/download) for help in setting these options.
WORK=/home/ryand#This should be changed to match your needs
PREFIX = $(WORK)/bin
CC = mpicc
HTSLIB=htslib/libhts.a #Use -lhtslib and set LIB_DIRS for dynamic linkage
INCLUDE_DIRS = -Ihtslib
LIB_DIRS =
OPTS = -Wall -g #-DDEBUG #-DNOTHROTTLE
#MPI = -lmpich -lmpl #This is usually appropriate for mpich2
#MPI = #This is appropriate for mvapich2
MPI = -lmpi #This is usually appropriate for openmpi

#Don't edit below here unless you know what you're doing!

OBJS = aux.o fastq.o genome.o slurp.o master.o common.o MPI_packing.o worker.o
HERD_OBJS = herd/fastq.o herd/master.o herd/MPI_packing.o herd/slurp.o herd/worker.o herd/writer.o

.SUFFIXES:.c .o

all: HTSlib align index extractor mbias markduplicates

.c.o:
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) $< -o $@

HTSlib:
	$(MAKE) -C htslib

markduplicates: HTSlib markduplicates.o
	$(CC) $(LIB_DIRS) -o bison_markduplicates markduplicates.o $(HTSLIB) -lpthread -lz

mbias: HTSlib mbias.o
	$(CC) $(LIB_DIRS) -o bison_mbias mbias.o $(HTSLIB) -lpthread -lz

index:
	$(CC) $(OPTS) -o bison_index index.c -lpthread

align: HTSlib $(OBJS) main.o
	$(CC) $(LIB_DIRS) -o bison main.o $(OBJS) $(HTSLIB) -lm -lpthread $(MPI) -lz

extractor: HTSlib common.o methylation_extractor.o
	$(CC) $(LIB_DIRS) -o bison_methylation_extractor common.o methylation_extractor.o $(HTSLIB) -lpthread -lz

#Don't compile herd by default
herd:  HTSlib $(OBJS) $(HERD_OBJS) herd/main.o
	$(CC) $(LIB_DIRS) -o bison_herd herd/main.o $(OBJS) $(HERD_OBJS) $(HTSLIB) -lm -lpthread $(MPI) -lz

#Auxiliary programs, don't compile by default
auxiliary: merge_CpGs bedGraph2methylKit make_reduced_genome aux_python_scripts CpG_coverage bedGraph2MOABS

aux_python_scripts:
	cp -f auxiliary/bedGraph2BSseq.py ./
	cp -f auxiliary/merge_bedGraphs.py ./
	cp -f auxiliary/bedGraph2MethylSeekR.py ./

CpG_coverage: common.o auxiliary/CpG_coverage.o
	$(CC) $(LIB_DIRS) -o bison_CpG_coverage common.o auxiliary/CpG_coverage.o

merge_CpGs: common.o auxiliary/merge_CpGs.o
	$(CC) $(LIB_DIRS) -o bison_merge_CpGs common.o auxiliary/merge_CpGs.o

bedGraph2methylKit: common.o auxiliary/bedGraph2methylKit.o
	$(CC) $(LIB_DIRS) -o bedGraph2methylKit common.o auxiliary/bedGraph2methylKit.o

bedGraph2MOABS: common.o auxiliary/bedGraph2MOABS.o
	$(CC) $(LIB_DIRS) -o bedGraph2MOABS common.o auxiliary/bedGraph2MOABS.o

make_reduced_genome:
	$(CC) $(OPTS) $(LIB_DIRS) auxiliary/make_reduced_genome.c -o make_reduced_genome

install :
	mv bison_* $(PREFIX)/ 
	chmod a+x Rscripts/*
	cp Rscripts/* $(PREFIX)/
	if [ -f bison ]; then mv bison $(PREFIX)/ ; fi;
	if [ -f bedGraph2methylKit ]; then mv bedGraph2methylKit $(PREFIX)/ ; fi;
	if [ -f bedGraph2MOABS ]; then mv bedGraph2MOABS $(PREFIX)/ ; fi;
	if [ -f bedGraph2BSseq.py ]; then chmod a+x bedGraph2BSseq.py ; mv bedGraph2BSseq.py $(PREFIX)/ ; fi;
	if [ -f bedGraph2MethylSeekR.py ]; then chmod a+x bedGraph2MethylSeekR.py ; mv bedGraph2MethylSeekR.py $(PREFIX)/ ; fi;
	if [ -f merge_bedGraphs.py ]; then chmod a+x merge_bedGraphs.py ; mv merge_bedGraphs.py $(PREFIX)/ ; fi;
	if [ -f check_accuracy ]; then mv check_accuracy $(PREFIX)/ ; fi;
	if [ -f make_reduced_genome ]; then mv make_reduced_genome $(PREFIX)/ ; fi;

clean:
	rm -f *.o bison bison_* bedGraph2methylKit bedGraph2MOABS check_accuracy make_reduced_genome *.py
	rm -f herd/*.o
	rm -f auxiliary/*.o

clean-all: clean
	make --directory=htslib clean
