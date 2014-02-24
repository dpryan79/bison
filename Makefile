WORK=/home/ryand#This should be changed to match your needs
PREFIX = $(WORK)/bin
CC = mpicc
INCLUDE_DIRS = -I$(WORK)/include #This should be were samtools was compiled -I/path/to/samtools/compilation
LIB_DIRS = -L$(WORK)/lib #As above, but -L/path/to/samtools/compilation
OPTS = -Wall -O3 #-DDEBUG #-DNOTHROTTLE -g
MPI = -lmpich -lmpl #This is usually appropriate for mpich2
#MPI = #This is appropriate for mvapich2
#MPI = -lmpi #This is usually appropriate for openmpi

#Don't edit below here unless you know what you're doing!

OBJS = aux.o fastq.o genome.o slurp.o master.o common.o MPI_packing.o worker.o
HERD_OBJS = herd/fastq.o herd/master.o herd/MPI_packing.o herd/slurp.o herd/worker.o herd/writer.o

.SUFFIXES:.c .o

all: align index extractor mbias markduplicates

.c.o:
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) $< -o $@

markduplicates:
	$(CC) $(OPTS) $(INCLUDE_DIRS) $(LIB_DIRS) -o bison_markduplicates markduplicates.c -lpthread -lbam -lz

mbias:
	$(CC) $(OPTS) $(INCLUDE_DIRS) $(LIB_DIRS) -o bison_mbias mbias.c -lpthread -lbam -lz

index:
	$(CC) $(OPTS) -o bison_index index.c -lpthread

align: $(OBJS)
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) main.c -o main.o
	$(CC) $(OPTS) $(OBJS) main.o -o bison $(LIB_DIRS) -lm -lpthread $(MPI) -lbam -lz

extractor:
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) common.c -o common.o
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) methylation_extractor.c -o methylation_extractor.o
	$(CC) $(OPTS) $(LIB_DIRS) common.o methylation_extractor.o -o bison_methylation_extractor -lpthread -lbam -lz

#Don't compile herd by default
herd:  $(OBJS) $(HERD_OBJS)
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) herd/main.c -o herd/main.o
	$(CC) $(OPTS) $(OBJS) $(HERD_OBJS) herd/main.o -o bison_herd $(LIB_DIRS) -lm -lpthread $(MPI) -lbam -lz

#Auxiliary programs, don't compile by default
auxiliary:	merge_CpGs bedGraph2methylKit make_reduced_genome aux_python_scripts CpG_coverage

aux_python_scripts:
	cp -f auxiliary/bedGraph2BSseq.py ./
	cp -f auxiliary/merge_bedGraphs.py ./

CpG_coverage:	common.o
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) auxiliary/CpG_coverage.c -o auxiliary/CpG_coverage.o
	$(CC) $(OPTS) $(LIB_DIRS) common.o auxiliary/CpG_coverage.o -o bison_CpG_coverage

merge_CpGs:	common.o
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) auxiliary/merge_CpGs.c -o auxiliary/merge_CpGs.o
	$(CC) $(OPTS) $(LIB_DIRS) common.o auxiliary/merge_CpGs.o -o bison_merge_CpGs

bedGraph2methylKit:common.o
	$(CC) -c $(OPTS) $(INCLUDE_DIRS) auxiliary/bedGraph2methylKit.c -o auxiliary/bedGraph2methylKit.o
	$(CC) $(OPTS) $(LIB_DIRS) common.o auxiliary/bedGraph2methylKit.o -o bedGraph2methylKit

make_reduced_genome:
	$(CC) $(OPTS) $(LIB_DIRS) auxiliary/make_reduced_genome.c -o make_reduced_genome

install :
	mv bison_* $(PREFIX)/ 
	chmod a+x Rscripts/*
	cp Rscripts/* $(PREFIX)/
	if [ -f bison ]; then mv bison $(PREFIX)/ ; fi;
	if [ -f bedGraph2methylKit ]; then mv bedGraph2methylKit $(PREFIX)/ ; fi;
	if [ -f bedGraph2BSseq.py ]; then chmod a+x bedGraph2BSseq.py ; mv bedGraph2BSseq.py $(PREFIX)/ ; fi;
	if [ -f merge_bedGraphs.py ]; then chmod a+x merge_bedGraphs.py ; mv merge_bedGraphs.py $(PREFIX)/ ; fi;
	if [ -f check_accuracy ]; then mv check_accuracy $(PREFIX)/ ; fi;
	if [ -f make_reduced_genome ]; then mv make_reduced_genome $(PREFIX)/ ; fi;

clean:
	rm -f *.o bison bison_* bedGraph2methylKit check_accuracy make_reduced_genome bedGraph2BSseq.py
	rm -f herd/*.o
	rm -f auxiliary/*.o