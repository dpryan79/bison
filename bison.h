#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>
#include <ctype.h>

#include <zlib.h>
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include <limits.h>
#include <inttypes.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#define MAXREAD 1024
#define MASTER 0
#define VERSION "0.4.0"
#define BT2BUF_SZ 256 * 1024
#define THROTTLE_CHECK_INTERVAL 100000 //When bison_herd auto-throttles, this specifies how frequently it should check whether it should do so (units are "reads")
#define version() printf("Bison, version %s\n", VERSION)

/******************************************
*
* MPI Send/Recv tags:
*
* 0: Workers should start
* 1: Header size (this could be removed)
* 2: packed header struct
* 3: Packed fastq struct
* 4: Unused (used to be packed read size)
* 5: Packed read
*
******************************************/

//Mutexes for thread i/o designation. A thread should not read/write until it's ID number is equal to these
FILE *zip1;
FILE *zip2;
FILE *unmapped1;
FILE *unmapped2;
htsFile *OUTPUT_BAM;
unsigned long long t_reads; //total number of reads
unsigned long long t_concordant; //Total concordant pairs
unsigned long long t_discordant; //Total discordant pairs
unsigned long long t_singletons; //Total singletons
unsigned long long m_reads_OT; //total number mapped to the OT strand
unsigned long long m_reads_OB;
unsigned long long m_reads_CTOT;
unsigned long long m_reads_CTOB;
unsigned long long t_CpG; //Total CpGs
unsigned long long m_CpG; //Methylated CpGs
unsigned long long t_CHG;
unsigned long long m_CHG;
unsigned long long t_CHH;
unsigned long long m_CHH;

//This is useful for single-node debugging
#ifdef DEBUG
int global_debug_taskid;
htsFile *fp1;
htsFile *fp2;
htsFile *fp3;
htsFile *fp4;
#endif

//Some people may find it useful to have the system throttle itself so as not to overwhelm the MPI buffer
unsigned long long *nwritten;

//Mutex for controlling access to the global metrics struct and the output files
pthread_mutex_t metrics_mutex;

typedef struct {
    char *chrom;
    unsigned long long offset;
    unsigned long long length;
} chromosome_struct;

typedef struct {
    int nchromosomes;
    unsigned long long max_genome; //The offset value in the read_genome() function will be used to keep track of how close we are
    chromosome_struct **chromosome;
    char *genome; //This will hold the genomic sequence in memory as a continuous string.
} chromosomes_struct;

typedef struct {
    char *FASTQ1;
    char *FASTQ2;
    char *FASTQ1CT;
    char *FASTQ1GA;
    char *FASTQ2CT;
    char *FASTQ2GA;
    char *unmapped1;
    char *unmapped2;
    char *genome_dir;
    char *basename;
    char *fai;
    char *odir;
    char *tmpdir;
    char *bowtie2_options;
    char *outname;
    char scoremin_type;
    int paired;
    int directional;
    int nthreads;
    int nmthreads;
    int buffer_size;
    int send_receive_buffer_size;
    int unmapped;
    int mode; //0 is --end-to-end (default), 1 is local
    int quiet; //0 or 1, the latter supresses all output to the screen
    int reorder; //0 or 1, latter reorders writing to match input, only meaningful in herd
    int n_compression_threads; //Default is 0
#ifndef NOTHROTTLE
    int reads_in_queue;
#endif
    float scoremin_intercept;
    float scoremin_coef;
    int isCRAM;
    char **argv;
    int argc;
    uint64_t maxMem;
    int sort; //To sort or not
} t_config;

typedef struct {
    int size; //Theres an effective size limit imposed by MPI of whatever int is
    void *packed;
} MPI_Header;

typedef struct {
    int size;
    void *packed; //the format is sizeof(bam1_t) followed by data, which is of size data_len
} MPI_read;

typedef struct {
    int size;
    void *packed; //format is: char *name1\0seq1\0qual1\0 followed by optional char *name2\0seq2\0qual2\0
} MPI_Fastq;

typedef struct {
    int max_name1; //current maximum length of memory for name1
    int max_seq1;
    int max_qual1;
    int max_name2;
    int max_seq2;
    int max_qual2;
    char *name1;
    char *seq1;
    char *qual1;
    char *name2;
    char *seq2;
    char *qual2;
} fastq;

//This is used as the input struct for slurp_fastq
typedef struct {
    int thread_id;
    char *fastq1;
    char *fastq2;
} slurp_fastq_struct;

struct packed_struct {
    void *packed; //If NULL, then finished
    struct packed_struct *next;
    struct packed_struct *previous; //Only used on last sentinel struct
    int state; //0    no next node (not ready)
               //1    has next node (ready)
};

//Global values
t_config config;
bam_hdr_t *global_header;
//This will be the global structure for pointers to chromosome_struct's holding the information for *genome
chromosomes_struct chromosomes;
char **fnames1, **fnames2; //This will hold the file names so that the writer thread knows what to rename things
unsigned long long *flengths; //This will hold the size of each file

//Linked-list of reads
struct packed_struct *node1, *node1_last_sentinel;
struct packed_struct *node2, *node2_last_sentinel;
struct packed_struct *node3, *node3_last_sentinel;
struct packed_struct *node4, *node4_last_sentinel;
//bison-herd
struct packed_struct **nodes, **last_sentinel_node;
struct packed_struct **fastq_nodes, **last_fastq_sentinel_node;
struct packed_struct **to_write_node, **to_write_sentinel_node;

/******************************************************************************
*
*   Take a fastq struct and convert it G->A, the conversion is in place
*
*   fastq *read, input struct
*   int which, which of the reads to convert
*
*******************************************************************************/
void convertGA(fastq *, int); //fastq.c

/******************************************************************************
*
*   Take a fastq struct and convert it C->T, the conversion is in place
*
*   fastq *read, input struct
*   int which, which of the reads to convert
*
*******************************************************************************/
void convertCT(fastq *, int ); //fastq.c

/******************************************************************************
*
*   Write an unmapped read to a gzipped fastq file.
*
*   FILE *fp: gzipped fastq file
*   bam1_t *read: read to write in fastq format
*
*******************************************************************************/
void write_unmapped(FILE *, bam1_t *); //fastq.c

/******************************************************************************
*
*   Read in the fastq file(s) sending the reads to the appropriate nodes and
*   also storing the unconverted reads in a linked list on the master node.
*
*   This will act as its own thread on the master node.
*
*   void *a is unused but required by pthreads
*
*******************************************************************************/
void * send_store_fastq(void *);

/******************************************************************************
*
*   Add an element to the end of a linked-list
*
*   struct packed_struct *last: last sentinel struct
*   void *packed: a packed read
*
*******************************************************************************/
void add_element(struct packed_struct *, void *);

/******************************************************************************
*
*   Remove an element from the start of a linked-list
*   is_ready(first, 0) must return 1!
*
*   struct packed_struct *: first sentinel struct
*
*******************************************************************************/
void remove_element(struct packed_struct *); //slurp.c
//As above, but the packed component can't have been updated to a bam1_t
void remove_raw_element(struct packed_struct *); //slurp.c

/******************************************************************************
*
*   Move an element from one linked-list to another.
*
*   struct packed_struct *source: source linked list
*   struct packed-struct *dest: destination sentinel node
*
*******************************************************************************/
void move_element(struct packed_struct *, struct packed_struct *);

/******************************************************************************
*
*   Is the first or second element ready?
*
*   struct packed_struct *first: first sentinel struct
*   int offset: 0 (first element) or 1 (second element)
*
*    returns 1 for element ready, or 0 otherwise
*
*******************************************************************************/
int is_ready(struct packed_struct *, int); //slurp.c

/******************************************************************************
*
*   Is the linked list finished?
*
*   struct packed_struct *: first sentinel struct
*
*   returns 1 for finished, 0 otherwise
*
*******************************************************************************/
int is_finished(struct packed_struct *); //slurp.c

/******************************************************************************
*
*   Add an elemnt to a node designating that the list is finished
*
*   struct packed_struct: last sentinel node
*
*******************************************************************************/
void add_finished(struct packed_struct *); //slurp.c

/******************************************************************************
*
*   Initialize a linked list, returning the last sentinel struct
*
*   struct packed_struct *: first sentinel struct
*
*   returns first sentinel struct
*
*******************************************************************************/
struct packed_struct *initialize_list(struct packed_struct *); //slurp.c

/******************************************************************************
*
*   Destroy a linked list of packed_structs
*
*   struct packed_struct *first: linked list to destroy
*
*******************************************************************************/
void destroy_list(struct packed_struct *); //slurp.c
//As above, but for lists where ->packed hasn't been converted to a bam1_t
void destroy_raw_list(struct packed_struct *); //slurp.c

/******************************************************************************
*
*   The MPI receiver thread on the main node
*
*   void *: NULL input
*
*   returns NULL
*
*******************************************************************************/
void *slurp(void *); //slurp.c
void *herd_slurp(void *); // herd/slurp.c

/******************************************************************************
*
*   Construct the output directory name, putting it in config.odir
*
*******************************************************************************/
void update_odir(); //fastq.c

/******************************************************************************
*
*   Given the name of a (possibly gzipped) fastq file, return the file name
*   with the .fastq.gz, .fq.gz, .fastq, or .fq extension removed.
*
*   char *file: filename
*
*   CAUTION, THE OUTPUT MUST BE free()d!
*
*******************************************************************************/
char * get_basename(char *); //fastq.c

/******************************************************************************
*
*   Invoke the C->T and G->A conversion threads of the fastq files (located in
*   the global config structure).
*
*   FLAGS: integer bit field denoting the conversions to make
*       0x8 fastq #1 C->T
*       0x4 fastq #1 G->A
*       0x2 fastq #2 C->T
*       0x1 fastq #2 G->A
*
*******************************************************************************/
void convert_fastq(int, unsigned int); //fastq.c

/******************************************************************************
*
*   Take the config.FASTQ1 and config.FASTQ2 filenames and use them to generate
*   the config.FASTQ1CT... filenames. These must subsequently be free()d, which
*   is done in the quit() function.
*
*   char *f1, config.FASTQ1
*   char *f2, config.FASTQ2
*   these are only really needed if there's more than one input file
*
*******************************************************************************/
void create_fastq_names(char *, char*); //fastq.c

/******************************************************************************
*
*   Read in all .fa and .fasta files within config.genome_dir. The sequences 
*   are concatenated onto chromosomes.genome. The global chromosomes structure
*   is modified with each new chromosome.
*
*   Note, chromosomes.genome (in fact, all of chromosomes) need to be free()d
*   The is performed by the quit() function.
*
*******************************************************************************/
void read_genome(); //common.c

/*******************************************************************************
*
*  Replace the @PG line in the header with the actual comman executed
*
*******************************************************************************/
bam_hdr_t *modifyHeader(bam_hdr_t *hdr, int argc, char *argv[]);

/******************************************************************************
*
*   Print metrics to STDOUT and a file.
*
*******************************************************************************/
void print_metrics(); //aux.c

/******************************************************************************
*
*   Return the number of worker nodes that will actually run.
*
*******************************************************************************/
int effective_nodes(); //aux.c

/******************************************************************************
*
*   quit, while performing some cleanup
*
*   int FLAG: What to free/close/etc.
*             0x1 things created by create_fastq_names()
*             0x2 things pthreads are closed and bam headers destroyed
*             In addition, the master node will free chromosomes.genome, close
*             the BAM file, and free everything in the chromosomes struct.
*             Also, everynode will free config.bowtie2_options
*
*   int rv: return value
*
*******************************************************************************/
void quit(int, int); //aux.c

/******************************************************************************
*
*   Take a BAM header and pack it into a single contiguous memory block. Store
*   the resulting block and its size in an MPI_Header structure.
*
*   THE RESULT MUST BE free()d
*
*   bam_header_t *header: The header to store
*
*******************************************************************************/
MPI_Header * pack_header(bam_hdr_t *); //MPI_packing.c

/******************************************************************************
*
*   Unpack a header packed into an initialized bam_header_t
*
*   bam_header_t *header: The header to unpack into
*   void *packed: The packed header
*
*******************************************************************************/
void unpack_header(bam_hdr_t *, void *); //MPI_packing.c

/******************************************************************************
*
*   Take a fastq struct and pack it for shipping
*
*   THE RESULT MUST BE free()d eventually
*
*   fastq *read: The read(s) to store
*   MPI_Fastq *output: the struct into which to pack things
*
*******************************************************************************/
MPI_Fastq * pack_fastq(fastq *); //MPI_packing.c

/******************************************************************************
*
*   Take unpack a packed fastq struct
*
*   THE RESULT MUST BE free()d
*
*   fastq *read: The fastq struct to unpack into
*   void *packed: The packed structure
*
*******************************************************************************/
fastq * unpack_fastq(fastq *, void *); //MPI_packing.c

/******************************************************************************
*
*   Unpack a packed read into an initialized bam1_t read.
*
*   bam1_t *read: The read to unpack into
*   void *packed: The packed read
*
*******************************************************************************/
bam1_t *unpack_read(bam1_t *, void *); //MPI_packing.c

/******************************************************************************
*
*   Take a BAM read and pack it into a single contiguous memory block. Store
*   the resulting block and its size in an MPI_Read structure.
*
*   THE RESULT MUST BE free()d
*
*   bam1_t *read: The read to store
*
*******************************************************************************/
MPI_read * pack_read(bam1_t *, MPI_read *); //MPI_packing.c

/******************************************************************************
*
*   Extract the next sequence line from a file stream.
*
*   char *seq: destination
*   FILE *fp: source
*
*   THE OUTPUT MUST BE free()d
*   This function is affected by the MAXREAD definition, above. If this value is
*   less than the longest read, things will break. It would be better to realloc
*   as needed.
*
*******************************************************************************/
void get_seq(char *, FILE *); //genome.c

/******************************************************************************
*
*   Reverse complement a sequence (in place)
*
*   char *seq: the sequence
*
*******************************************************************************/
void reverse_complement(char *); //common.c

/******************************************************************************
*
*   Determine the appropriate offset in chromosomes.genome
*
*   char *chrom: Chromosome name
*   int32_t pos: 0-based position on Chromosome. This is read->core.pos
*
*******************************************************************************/
unsigned long long genome_offset(char*, int32_t); //common.c

/******************************************************************************
*
*   Return the length of a given chromosome.
*
*   char *chrom: the chromosome of interest
*
*******************************************************************************/
unsigned long long genome_chrom_length(char *); //genome.c

/******************************************************************************
*
*   Return a pointer to the chromosome name onto which a read maps.
*
*   bam1_t *read: The read in question
*
*******************************************************************************/
char *lookup_chrom(bam1_t *); //common.c

/******************************************************************************
*
*   Return a base and another 2 bases on one of its sides. This is needed for
*   making methylation calls. If this span goes off the edge of a chromosome,
*   N's will be used.
*
*   unsigned long long offset: from genome_offset
*   unsigned long long position: converted read->core.pos
*   int change: Direction of the context (- is backwards)
*   unsigned long long chrom_length: from genome_chrom_length
*
*   The output needs to be free()d
*
*******************************************************************************/
char* get_genomic_context(unsigned long long, unsigned long long, int, unsigned long long); //genome.c

/*******************************************************************************
*
*  Create a position array to account for any InDels
*  This function assumes that the first base is not marked as an InDel or
*  clipped in any way. If that occurs then things will break.
*
*  The output needs to be free()d
*
*******************************************************************************/
unsigned long long *calculate_positions(bam1_t *); //common.c

/*******************************************************************************
*
*   The master node function.
*
*   void *a: Actually an int*, the thread_id
*
*******************************************************************************/
void * master_processer_thread(void*); //master.c
void * herd_master_processer_thread(void*); //master.c under herd/

/******************************************************************************
*
*   Given a set of single-end reads, determine which one, if any, aligns best.
*   Then, add the various XM/XX/etc. tags and prepare the read for writing. The
*   final read will always be stored in read1. Return the worker node number
*   producing the best alignment (or 0).
*
*   bam1_t *readN: Unpacked reads from the worker nodes
*   char *seq: The unconverted fastq read
*
*******************************************************************************/
int32_t process_single(bam1_t *, bam1_t *, bam1_t *, bam1_t *, char *); //master.c

/******************************************************************************
*
*   Like process_single, but for paired_end reads. The bam1_t**s hold the
*   buffered reads. i denotes the read#1 of interest (read #2 is the next read)
*
*******************************************************************************/
int32_t process_paired(bam1_t **, bam1_t **, bam1_t **, bam1_t **, char **); //master.c

/*******************************************************************************
*
*   Update a packed read so that it's a proper bam1_t and return a pointer
*
*   struct packed_struct *first: first sentinel node
*   int offset: Return the read from the first (0) or second (1) element
*
*   returns a pointer to a bam1_t read
*
*******************************************************************************/
bam1_t *update_read(struct packed_struct *, int); //master.c

/******************************************************************************
*
*   This function will run as its own thread and process the linked lists
*   output from the master processor threads, writing them in order to a BAM
*   file. This will also write all of the other output (aside from metrics).
*   Furthermore, this provides a readout of the current number of reads
*   processed.
*
*   Output is NULL, as is the input (needed by pthreads).
*
*******************************************************************************/
void * bam_writer(void *); //writer.c

/******************************************************************************
*
*   This receives the reads, converts them, and writes them to the FIFO(s)
*   
*   void *a: a pointer to a struct with the following components:
*
*   int thread_id: the thread_id
*   char *fastq1: FIFO from which bowtie2 can get read1
*   char *fastq2: FIFO from which bowtie2 can get read2 (if it exists)
*
*******************************************************************************/
void * slurp_fastq(void *); //worker.c

/******************************************************************************
*
*   The main worker node function.
*
*   int thread_id: the thread_id
*
*******************************************************************************/
void worker_node(int); //worker.c

/******************************************************************************
*
*   The main worker node function.
*
*   int thread_id: the thread_id
*   char *fastq1: FIFO from which bowtie2 can get read1
*   char *fastq2: FIFO from which bowtie2 can get read2 (if it exists)
*
*******************************************************************************/
void herd_worker_node(int, char *, char *); //worker.c under herd/

/******************************************************************************
*
*   Open a sam file for reading via popen
*
*   char *cmd: The command given to popen, the mode is always "r".
*
*******************************************************************************/
htsFile * sam_popen(char *); //aux.c

/******************************************************************************
*
*   Close a SAM file that was opened with sam_popen
*
*   tamFile fp: The file pointer struct returned from sam_popen
*
*******************************************************************************/
void sam_pclose(samFile *fp);

//Sorting stuff
typedef struct {
    uint32_t l,m; //Current and max buffer length
    uint64_t curMem, maxMem; //Current and maximum permitted memory usage
    int offset; //File number offset for the temp files
    char *opref; //output file name prefix
    bam1_t **buf; //alignment buffer
} alignmentBuffer;

/******************************************************************************
*
*   Push an alignment onto an alignmentBuffer struct 
*
*   alignmentBuffer *buf: the alignmentBuffer
*   bam1_t *b: the alignment
*
*   returns the modified buffer. If we would have hit maxMem, then the buffer
*   is written out to temp files as appropriate.
*
*******************************************************************************/
alignmentBuffer *pushAlignmentBuffer(alignmentBuffer *buf, bam1_t *b);

/******************************************************************************
*
*   Merge temp files created by an alignment buffer. This will write any
*   remaining alignments to temp files.
*
*   alignmentBuffer *buf: the alignmentBuffer
*
*   Note that the buffer is modified so it can be resused.
*
*******************************************************************************/
void mergeTemp(alignmentBuffer *buf);

/*******************************************************************************
*
*   Parse a memory amount specified like 2G, or 1.5M
*
*   char *str: The string to parse
*
*   Returns the number of bytes represented by the string
*******************************************************************************/
uint64_t str2Mem(char *s);
