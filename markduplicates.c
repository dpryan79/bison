#include <assert.h>
#include <inttypes.h>
#include <bam.h>

#define WORD_OFFSET(b) b/32
#define BIT_OFFSET(b) b%32

typedef struct {
    int32_t tid, start1, start2, stop1, stop2;
    int strand, MAPQ;
    unsigned read_number;
} alignment;

typedef struct {
    uint64_t nelements;
    int threadid;
    alignment *alignments;
} qsort_func_struct;

/*
    Sort a list of alignments, they'll be ordered as follows (always low to high):
    (1) tid (chromosome index ID)
    (2) start1 (5' position of read #1)
    (3) stop1 (3' position of read #1)
    (4) start2 (5' position of read #2)
    (5) stop2 (3' position of read #2)
    (6) strand (0 OT, 1 OB, 2 CTOT, 3 CTOB)
    (7) MAPQ (high to low)
*/
int comp_func(const void *a, const void *b) {
    alignment *a1 = (alignment*) a;
    alignment *a2 = (alignment*) b;

    if(a1->tid < a2->tid) return -1;
    else if(a1->tid > a2->tid) return 1;
    else {
        if(a1->start1 < a2->start1) return -1;
        else if(a1->start1 > a2->start1) return 1;
        else {
            if(a1->stop1 < a2->stop1) return -1;
            else if(a1->stop1 > a2->stop1) return 1;
            else {
                if(a1->start2 < a2->start2) return -1;
                else if(a1->start2 > a2->start2) return 1;
                else {
                    if(a1->stop2 < a2->stop2) return -1;
                    else if(a1->stop2 > a2->stop2) return 1;
                    else {
                        if(a1->strand < a2->strand) return -1;
                        else if(a1->strand > a2->strand) return 1;
                        else { //This will be the other way around
                            if(a1->MAPQ > a2->MAPQ) return -1; 
                            else if(a1->MAPQ < a2->MAPQ) return 1; 
                            else return 0;
                        }
                    }
                }
            }
        }
    }
}

//This is the same as comp_func(), except that MAPQ is ignored
int comp_func2(const void *a, const void *b) {
    alignment *a1 = (alignment*) a;
    alignment *a2 = (alignment*) b;

    if(a1->tid < a2->tid) return -1;
    else if(a1->tid > a2->tid) return 1;
    else {
        if(a1->start1 < a2->start1) return -1;
        else if(a1->start1 > a2->start1) return 1;
        else {
            if(a1->stop1 < a2->stop1) return -1;
            else if(a1->stop1 > a2->stop1) return 1;
            else {
                if(a1->start2 < a2->start2) return -1;
                else if(a1->start2 > a2->start2) return 1;
                else {
                    if(a1->stop2 < a2->stop2) return -1;
                    else if(a1->stop2 > a2->stop2) return 1;
                    else {
                        if(a1->strand < a2->strand) return -1;
                        else if(a1->strand > a2->strand) return 1;
                        else return 0;
                    }
                }
            }
        }
    }
}

void *qsort_func(void *a) {
    uint64_t total_pairs = ((qsort_func_struct *) a)->nelements;
    int thread_id = ((qsort_func_struct *) a)->threadid;
    uint64_t nelements = total_pairs/(thread_id+1);
    uint64_t offset = thread_id*nelements;
    alignment *alignments = ((qsort_func_struct *) a)->alignments;
    void *p = (void *) (alignments+offset);

    if(thread_id == 5) nelements += total_pairs % thread_id;
    qsort(p, (size_t) nelements, sizeof(alignment), &comp_func);

    return NULL;
}

//Set 1 at the given offset
void set_bit(uint32_t *map, uint64_t n) { 
    map[WORD_OFFSET(n)] |= (1 << BIT_OFFSET(n));
}

//Get the value at a given offset
int get_bit(uint32_t *map, uint64_t n) {
    uint32_t bit = map[WORD_OFFSET(n)] & (1 << BIT_OFFSET(n));
    return bit != 0; 
}

uint64_t mark_dups(alignment *alignments, uint32_t *bitmap, uint64_t total_pairs) {
    uint64_t i, ndups = 0;
    void *cur_alignment = (void *) alignments;

    for(i=1;i<total_pairs;i++) {
        if(comp_func2(cur_alignment, (void *) (alignments+i)) == 0) {
            set_bit(bitmap, alignments[i].read_number);
            ndups++;
        } else {
            cur_alignment = (void *) (alignments+i);
        }
    }
    return ndups;
}
        

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] input.bam output.bam\n", prog);
    printf("\n\
This program will parse a BAM file produced by bison and mark likely PCR\n\
duplicates in a new file. A PCR duplicate is defined as two reads/read pairs\n\
having the same start and stop coordinates on the some strand of the same\n\
chromosome. The read or pair with the best MAPQ score will be kept.\n\
\n\
There are better ways to go about this (both in terms of algorithms and in the\n\
information used to determine a PCR duplicate), but this will normally suffice.\n\
\n\
N.B., this program does not currently support coordinate-sorted BAM files or\n\
paired-end files containing discordant or other mixed alignments where the reads\n\
are on different chromosomes or different strand.\n");
    printf("\nOptions:\n\
\n\
    -s INT Initial array size used to hold mapping coordinates (measured in\n\
           reads/read pairs). The default is 10 million. In an ideal world, this\n\
           would be at least as large as the number of reads you have, but the\n\
           array will grow as needed.\n\
\n\
    -g INT How much the array should grow (measured in reads/read pairs) when\n\
           it's full). The default is 1000000.\n");
}

inline int get_strand(bam1_t *read) {
    char *XG = bam_aux2Z(bam_aux_get(read, "XG"));
    char *XR = bam_aux2Z(bam_aux_get(read, "XR"));

    if(*XG == 'C') { //OT or CTOT
        if(*XR == 'C') return 0; //OT
        else return 2; //CTOT
    } else {
        if(*XR == 'C') return 1; //OB
        else return 3; //CTOB
    }
}

int main(int argc, char *argv[]) {
    bam1_t *read;
    bam_header_t *header;
    bamFile fp = NULL, of = NULL;
    uint64_t total_pairs = 0, max_length = 10000000;
    uint64_t bitmap_length = 0, cur_read = 0, ndups = 0;
    uint64_t grow_size = 1000000;
    alignment *alignments;
    uint32_t *bitmap;
    int i;
    char *iname = NULL, *oname = NULL;

    if(argc < 3) {
        usage(argv[0]);
        return 0;
    }

    for(i=1; i<argc; i++) {
        if(strcmp("-h", argv[i]) == 0 || strcmp("--help", argv[i]) == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp("-s", argv[i]) == 0) {
            max_length = (uint64_t) strtoull(argv[++i], NULL, 10);
        } else if(strcmp("-g", argv[i]) == 0) {
            grow_size = (uint64_t) strtoull(argv[++i], NULL, 10);
        } else if(iname == NULL) {
            iname = argv[i];
            fp = bam_open(iname, "r");
        } else if(oname == NULL) {
            oname = argv[i];
            of = bam_open(oname, "w");
        } else {
            printf("Unrecognized option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }
    if(iname == NULL || oname == NULL) {
        printf("Either the input or output file name were not specified!\n");
        usage(argv[0]);
        return 2;
    }

    //Set everything up
    header = bam_header_read(fp);
    read = bam_init1();
    alignments = malloc(sizeof(alignment)*max_length);
    assert(alignments != NULL);

    //Read in the alignments
    while(bam_read1(fp, read) > 1) {
        alignments[total_pairs].tid = read->core.tid;
        alignments[total_pairs].start1 = read->core.pos;
        alignments[total_pairs].stop1 = bam_calend(&(read->core), bam1_cigar(read));
        alignments[total_pairs].MAPQ = read->core.qual;
        alignments[total_pairs].strand = get_strand(read);
        if(read->core.flag & BAM_FPAIRED) {
            assert(bam_read1(fp,read)>1);
            alignments[total_pairs].start2 = read->core.pos;
            alignments[total_pairs].stop2 = bam_calend(&(read->core), bam1_cigar(read));
        } else {
            alignments[total_pairs].start2 = 0;
            alignments[total_pairs].stop2 = 0;
        }
        alignments[total_pairs].read_number = total_pairs;
        total_pairs++;

        //Lengthen the array as needed
        if(max_length-total_pairs == 0) {
            max_length += grow_size;
            alignments = realloc(alignments, max_length * sizeof(alignment));
            assert(alignments != NULL);
        }
    }
    bam_close(fp);

    //create bitmap
    bitmap_length = max_length/32;
    bitmap_length += (max_length % 32 > 0) ? 1 : 0;
    bitmap = calloc(bitmap_length, sizeof(uint32_t));

    //Sort
    qsort((void *) alignments, (size_t) total_pairs, sizeof(alignment), &comp_func);

    //Mark duplicates in bitmap
    ndups = mark_dups(alignments, bitmap, total_pairs);
    free(alignments);
    printf("There were %"PRIu64" duplicates from %"PRIu64" total reads or pairs\n", ndups, total_pairs);

    //reopen file, iterate through and change flags as appropriate
    fp = bam_open(iname, "r");
    header = bam_header_read(fp);
    of = bam_open(oname, "w");
    bgzf_mt(of, 4, 256); //This should be user configurable
    bam_header_write(of, header);

    while(bam_read1(fp, read) > 1) {
        if(get_bit(bitmap, cur_read)) read->core.flag = read->core.flag | BAM_FDUP;
        bam_write1(of, read);
        if(read->core.flag & BAM_FPAIRED) {
            assert(bam_read1(fp, read) > 1);
            if(get_bit(bitmap, cur_read)) read->core.flag = read->core.flag | BAM_FDUP;
            bam_write1(of, read);
        }
        cur_read++;
    }

    //Clean up
    bam_close(fp);
    bam_close(of);
    bam_destroy1(read);
    free(bitmap);

    return 0;
}
