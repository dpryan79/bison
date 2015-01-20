#include "bison.h"

#define WORD_OFFSET(b) (b)/32
#define BIT_OFFSET(b) (b)%32
#define GET_STRAND(b) (b) >> 30
#define MIN_PHRED 5 //Minimum phred score that a basecall must have to be included in the methylation call bitmap
#define min(a,b) (a<=b) ? a : b

//Lookup table for popcount, 16 bits so 4 function calls/value
unsigned bitcounts_table[65536];

//The following is only used for filling in the table above
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

unsigned popcount(uint64_t x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

void fill_table() {
    uint64_t i;
    for(i=0;i<65536;i++) bitcounts_table[i] = popcount(i);
}

typedef struct {
    int32_t tid1, tid2, start1, start2, total_phred; //First 2 bits of tid denote strand: OT 00, OB 01, CTOT 10, CTOB 11, so ~1 billion contigs
    uint64_t meth[8], unmeth[8];
    uint32_t read_number;
} alignment;

typedef struct {
    uint64_t nelements;
    int threadid;
    alignment *alignments;
} qsort_func_struct;

//Get the total phred score
inline int32_t total_phred(bam1_t *read) {
    int32_t l = read->core.l_qseq;
    int32_t i=0, total=0;
    uint8_t *qual = bam_get_qual(read);

    for(i=0; i<l; i++) {
        total += *(qual+i);
    }
    return total;
}

//Get the total number of bits set in a uint64_t
int lookup_popcount(uint64_t c) {
    int total, i;

    i = c & 0x000000000000FFFF;
    total = bitcounts_table[i];
    i = (c & 0x00000000FFFF0000) >> 16;
    total += bitcounts_table[i];
    i = (c & 0x0000FFFF00000000) >> 48;
    total += bitcounts_table[i];
    i = (c & 0xFFFF000000000000) >> 48;
    total += bitcounts_table[i];
    return total;
}

/*
    Sort a list of alignments, they'll be ordered as follows (always low to high):
    (1) tid1 (chromosome index ID), which is also strand
    (2) tid2 (chromosome index ID), which is also strand
    (3) start1 (5' position of read #1)
    (4) start2 (5' position of read #2)
    (5) strand (0 OT, 1 OB, 2 CTOT, 3 CTOB)
    (6) total calls, aka Hamming weight (high to low)
    (7) total phred (high to low)
*/
int comp_func(const void *a, const void *b) {
    alignment *a1 = (alignment*) a;
    alignment *a2 = (alignment*) b;
    int ncalls1, ncalls2;

    if(a1->tid1 < a2->tid1) return -1; //tid1, strand
    else if(a1->tid1 > a2->tid1) return 1;
    else {
        if(a1->tid2 < a2->tid2) return -1; //tid2, strand
        else if(a1->tid2 > a2->tid2) return 1;
        else {
            if(a1->start1 < a2->start1) return -1; //start1
            else if(a1->start1 > a2->start1) return 1;
            else {
                if(a1->start2 < a2->start2) return -1; //start2
                else if(a1->start2 > a2->start2) return 1;
                else {
                    ncalls1 = lookup_popcount(a1->meth[0]) + lookup_popcount(a1->meth[1]) \
                        + lookup_popcount(a1->meth[2]) + lookup_popcount(a1->meth[3]) \
                        + lookup_popcount(a1->meth[4]) + lookup_popcount(a1->meth[5]) \
                        + lookup_popcount(a1->meth[6]) + lookup_popcount(a1->meth[7]) \
                        + lookup_popcount(a1->unmeth[0]) + lookup_popcount(a1->unmeth[1]) \
                        + lookup_popcount(a1->unmeth[3]) + lookup_popcount(a1->unmeth[3]) \
                        + lookup_popcount(a1->unmeth[4]) + lookup_popcount(a1->unmeth[5]) \
                        + lookup_popcount(a1->unmeth[6]) + lookup_popcount(a1->unmeth[7]);
                    ncalls2 = lookup_popcount(a2->meth[0]) + lookup_popcount(a2->meth[1]) \
                        + lookup_popcount(a2->meth[2]) + lookup_popcount(a2->meth[3]) \
                        + lookup_popcount(a2->meth[4]) + lookup_popcount(a2->meth[5]) \
                        + lookup_popcount(a2->meth[6]) + lookup_popcount(a2->meth[7]) \
                        + lookup_popcount(a2->unmeth[0]) + lookup_popcount(a2->unmeth[1]) \
                        + lookup_popcount(a2->unmeth[3]) + lookup_popcount(a2->unmeth[3]) \
                        + lookup_popcount(a2->unmeth[4]) + lookup_popcount(a2->unmeth[5]) \
                        + lookup_popcount(a2->unmeth[6]) + lookup_popcount(a2->unmeth[7]);
                    if(ncalls1 > ncalls2) return -1; //total calls
                    else if(ncalls1 < ncalls2) return 1;
                    else {
                        if(a1->total_phred > a2->total_phred) return -1; //sum of phred scores, this is backwords
                        else if(a1->total_phred < a2->total_phred) return 1;
                        else return 0;
                    }
                }
            }
        }
    }
}

//This is the same as comp_func(), with only tid1, tid2, start1, start2 and strand used
int comp_func2(alignment *a1, alignment *a2) {
    if(a1->tid1 < a2->tid1) return -1;
    else if(a1->tid1 > a2->tid1) return 1;
    else {
        if(a1->tid2 < a2->tid2) return -1;
        else if(a1->tid2 > a2->tid2) return 1;
        else {
            if(a1->start1 < a2->start1) return -1;
            else if(a1->start1 > a2->start1) return 1;
            else {
                if(a1->start2 < a2->start2) return -1;
                else if(a1->start2 > a2->start2) return 1;
                else return 0;
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

//Return the number of reads to compare for PCR duplicates, based on tid,
//start1, and start2
uint64_t get_bin_limits(alignment *alignments, uint64_t total) {
    uint64_t i, n = 1;
    alignment *next_alignment;

    for(i=1; i<total; i++) {
        next_alignment = alignments+i;
        if(comp_func2(alignments, next_alignment) != 0) break;
        n++;
    }
    return n;
}

//Return the edit distance between two methylation call string. The edit distance is unaffected by call/no-call differences in the same position
/*
    The general idea is:
    (1) XOR the METH calls
    (2) XOR the UNMETH calls
    (3) AND the results
    (4) the popcount of the result is the edit distance
*/
int edit_distance(alignment *a1, alignment *a2) {
    int distance = 0, i;
    uint64_t diff1, diff2;

    for(i=0; i<8; i++) {
        diff1 = a1->meth[i] ^ a2->meth[i];
        diff2 = a1->unmeth[i] ^ a2->unmeth[i];
        distance += lookup_popcount(diff1 & diff2);
    }

    return distance;
}

//Mark the duplicates in a bin
uint64_t mark_bin(alignment *alignments, uint64_t i, uint64_t bin_width, uint32_t *bitmap) {
    uint64_t j, k, ndups=0;

    for(j=i; j<i+bin_width-1; j++) { //There's no point in starting on the last one
        //If this read is already a duplicate, then ignore
        if(get_bit(bitmap, alignments[j].read_number) == 0) {
            for(k=j+1; k<j+bin_width; k++) {
                if(get_bit(bitmap, alignments[k].read_number) == 0) {
                    if(edit_distance(alignments+j, alignments+k) == 0) {
                        set_bit(bitmap, alignments[k].read_number);
                        ndups++;
                    }
                }
            }
        }
    }
    return ndups;
}
    
uint64_t mark_dups(alignment *alignments, uint32_t *bitmap, uint64_t total_pairs) {
    uint64_t i = 0, ndups = 0, bin_width;

    while(i<total_pairs) {
        /*
            Process the alignments in bins of possible duplicates (by tid1, tid2, start1, start2)
            Within each bin, mark duplicates.
        */
        bin_width = get_bin_limits(alignments+i, total_pairs-i);
        if(bin_width > 1) {
            ndups += mark_bin(alignments, i, bin_width, bitmap);
        }
        i+=bin_width;
    }
    return ndups;
}

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] input.bam output.bam\n", prog);
    printf("\n\
This program will parse a BAM file produced by bison and mark likely PCR\n\
duplicates in a new file. A PCR duplicate is defined as two reads/read pairs\n\
having the same start coordinates on the same strand of the same\n\
chromosome/contig and having a methylation call edit-distance of 0 (ignoring\n\
instances of a no-call in one read and a call in the other). The read or pair\n\
with the greater sum of base-call phred scores is left unmarked. Note that\n\
soft-clipped bases are included in computing the 5' coordinate of the read,\n\
as is also the case with picard markDuplicates.\n\
\n\
N.B., only methylation calls with a phred score of at least 5 are compared.\n\
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
    -@ INT Number of compression threads to use when writing the BAM file.\n\
\n\
    -g INT How much the array should grow (measured in reads/read pairs) when\n\
           it's full). The default is 1000000.\n");
}

inline int32_t get_strand(bam1_t *read) {
    char *XG = bam_aux2Z(bam_aux_get(read, "XG"));
    char *XR = bam_aux2Z(bam_aux_get(read, "XR"));

    if(read->core.flag & BAM_FREAD1) {
        if(*XG == 'C') { //OT or CTOT
            if(*XR == 'C') return 0; //OT
            else return 2; //CTOT
        } else {
            if(*XR == 'C') return 1; //OB
            else return 3; //CTOB
        }
    } else { //For read #2, things are a bit different
        if(*XG == 'C') { //OT or CTOT
            if(*XR == 'C') return 2; //CTOT
            else return 0; //OT
        } else {
            if(*XR == 'C') return 3; //CTOB
            else return 1; //OB
        }
    }
}

inline int32_t get_pos(bam1_t *read) {
    int32_t pos = read->core.pos;
    uint32_t *CIGAR = bam_get_cigar(read);
    int32_t op, nops;

    if(read->core.n_cigar == 1) return pos;

    op = CIGAR[0] & 0xf;
    nops = CIGAR[0] >> 4;
    if(op == BAM_CHARD_CLIP) { //The only thing that can precede an S is an H
        op = CIGAR[1] & 0xf;
        nops = CIGAR[1] >> 4;
        if(op == BAM_CSOFT_CLIP) pos -= nops;
    } else if(op == BAM_CSOFT_CLIP) {
        pos -= nops;
    }
    if(pos<0) pos=0;
    return pos;
}

//return the tid, such that first 2 bits denote the strand
inline int32_t get_tid(bam1_t *read) {
    int32_t strand = get_strand(read);

    strand = strand << 30;
    assert(read->core.tid < 0x3FFFFFFF);
    return strand | read->core.tid;
}

//Given an XM string, set 4 64bit bit-fields denoting methylated/unmethylated positions
//This is only sufficient for reads up to 256 bases
//Only include things if a call has a phred score of at least MIN_PHRED
//Soft-clipping will offset things
void calls2uint(bam1_t *read, uint64_t *meth, uint64_t *unmeth) {
    int32_t i, l = read->core.l_qseq, start, end, offset = get_pos(read)-read->core.pos;
    uint8_t *p = bam_aux_get(read, "XM");
    char *XM = NULL;
    uint8_t *qual = bam_get_qual(read);
    int word_offset = 0;

    if(read->core.flag & BAM_FPAIRED) {
        word_offset = (read->core.flag & BAM_FREAD1) ? 0 : 4;
    }

    *meth = 0;
    *unmeth = 0;
    if(p != NULL) {
        XM = bam_aux2Z(p);
        //chunk #1
        start = 0;
        end = min(l-offset, 64-offset);
        for(i=start; i<end; i++) {
            switch(*(XM+i)) {
                case 'Z' : case 'X' : case 'H' :
                    if(*(qual+i) >= MIN_PHRED) meth[0+word_offset] |= (1 << BIT_OFFSET(i+offset)); //This is just set_bit() with a 64bit map
                    break;
                case 'z' : case 'x' : case 'h' :
                    if(*(qual+i) >= MIN_PHRED) unmeth[0+word_offset] |= (1 << BIT_OFFSET(i+offset));
                    break;
            }
        }
        if(l+offset > 64) {
            //chunk #2
            start = 64-offset;
            end = min(l-offset, 128-offset);
            for(i=start; i<end; i++) {
                switch(*(XM+i)) {
                    case 'Z' : case 'X' : case 'H' :
                        if(*(qual+i) >= MIN_PHRED) meth[1+word_offset] |= (1 << BIT_OFFSET(i+offset)); //This is just set_bit() with a 64bit map
                        break;
                    case 'z' : case 'x' : case 'h' :
                        if(*(qual+i) >= MIN_PHRED) unmeth[1+word_offset] |= (1 << BIT_OFFSET(i+offset));
                        break;
                }
            }
        }
        if(l+offset > 128) {
            //chunk #3
            start = 128-offset;
            end = min(l-offset, 192-offset);
            for(i=start; i<end; i++) {
                switch(*(XM+i)) {
                    case 'Z' : case 'X' : case 'H' :
                        if(*(qual+i) >= MIN_PHRED) meth[2+word_offset] |= (1 << BIT_OFFSET(i+offset)); //This is just set_bit() with a 64bit map
                        break;
                    case 'z' : case 'x' : case 'h' :
                        if(*(qual+i) >= MIN_PHRED) unmeth[2+word_offset] |= (1 << BIT_OFFSET(i+offset));
                        break;
                }
            }
        }
        if(l+offset > 192) {
            //chunk #2
            start = 192-offset;
            end = min(l-offset, 256-offset);
            for(i=start; i<end; i++) {
                switch(*(XM+i)) {
                    case 'Z' : case 'X' : case 'H' :
                        if(*(qual+i) >= MIN_PHRED) meth[3+word_offset] |= (1 << BIT_OFFSET(i+offset)); //This is just set_bit() with a 64bit map
                        break;
                    case 'z' : case 'x' : case 'h' :
                        if(*(qual+i) >= MIN_PHRED) unmeth[3+word_offset] |= (1 << BIT_OFFSET(i+offset));
                        break;
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    bam1_t *read;
    bam_hdr_t *header;
    htsFile *fp = NULL, *of = NULL;
    uint64_t total_pairs = 0, max_length = 10000000;
    uint64_t bitmap_length = 0, cur_read = 0, ndups = 0;
    uint64_t grow_size = 1000000;
    alignment *alignments;
    uint32_t *bitmap;
    int i, n_compression_threads = 1;
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
        } else if(strcmp("-@", argv[i]) == 0) {
            n_compression_threads = atoi(argv[++i]);
        } else if(iname == NULL) {
            iname = argv[i];
            fp = sam_open(iname, "rb");
        } else if(oname == NULL) {
            oname = argv[i];
        } else {
            fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }
    if(iname == NULL || oname == NULL) {
        fprintf(stderr, "Either the input or output file name were not specified!\n");
        usage(argv[0]);
        return 2;
    }

    //Set everything up
    header = sam_hdr_read(fp);
    read = bam_init1();
    alignments = malloc(sizeof(alignment)*max_length);
    assert(alignments != NULL);

    //Initialize lookup table
    fill_table();

    //Read in the alignments
    while(sam_read1(fp, header, read) > 1) {
        if(read->core.flag & BAM_FUNMAP) continue; //We just skip over unmapped reads

        alignments[total_pairs].tid1 = get_tid(read);
        alignments[total_pairs].start1 = get_pos(read);
        alignments[total_pairs].total_phred = total_phred(read);
        for(i=0;i<8;i++) { //Ensure that everything starts at 0!
            alignments[total_pairs].meth[i] = 0;
            alignments[total_pairs].unmeth[i] = 0;
        }
        alignments[total_pairs].start2 = 0;
        calls2uint(read, alignments[total_pairs].meth, alignments[total_pairs].unmeth);
        if(read->core.flag & BAM_FPAIRED && !(read->core.flag & BAM_FMUNMAP)) {
            assert(sam_read1(fp,header,read)>1);
            alignments[total_pairs].tid2 = get_tid(read);
            alignments[total_pairs].start2 = get_pos(read);
            alignments[total_pairs].total_phred += total_phred(read);
            calls2uint(read, alignments[total_pairs].meth, alignments[total_pairs].unmeth);
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
    sam_close(fp);

    //create bitmap
    bitmap_length = max_length/32;
    bitmap_length += (max_length % 32 > 0) ? 1 : 0;
    bitmap = calloc(bitmap_length, sizeof(uint32_t));

    //Sort
    qsort((void *) alignments, (size_t) total_pairs, sizeof(alignment), &comp_func);

    //Mark duplicates in bitmap
    ndups = mark_dups(alignments, bitmap, total_pairs);
    free(alignments);
    fprintf(stderr, "There were %"PRIu64" duplicates from %"PRIu64" total reads or pairs\n", ndups, total_pairs);

    //reopen file, iterate through and change flags as appropriate
    fp = sam_open(iname, "rb");
    header = sam_hdr_read(fp);
    of = sam_open(oname, "wb");
    if(n_compression_threads > 1) hts_set_threads(of, n_compression_threads);
    sam_hdr_write(of, header);

    while(sam_read1(fp, header, read) > 1) {
        //Only set as duplicate if it's mapped
        if(!(read->core.flag & BAM_FUNMAP)) if(get_bit(bitmap, cur_read)) read->core.flag = read->core.flag | BAM_FDUP;
        sam_write1(of, header, read);
        if(read->core.flag & BAM_FPAIRED) {
            assert(sam_read1(fp, header, read) > 1);
            if(!(read->core.flag & BAM_FUNMAP)) if(get_bit(bitmap, cur_read)) read->core.flag = read->core.flag | BAM_FDUP;
            sam_write1(of, header, read);
        }
        if(!(read->core.flag & BAM_FUNMAP)) cur_read++;
    }

    //Clean up
    sam_close(fp);
    sam_close(of);
    bam_destroy1(read);
    free(bitmap);

    return 0;
}
