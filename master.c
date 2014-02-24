#include "bison.h"
#include <math.h>
#include <sys/time.h>

typedef struct {
    unsigned long long t_reads; //total reads
    unsigned long long m_reads_OT; //reads mapped to the OT
    unsigned long long m_reads_OB;
    unsigned long long m_reads_CTOT;
    unsigned long long m_reads_CTOB;
    unsigned long long t_CpG; //Total CpGs
    unsigned long long m_CpG; //Methylated CpGs
    unsigned long long t_CHG;
    unsigned long long m_CHG;
    unsigned long long t_CHH;
    unsigned long long m_CHH;
} metrics_struct;

/******************************************************************************
*
*   Update the CpG/CHG/CHH metrics according to the methylation calls in a read
*
*******************************************************************************/
void update_counts(bam1_t *read, metrics_struct *metrics) {
    char *XM = bam_aux2Z(bam_aux_get(read, "XM"));
    char base;
    int i;

    for(i=0; i<read->core.l_qseq; i++) {
        base = *(XM+i);
        if(base != '.') {
            if(base == 'Z') {
                metrics->t_CpG++;
                metrics->m_CpG++;
            } else if(base == 'z') {
                metrics->t_CpG++;
            } else if(base == 'X') {
                metrics->t_CHG++;
                metrics->m_CHG++;
            } else if(base == 'x') {
                metrics->t_CHG++;
            } else if(base == 'H') {
                metrics->t_CHH++;
                metrics->m_CHH++;
            } else if(base == 'h') {
                metrics->t_CHH++;
            }
        }
    }
}

/******************************************************************************
*
*   Return the alignment score or -MAX_INT if unaligned
*
*   bam1_t *read: the read in question
*
*******************************************************************************/
int get_AS(bam1_t *read) {
    int AS = INT_MIN>>2;
    uint8_t *p = bam_aux_get(read, "AS");

    if(read->core.flag & BAM_FUNMAP) return AS;
    if(p != NULL) AS = bam_aux2i(p);
    return AS;
}

/******************************************************************************
*
*   Calculate the minimum score for a given readlength
*
*   int32_t rlen: a read length
*
*******************************************************************************/
inline int scoreMin(int32_t rlen) {
    //Return different values, depending on --score-min
    if(config.scoremin_type == 'L') {
        return (config.scoremin_intercept + config.scoremin_coef * rlen);
    } else if(config.scoremin_type == 'S') {
        return (config.scoremin_intercept + config.scoremin_coef * sqrt((float) rlen));
    } else if(config.scoremin_type == 'G') {
        return (config.scoremin_intercept + config.scoremin_coef * log((float) rlen));
    } else { //'C'
        return (config.scoremin_intercept + config.scoremin_coef);
    }
}

/******************************************************************************
*
*   Return the secondary alignment score or -MAX_INT if unaligned
*
*   bam1_t *read: the read in question
*
*******************************************************************************/
int get_XS(bam1_t *read) {
    int XS = INT_MIN>>2;
    uint8_t *p = bam_aux_get(read, "XS");

    if(read->core.flag & BAM_FUNMAP) return XS;
    if(p != NULL) XS = bam_aux2i(p);
    return XS;
}

/******************************************************************************
*
*   Calculate a MAPQ, given AS, XS, and the minimum score (ala bowtie2)
*
*******************************************************************************/
int calc_MAPQ_BT2(int AS, int XS, int scMin) {
    int diff, bestOver, bestdiff;
    diff = abs(scMin); //Range of possible alignment scores
    bestOver = AS-scMin; //Shift alignment score range, so worst score is 0
    
    //This seems like an odd way to calculate this!

    //The method depends on config.mode
    bestdiff = (int) abs(abs((float) AS)-abs((float) XS)); //Absolute distance between alignment scores
    if(config.mode == 0) { //--end-to-end (default)
        if(XS < scMin) {
            if(bestOver >= diff * (double) 0.8f) return 42;
            else if(bestOver >= diff * (double) 0.7f) return 40;
            else if(bestOver >= diff * (double) 0.6f) return 24;
            else if(bestOver >= diff * (double) 0.5f) return 23;
            else if(bestOver >= diff * (double) 0.4f) return 8;
            else if(bestOver >= diff * (double) 0.3f) return 3;
            else return 0;
        } else {
            if(bestdiff >= diff * (double) 0.9f) {
                if(bestOver == diff) {
                    return 39;
                } else {
                    return 33;
                }
            } else if(bestdiff >= diff * (double) 0.8f) {
                if(bestOver == diff) {
                    return 38;
                } else {
                    return 27;
                }
            } else if(bestdiff >= diff * (double) 0.7f) {
                if(bestOver == diff) {
                    return 37;
                } else {
                    return 26;
                }
            } else if(bestdiff >= diff * (double) 0.6f) {
                if(bestOver == diff) {
                    return 36;
                } else {
                    return 22;
                }
            } else if(bestdiff >= diff * (double) 0.5f) {
                if(bestOver == diff) {
                    return 35;
                } else if(bestOver >= diff * (double) 0.84f) {
                    return 25;
                } else if(bestOver >= diff * (double) 0.68f) {
                    return 16;
                } else {
                    return 5;
                }
            } else if(bestdiff >= diff * (double) 0.4f) {
                if(bestOver == diff) {
                    return 34;
                } else if(bestOver >= diff * (double) 0.84f) {
                    return 21;
                } else if(bestOver >= diff * (double) 0.68f) {
                    return 14;
                } else {
                    return 4;
                }
            } else if(bestdiff >= diff * (double) 0.3f) {
                if(bestOver == diff) {
                    return 32;
                } else if(bestOver >= diff * (double) 0.88f) {
                    return 18;
                } else if(bestOver >= diff * (double) 0.67f) {
                    return 15;
                } else {
                    return 3;
                }
            } else if(bestdiff >= diff * (double) 0.2f) {
                if(bestOver == diff) {
                    return 31;
                } else if(bestOver >= diff * (double) 0.88f) {
                    return 17;
                } else if(bestOver >= diff * (double) 0.67f) {
                    return 11;
                } else {
                    return 0;
                }
            } else if(bestdiff >= diff * (double) 0.1f) {
                if(bestOver == diff) {
                    return 30;
                } else if(bestOver >= diff * (double) 0.88f) {
                    return 12;
                } else if(bestOver >= diff * (double) 0.67f) {
                    return 7;
                } else {
                    return 0;
                }
            } else if(bestdiff > 0) {
                if(bestOver >= diff * (double)0.67f) {
                    return 6;
                } else {
                    return 2;
                }
            } else {
                if(bestOver >= diff * (double)0.67f) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    } else { //--local
        if(XS < scMin) {
            if(bestOver >= diff * (double) 0.8f) return 44;
            else if(bestOver >= diff * (double) 0.7f) return 42;
            else if(bestOver >= diff * (double) 0.6f) return 41;
            else if(bestOver >= diff * (double) 0.5f) return 36;
            else if(bestOver >= diff * (double) 0.4f) return 28;
            else if(bestOver >= diff * (double) 0.3f) return 24;
            else return 22;
        } else {
            if(bestdiff >= diff * (double) 0.9f) return 40;
            else if(bestdiff >= diff * (double) 0.8f) return 39;
            else if(bestdiff >= diff * (double) 0.7f) return 38;
            else if(bestdiff >= diff * (double) 0.6f) return 37;
            else if(bestdiff >= diff * (double) 0.5f) {
                if     (bestOver == diff)       return 35;
                else if(bestOver >= diff * (double) 0.5f) return 25;
                else                            return 20;
            } else if(bestdiff >= diff * (double) 0.4f) {
                if     (bestOver == diff)       return 34;
                else if(bestOver >= diff * (double) 0.5f) return 21;
                else                            return 19;
            } else if(bestdiff >= diff * (double) 0.3f) {
                if     (bestOver == diff)       return 33;
                else if(bestOver >= diff * (double) 0.5f) return 18;
                else                            return 16;
            } else if(bestdiff >= diff * (double) 0.2f) {
                if     (bestOver == diff)       return 32;
                else if(bestOver >= diff * (double) 0.5f) return 17;
                else                            return 12;
            } else if(bestdiff >= diff * (double) 0.1f) {
                if     (bestOver == diff)       return 31;
                else if(bestOver >= diff * (double) 0.5f) return 14;
                else                            return 9;
            } else if(bestdiff > 0) {
                if(bestOver >= diff * (double) 0.5f)      return 11;
                else                            return 2;
            } else {
                if(bestOver >= diff * (double) 0.5f)      return 1;
                else                            return 0;
            }
        }
    }
}

/******************************************************************************
*
*   Determine whether the alignment is actually unique by comparing the AS and
*   XS auxiliary tags.
*
*   bam1_t *read: The read to look at
*
*******************************************************************************/
int unique_alignment(bam1_t *read) {
    int AS, XS;

    AS = bam_aux2i(bam_aux_get(read, "AS"));
    if(bam_aux_get(read, "XS") == 0) return 1;
    XS = bam_aux2i(bam_aux_get(read, "XS"));
    if(AS > XS) return 1;
    return 0;
}

/******************************************************************************
*
*   Replace the stored sequence in a read.
*
*   bam1_t *read: The read whose sequence will be replaced
*   char *seq: Sequence to coopy into read.
*
*   If read is reverse complemented, the same will be done to seq.
*
*******************************************************************************/
void swap_sequence(bam1_t *read, char *seq) {
    uint8_t *sequence = bam1_seq(read), val;
    char *seq2 = strdup(seq);
    int i, j;

    //Do we need to reverse complement?
    if(read->core.flag & BAM_FREVERSE) reverse_complement(seq2);
    for(i=0, j=0; i<strlen(seq2); i+=2, j++) {
        if(*(seq2+i) == 'A') {
            if(*(seq2+i+1) == 'A') val = 17;
            else if(*(seq2+i+1) == 'C') val = 18;
            else if(*(seq2+i+1) == 'G') val = 20;
            else if(*(seq2+i+1) == 'T') val = 24;
            else if(*(seq2+i+1) == 'N') val = 31;
            else val = 16;
        } else if(*(seq2+i) == 'C') {
            if(*(seq2+i+1) == 'A') val = 33;
            else if(*(seq2+i+1) == 'C') val = 34;
            else if(*(seq2+i+1) == 'G') val = 36;
            else if(*(seq2+i+1) == 'T') val = 40;
            else if(*(seq2+i+1) == 'N') val = 47;
            else val = 32;
        } else if(*(seq2+i) == 'G') {
            if(*(seq2+i+1) == 'A') val = 65;
            else if(*(seq2+i+1) == 'C') val = 66;
            else if(*(seq2+i+1) == 'G') val = 68;
            else if(*(seq2+i+1) == 'T') val = 72;
            else if(*(seq2+i+1) == 'N') val = 79;
            else val = 64;
        } else if(*(seq2+i) == 'T') {
            if(*(seq2+i+1) == 'A') val = 129;
            else if(*(seq2+i+1) == 'C') val = 130;
            else if(*(seq2+i+1) == 'G') val = 132;
            else if(*(seq2+i+1) == 'T') val = 136;
            else if(*(seq2+i+1) == 'N') val = 143;
            else val = 128;
        } else {
            if(*(seq2+i+1) == 'A') val = 241;
            else if(*(seq2+i+1) == 'C') val = 242;
            else if(*(seq2+i+1) == 'G') val = 244;
            else if(*(seq2+i+1) == 'T') val = 248;
            else if(*(seq2+i+1) == 'N') val = 255;
            else val = 240;
        }
        *(sequence+j) = val;
    }
    free(seq2);
}

/******************************************************************************
*
*   Return the XM string that will be appended to a read.
*
*   bam1_t *read; the read in question
*   char *XG: The XG tag, indicating which coversion to pay attention to.
*
*   THE OUTPUT MUST BE free()d
*******************************************************************************/
char *callXM(bam1_t *read, char *XG) {
    char *chrom = lookup_chrom(read);
    unsigned long long offset = genome_offset(chrom, 0), current_position;
    unsigned long long chrom_end = genome_chrom_length(chrom);
    unsigned long long *genomic_position = calculate_positions(read);

    char *read_seq = calloc(1+read->core.l_qseq, sizeof(char));
    char *XM = calloc(1+read->core.l_qseq, sizeof(char));
    char genome_base, read_base, *bases;
    int i;
    uint8_t b;

    //Extract the read sequence
    for(i=0; i<read->core.l_qseq; i++) {
        b = bam1_seqi(bam1_seq(read), i);
        if(b == 1) {
            *(read_seq+i) = 'A';
        } else if(b == 2) {
            *(read_seq+i) = 'C';
        } else if(b == 4) {
            *(read_seq+i) = 'G';
        } else if(b == 8) {
            *(read_seq+i) = 'T';
        } else if(b == 15) {
            *(read_seq+i) = 'N';
        }
        current_position = *(genomic_position+i);
    }

    for(i=0; i<read->core.l_qseq; i++) {
        current_position = *(genomic_position+i);
        if(current_position == ULLONG_MAX) {
            *(XM+i) = '.';
            continue;
        }
        genome_base = toupper(*(chromosomes.genome+offset+current_position));
        read_base = toupper(*(read_seq+i));
        if(read_base != genome_base) {
            //Mismatches to the top and bottom strands are treated differently
            if(*XG == 'C') { //OT or CTOT
                if(genome_base == 'C' && read_base == 'T') {
                    bases = get_genomic_context(offset, current_position, 2, chrom_end);
                    if(*(bases+1) == 'G') {
                        //Unmethylated CpG
                        *(XM+i) = 'z';
                    } else if(*(bases+2) == 'G') {
                        //Unmethylated CHG
                        *(XM+i) = 'x';
                    } else {
                        //Unmethylated CHH
                        *(XM+i) = 'h';
                    }
                    free(bases);
                } else {
                    //Just a mismatch
                    *(XM+i) = '.';
                }
            } else { //OB or CTOB
                if(genome_base == 'G' && read_base == 'A') {
                    bases = get_genomic_context(offset, current_position, -2, chrom_end);
                    if(*(bases+1) == 'C') {
                        //Unmethylated CpG
                        *(XM+i) = 'z';
                    } else if(*(bases+0) == 'C') {
                        //Unmethylated CHG
                        *(XM+i) = 'x';
                    } else {
                        //Unmethylated CHH
                        *(XM+i) = 'h';
                    }
                    free(bases);
                } else {
                    *(XM+i) = '.';
                }
            }
        } else {
            if(*XG == 'C') { //OT or CTOT
                if(genome_base == 'C') {
                    bases = get_genomic_context(offset, current_position, 2, chrom_end);
                    if(*(bases+1) == 'G') {
                        //Methylated CpG
                        *(XM+i) = 'Z';
                    } else if(*(bases+2) == 'G') {
                        //Methylated CHG
                        *(XM+i) = 'X';
                    } else {
                        //Methylated CHH
                        *(XM+i) = 'H';
                    }
                    free(bases);
                } else {
                    *(XM+i) = '.';
                }
            } else { //OB or CTOB
                if(genome_base == 'G') {
                    bases = get_genomic_context(offset, current_position, -2, chrom_end);
                    if(*(bases+1) == 'C') {
                        //Methylated CpG
                        *(XM+i) = 'Z';
                    } else if(*(bases+0) == 'C') {
                        //Methylated CHG
                        *(XM+i) = 'X';
                    } else {
                        //Methylated CHH
                        *(XM+i) = 'H';
                    }
                    free(bases);
                } else {
                    *(XM+i) = '.';
                }
            }
        }
    }

    free(read_seq);
    free(genomic_position);
    return XM;
}

/******************************************************************************
*
*   As with callXM, but return the mismatches with the reference.
*
*   bam_t *read; the read in question
*   char *XM: output from callXM
*   char *XG: The XG tag, indicating which coversion to pay attention to.
*
*   THE OUTPUT MUST BE fre()d
*   The length of XX is currently limited to MAXREAD!!!
*
*******************************************************************************/
char *callXX(bam1_t *read, char *XM, char *XG) {
    char *chrom = lookup_chrom(read);
    unsigned long long offset = genome_offset(chrom, 0), current_position;
    unsigned long long *genomic_position = calculate_positions(read);
    uint8_t base, NM = 0;

    char *read_seq = calloc(1+read->core.l_qseq, sizeof(char));
    char *XX = calloc(MAXREAD, sizeof(char));
    int i, good = 0;

    //Extract the read sequence
    for(i=0; i<read->core.l_qseq; i++) {
        base = bam1_seqi(bam1_seq(read), i);
        if(base == 1) *(read_seq+i) = 'A';
        else if(base == 2) *(read_seq+i) = 'C';
        else if(base == 4) *(read_seq+i) = 'G';
        else if(base == 8) *(read_seq+i) = 'T';
        else *(read_seq+i) = 'N';
        current_position = *(genomic_position+i);
    }

    //Create the XM string
    for(i=0; i<strlen(XM); i++) {
        if(*(XM+i) != '.') {
            //unlike bismark, we don't count methylation changes as mismatches
            good++;
        } else {
            current_position = *(genomic_position+i);
            base = toupper(*(chromosomes.genome+offset+current_position));
            if(base != *(read_seq+i)) {
                NM++;
                //If the read starts with a mismatch, the XX string should start with a 0
                if(i == 0) {
                    sprintf(XX, "0%c", *(read_seq+i));
                } else if(good) {
                    sprintf(XX, "%s%i%c", XX, good, *(read_seq+i));
                    good = 0;
                } else {
                    sprintf(XX, "%s%c", XX, *(read_seq+i));
                }
            } else {
                good++;
            }
        }
    }
    if(good) sprintf(XX, "%s%i", XX, good);

    //Update the NM tag
    *(bam_aux_get(read, "NM")+1) = NM;

    free(read_seq);
    free(genomic_position);

    return XX;
}

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
int process_single(bam1_t *read1, bam1_t *read2, bam1_t *read3, bam1_t *read4, char *seq) {
    int AS1=0, AS2=0, AS3=0, AS4=0;
    bam1_t *tmp_read = NULL;
    char *XM, *XX, XG[] = "CT", XR[] = "CT";
    int best_node = 0;
    kstring_t *kXM = (kstring_t *) calloc(1, sizeof(kstring_t));
    kstring_t *kXX = (kstring_t *) calloc(1, sizeof(kstring_t));
    //For recalculating the MAPQ ala bowtie2 v2 MAPQ calculator
    int XS, scMin, MAPQ = 0, AS=0;

    //Determine the read with the highest alignment score
    AS1 = get_AS(read1);
    AS2 = get_AS(read2);
    if(!config.directional) {
        AS3 = get_AS(read3);
        AS4 = get_AS(read4);
    }
    if(config.directional) {
        if(AS1 > AS2) {
            sprintf(XR, "CT");
            sprintf(XG, "CT");
            if(!(read1->core.flag & BAM_FUNMAP)) {
                tmp_read = read1;
                best_node = 1;
            }
        } else if(AS2 > AS1) {
            sprintf(XR, "CT");
            sprintf(XG, "GA");
            if(!(read2->core.flag & BAM_FUNMAP)) {
                tmp_read = read2;
                best_node = 2;
            }
        }
    } else {
        if(AS1 > AS2 && AS1 > AS3 && AS1 > AS4) { //OT
            sprintf(XR, "CT");
            sprintf(XG, "CT");
            if(!(read1->core.flag & BAM_FUNMAP)) {
                tmp_read = read1;
                best_node = 1;
            }
        } else if(AS2 > AS1 && AS2 > AS3 && AS2 > AS4) { //OB
            sprintf(XR, "CT");
            sprintf(XG, "GA");
            if(!(read2->core.flag & BAM_FUNMAP)) {
                tmp_read = read2;
                best_node = 2;
            }
        } else if(AS3 > AS1 && AS3 > AS2 && AS3 > AS4) { //CTOT
            sprintf(XR, "GA");
            sprintf(XG, "CT");
            if(!(read3->core.flag & BAM_FUNMAP)) {
                tmp_read = read3;
                best_node = 3;
            }
        } else if(AS4 > AS1 && AS4 > AS2 && AS4 > AS3) { //CTOB
            sprintf(XR, "GA");
            sprintf(XG, "GA");
            if(!(read4->core.flag & BAM_FUNMAP)) {
                tmp_read = read4;
                best_node = 4;
            }
        }
    }

    //If there is no best score (tmp_read == NULL), mark read1 as unmapped
    if(tmp_read == NULL) {
        swap_sequence(read1, seq);
        read1->core.flag = read1->core.flag | 0x4;
        best_node = 1;
    } else {
        swap_sequence(tmp_read, seq);
        XM = callXM(tmp_read, XG);
        XX = callXX(tmp_read, XM, XG);
        //append the tags
        kputs(XX, kXX);
        kputs(XM, kXM);

        bam_aux_del(tmp_read, bam_aux_get(tmp_read, "XM"));
        bam_aux_del(tmp_read, bam_aux_get(tmp_read, "XG"));
        bam_aux_append(tmp_read, "XX", 'Z', kXX->l + 1, (uint8_t*) kXX->s);
        bam_aux_append(tmp_read, "XM", 'Z', kXM->l + 1, (uint8_t*) kXM->s);
        bam_aux_append(tmp_read, "XR", 'Z', 3, (uint8_t*) XR);
        bam_aux_append(tmp_read, "XG", 'Z', 3, (uint8_t*) XG);
        free(kXX->s);
        free(kXM->s);
        free(XM);
        free(XX);

        //Recalculate MAPQ and replace the XS score
        scMin = scoreMin(tmp_read->core.l_qseq);
        XS = get_XS(tmp_read);
        if(best_node == 1) {
            AS = AS1;
            if(AS2 > XS) XS = AS2;
            if(!config.directional) {
                if(AS3 > XS) XS = AS3;
                if(AS4 > XS) XS = AS4;
            }
        }
        if(best_node == 2) {
            AS = AS2;
            if(AS1 > XS) XS = get_AS(read2);
            if(!config.directional) {
                if(AS3 > XS) XS = AS3;
                if(AS4 > XS) XS = AS4;
            }
        }
        if(best_node == 3) {
            AS = AS3;
            if(AS1 > XS) XS = AS1;
            if(AS2 > XS) XS = AS2;
            if(AS4 > XS) XS = AS4;
        }
        if(best_node == 4) {
            AS = AS4;
            if(AS1 > XS) XS = AS1;
            if(AS2 > XS) XS = AS2;
            if(AS3 > XS) XS = AS3;
        }
        MAPQ = calc_MAPQ_BT2(AS, XS, scMin);
        MAPQ = (MAPQ < tmp_read->core.qual) ? MAPQ : tmp_read->core.qual;
        tmp_read->core.qual = MAPQ;
        if(XS >= scMin) {
            //replace/add the XS tag
            if(bam_aux_get(tmp_read, "XS")) bam_aux_del(tmp_read, bam_aux_get(tmp_read, "XS"));
            bam_aux_append(tmp_read, "XS", 'i', 4, (uint8_t*) &XS);
        }
    }
    free(kXX);
    free(kXM);
    return best_node;
}

/******************************************************************************
*
*   Like process_single, but for paired_end reads. The bam1_t**s hold the
*   buffered reads. i denotes the read#1 of interest (read #2 is the next read)
*
*******************************************************************************/
int process_paired(bam1_t **read1, bam1_t **read2, bam1_t **read3, bam1_t **read4, char **seq) {
    int AS1=0, AS2=0, AS3=0, AS4=0;
    bam1_t *tmp_read1 = NULL, *tmp_read2 = NULL;
    char *XM1, *XM2, *XX1, *XX2, XG[] = "CT", XR1[] = "CT", XR2[] = "CT";
    kstring_t *kXM1 = (kstring_t *) calloc(1, sizeof(kstring_t));
    kstring_t *kXM2 = (kstring_t *) calloc(1, sizeof(kstring_t));
    kstring_t *kXX1 = (kstring_t *) calloc(1, sizeof(kstring_t));
    kstring_t *kXX2 = (kstring_t *) calloc(1, sizeof(kstring_t));
    int best_node = 0;
    //For MAPQ/XS replacement
    int MAPQ, XS1, XS2, scMin1, scMin2;

    //Determine the read with the highest alignment score
    AS1 = get_AS(*(read1)) + get_AS(*(read1+1));
    AS2 = get_AS(*(read2)) + get_AS(*(read2+1));
    if(!config.directional) {
        AS3 = get_AS(*(read3)) + get_AS(*(read3+1));
        AS4 = get_AS(*(read4)) + get_AS(*(read4+1));
    }
    if(config.directional) {
        if(AS1 > AS2) { //OT
            sprintf(XR1, "CT");
            sprintf(XR2, "GA");
            sprintf(XG, "CT");
            if(!((*(read1))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read1);
                tmp_read2 = *(read1+1);
                best_node = 1;
            }
        } else if(AS2 > AS1) { //OB
            sprintf(XR1, "CT");
            sprintf(XR2, "GA");
            sprintf(XG, "GA");
            if(!((*(read2))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read2);
                tmp_read2 = *(read2+1);
                best_node = 2;
            }
        }
    } else {
        if(AS1 > AS2 && AS1 > AS3 && AS1 > AS4) { //OT
            sprintf(XR1, "CT");
            sprintf(XR2, "GA");
            sprintf(XG, "CT");
            if(!((*(read1))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read1);
                tmp_read2 = *(read1+1);
                best_node = 1;
            }
        } else if(AS2 > AS1 && AS2 > AS3 && AS2 > AS4) { //OB
            sprintf(XR1, "CT");
            sprintf(XR2, "GA");
            sprintf(XG, "GA");
            if(!((*(read2))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read2);
                tmp_read2 = *(read2+1);
                best_node = 2;
            }
        } else if(AS3 > AS1 && AS3 > AS2 && AS3 > AS4) { //CTOT
            sprintf(XR1, "GA");
            sprintf(XR2, "CT");
            sprintf(XG, "CT");
            if(!((*(read3))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read3);
                tmp_read2 = *(read3+1);
                best_node = 3;
            }
        } else if(AS4 > AS1 && AS4 > AS2 && AS4 > AS3) { //CTOB
            sprintf(XR1, "GA");
            sprintf(XR2, "CT");
            sprintf(XG, "GA");
            if(!((*(read4))->core.flag & BAM_FUNMAP)) {
                tmp_read1 = *(read4);
                tmp_read2 = *(read4+1);
                best_node = 4;
            }
        }
    }

    //If there is no best score (tmp_read == NULL), mark read1 as unmapped
    if(tmp_read1 == NULL) {
        swap_sequence(*(read1), *(seq));
        swap_sequence(*(read1+1), *(seq+1));
        (*(read1))->core.flag = (*(read1))->core.flag | 0x4;
        (*(read1+1))->core.flag = (*(read1+1))->core.flag | 0x4;
        best_node = 1;
    } else {
        swap_sequence(tmp_read1, *(seq));
        swap_sequence(tmp_read2, *(seq+1));
        XM1 = callXM(tmp_read1, XG);
        XX1 = callXX(tmp_read1, XM1, XG);
        XM2 = callXM(tmp_read2, XG);
        XX2 = callXX(tmp_read2, XM2, XG);

        kputs(XX1, kXX1);
        kputs(XX2, kXX2);
        kputs(XM1, kXM1);
        kputs(XM2, kXM2);

        bam_aux_del(tmp_read1, bam_aux_get(tmp_read1, "XM"));
        bam_aux_del(tmp_read2, bam_aux_get(tmp_read2, "XM"));
        bam_aux_del(tmp_read1, bam_aux_get(tmp_read1, "XG"));
        bam_aux_del(tmp_read2, bam_aux_get(tmp_read2, "XG"));

        bam_aux_append(tmp_read1, "XX", 'Z', kXX1->l + 1, (uint8_t*) kXX1->s);
        bam_aux_append(tmp_read2, "XX", 'Z', kXX2->l + 1, (uint8_t*) kXX2->s);

        bam_aux_append(tmp_read1, "XM", 'Z', kXM1->l + 1, (uint8_t*) kXM1->s);
        bam_aux_append(tmp_read2, "XM", 'Z', kXM2->l + 1, (uint8_t*) kXM2->s);

        bam_aux_append(tmp_read1, "XR", 'Z', 3, (uint8_t*) XR1);
        bam_aux_append(tmp_read2, "XR", 'Z', 3, (uint8_t*) XR2);

        bam_aux_append(tmp_read1, "XG", 'Z', 3, (uint8_t*) XG);
        bam_aux_append(tmp_read2, "XG", 'Z', 3, (uint8_t*) XG);
        free(kXX1->s);
        free(kXX2->s);
        free(kXM1->s);
        free(kXM2->s);
        free(XM1);
        free(XM2);
        free(XX1);
        free(XX2);

        //Recalculate MAPQ and replace the XS score
        scMin1 = scoreMin(tmp_read1->core.l_qseq);
        scMin2 = scoreMin(tmp_read2->core.l_qseq);
        XS1 = get_XS(tmp_read1);
        XS2 = get_XS(tmp_read2);
        if(XS2 < scMin2) {
            if(XS1 >= scMin1) XS2 = get_AS(tmp_read2);
        } else if(XS1 < scMin1) {
            if(XS2 >= scMin2) XS1 = get_AS(tmp_read1);
        }
        if(best_node == 1) {
            if(AS2 > XS1+XS2) {
                XS1 = get_AS(*(read2)); 
                XS2 = get_AS(*(read2+1)); 
            }
            if(!config.directional) {
                if(AS3 > XS1+XS2) {
                    XS1 = get_AS(*(read3)); 
                    XS2 = get_AS(*(read3+1)); 
                }
                if(AS4 > XS1+XS2) {
                    XS1 = get_AS(*(read4)); 
                    XS2 = get_AS(*(read4+1)); 
                }
            }
        }
        if(best_node == 2) {
            if(AS1 > XS1+XS2) {
                XS1 = get_AS(*(read1)); 
                XS2 = get_AS(*(read1+1)); 
            }
            if(!config.directional) {
                if(AS3 > XS1+XS2) {
                    XS1 = get_AS(*(read3)); 
                    XS2 = get_AS(*(read3+1)); 
                }
                if(AS4 > XS1+XS2) {
                    XS1 = get_AS(*(read4)); 
                    XS2 = get_AS(*(read4+1)); 
                }
            }
        }
        if(best_node == 3) {
            if(AS1 > XS1+XS2) {
                XS1 = get_AS(*(read1)); 
                XS2 = get_AS(*(read1+1)); 
            }
            if(AS2 > XS1+XS2) {
                XS1 = get_AS(*(read2)); 
                XS2 = get_AS(*(read2+1)); 
            }
            if(AS4 > XS1+XS2) {
                XS1 = get_AS(*(read4)); 
                XS2 = get_AS(*(read4+1)); 
            }
        }
        if(best_node == 4) {
            if(AS1 > XS1+XS2) {
                XS1 = get_AS(*(read1)); 
                XS2 = get_AS(*(read1+1)); 
            }
            if(AS2 > XS1+XS2) {
                XS1 = get_AS(*(read2)); 
                XS2 = get_AS(*(read2+1)); 
            }
            if(AS3 > XS1+XS2) {
                XS1 = get_AS(*(read3)); 
                XS2 = get_AS(*(read3+1)); 
            }
        }
        MAPQ = calc_MAPQ_BT2(get_AS(tmp_read1)+get_AS(tmp_read2), XS1+XS2, scMin1+scMin2);
        MAPQ = (MAPQ < tmp_read1->core.qual) ? MAPQ : tmp_read1->core.qual; //Otherwise, a mapping can get worse but a score better!
        tmp_read1->core.qual = MAPQ;
        tmp_read2->core.qual = MAPQ;
        //replace/add the XS tag
        if(XS1 >= scMin1) {
            if(bam_aux_get(tmp_read1, "XS")) bam_aux_del(tmp_read1, bam_aux_get(tmp_read1, "XS"));
            bam_aux_append(tmp_read1, "XS", 'i', 4, (uint8_t*) &XS1);
        }
        if(XS2 >= scMin2) {
            if(bam_aux_get(tmp_read2, "XS")) bam_aux_del(tmp_read2, bam_aux_get(tmp_read2, "XS"));
            bam_aux_append(tmp_read2, "XS", 'i', 4, (uint8_t*) &XS2);
        }
    }
    free(kXX1);
    free(kXX2);
    free(kXM1);
    free(kXM2);
    return best_node;
}

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
bam1_t * update_read(struct packed_struct *first, int offset) {
    bam1_t *pbam1_t;
    uint8_t *data;
    bam1_t *new_copy = bam_init1();

    if(offset == 0) {
        pbam1_t = (bam1_t *) first->next->packed;
    } else {
        pbam1_t = (bam1_t *) first->next->next->packed;
    }
    data = (uint8_t *) (pbam1_t+1);
    pbam1_t->data = data;
    bam_copy1(new_copy, pbam1_t);
    free(pbam1_t);
    if(offset == 0) {
        first->next->packed = (void *) new_copy;
    } else {
        first->next->next->packed = (void *) new_copy;
    }
    return new_copy;
}

/*******************************************************************************
*
*   The master node function.
*
*   void *a: Actually a int*, the thread_id
*
*******************************************************************************/
void * master_processer_thread(void *a) {
    int thread_id = 0, best_node, i;
    int times = (config.paired) ? 2 : 1;
    char **seq = malloc(sizeof(char *) * 2);
    *(seq) = malloc(sizeof(char)*MAXREAD);
    *(seq+1) = malloc(sizeof(char)*MAXREAD);
    bam1_t **node1_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node2_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node3_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node4_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **best_read = NULL;
    time_t now;

    //Metrics
    metrics_struct *metrics = malloc(sizeof(metrics_struct));
    metrics->t_reads = 0;
    metrics->m_reads_OT = 0;
    metrics->m_reads_OB = 0;
    metrics->m_reads_CTOT = 0;
    metrics->m_reads_CTOB = 0;
    metrics->t_CpG = 0;
    metrics->m_CpG = 0;
    metrics->t_CHG = 0;
    metrics->m_CHG = 0;
    metrics->t_CHH = 0;
    metrics->m_CHH = 0;

    //Process read i/o
    while(1) {
        while(!is_ready(node1, 0));
        if(is_finished(node1)) break;
        *(node1_read) = update_read(node1, 0);
        if(config.paired) {
            while(!is_ready(node1, 1));
            *(node1_read+1) = update_read(node1, 1);
        }
        while(!is_ready(node2, 0));
        *(node2_read) = update_read(node2, 0);
        if(config.paired) {
            while(!is_ready(node2, 1));
            *(node2_read+1) = update_read(node2, 1);
        }
        if(!config.directional) {
            while(!is_ready(node3, 0));
            *node3_read = update_read(node3, 0);
            if(config.paired) {
                while(!is_ready(node3, 1));
                *(node3_read+1) = update_read(node3, 1);
            }
            while(!is_ready(node4, 0));
            *node4_read = update_read(node4, 0);
            if(config.paired) {
                while(!is_ready(node4, 1));
                *(node4_read+1) = update_read(node4, 1);
            }
        }
        metrics->t_reads++;

        //Give some output, it's a bit misleading as the count is actually only for this thread and it'll only display for thread 0.
        if(!config.quiet) {
            if(thread_id == 0) {
                if((metrics->t_reads) % 100000 == 0) {
                    now = time(NULL);
                    printf("%llu reads %s", metrics->t_reads, ctime(&now)); fflush(stdout);
                }
            }
        }

        get_seq(*seq, zip1);
        if(config.paired) get_seq(*(seq+1), zip2);

        //Process the reads
        if(!config.paired) {
            best_node = process_single(*node1_read, *node2_read, *node3_read, *node4_read, *seq); //Output is stored in read1
        } else {
            best_node = process_paired(node1_read, node2_read, node3_read, node4_read, seq); //Output is stored in read
        }
        if(best_node == 1) {
            best_read = node1_read;
            if(!((*best_read)->core.flag & BAM_FUNMAP)) metrics->m_reads_OT++;
        } else if(best_node == 2) {
            best_read = node2_read;
            if(!((*best_read)->core.flag & BAM_FUNMAP)) metrics->m_reads_OB++;
        } else if(best_node == 3) {
            best_read = node3_read;
            if(!((*best_read)->core.flag & BAM_FUNMAP)) metrics->m_reads_CTOT++;
        } else if(best_node == 4) {
            best_read = node4_read;
            if(!((*best_read)->core.flag & BAM_FUNMAP)) metrics->m_reads_CTOB++;
        }

        //Store the reads
        if(!config.paired) {
            if(!((*(best_read))->core.flag & BAM_FUNMAP)) {
                bam_write1(OUTPUT_BAM, *(best_read));
                update_counts(*(best_read), metrics);
            } else {
                if(config.unmapped) write_unmapped(unmapped1, *(best_read));
            }
        } else {
            if(!((*(best_read))->core.flag & BAM_FUNMAP) && !((*(best_read+1))->core.flag & BAM_FUNMAP)) {
                bam_write1(OUTPUT_BAM, *(best_read));
                update_counts(*(best_read), metrics);
                bam_write1(OUTPUT_BAM, *(best_read+1));
                update_counts(*(best_read+1), metrics);
            } else {
                if(config.unmapped) {
                    write_unmapped(unmapped1, *(best_read));
                    write_unmapped(unmapped2, *(best_read+1));
                }
            }
        }

        //Remove the processed reads
        for(i=0; i<times; i++){
            remove_element(node1);
            remove_element(node2);
            if(!config.directional) {
                remove_element(node3);
                remove_element(node4);
            }
        }
    }

    

    //Update the global metrics
    t_reads += metrics->t_reads;
    m_reads_OT += metrics->m_reads_OT;
    m_reads_OB += metrics->m_reads_OB;
    m_reads_CTOT += metrics->m_reads_CTOT;
    m_reads_CTOB += metrics->m_reads_CTOB;
    t_CpG += metrics->t_CpG;
    m_CpG += metrics->m_CpG;
    t_CHG += metrics->t_CHG;
    m_CHG += metrics->m_CHG;
    t_CHH += metrics->t_CHH;
    m_CHH += metrics->m_CHH;

    //Clean up
    free(*(seq)); free(*(seq+1)); free(seq);
    free(metrics);
    bam_header_destroy(global_header);
    free(node1_read);
    free(node2_read);
    free(node3_read);
    free(node4_read);
    destroy_list(node1);
    destroy_list(node2);
    destroy_list(node3);
    destroy_list(node4);
    return NULL;
}
