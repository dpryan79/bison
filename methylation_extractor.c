#include "bison.h"
#include "sam.h"

//Eh, this is simple enough for a small program
int storeCpG, storeCHG, storeCHH, min_Phred;

//Inclusion bounds
int OT[4], OB[4], CTOT[4], CTOB[4];

//This will hold the output file handles (some of which can be NULL)
struct of_struct {
    FILE *CpG;
    FILE *CHG;
    FILE *CHH;
};

//This stores an individual methylation call
typedef struct {
    int32_t tid;
    int32_t start;
    _Bool strand; //+/- == 1/0
    unsigned int type; //This would normally be 0 (unmethylated) or 1 (methylated)
} Site;

//This struct hold an array of methylation calls that will need to be sorted
typedef struct {
    Site *CpG;
    Site *CHG;
    Site *CHH;
    int num_CpG;
    int max_CpG;
    int num_CHG;
    int max_CHG;
    int num_CHH;
    int max_CHH;
    int only_CpG;
    int only_CHG;
    int only_CHH;
} Sites;

struct list_struct {
    int32_t tid;
    int32_t pos; //negative positions are - strand, otherwise, + strand
    unsigned int n_methylated;
    unsigned int n_unmethylated;
    struct list_struct *next;
};

//Linked lists holding the final methylation calls
struct list_struct *CpGlist, *CHGlist, *CHHlist;

//Initialize a linked list
struct list_struct* init_list() {
    struct list_struct *output = calloc(1, sizeof(struct list_struct));
    struct list_struct *next = calloc(1, sizeof(struct list_struct));
    output->next = next;
    output->tid = -1;
    output->pos = -1;

    next->next = NULL;;
    next->tid = INT_MAX;
    next->pos = INT_MAX;

    return output;
}

//Destroy the linked list
void destroy_methyl_list(struct list_struct *list) {
    struct list_struct *next = list->next;
    struct list_struct *current = list;

    while(next != NULL) {
        next = current->next;
        free(current);
        current = next;
   }
}

//Insert a new methylation call into the linked list
struct list_struct* insert_call(struct list_struct *current, Site *site) {
    struct list_struct *next = current->next;
    struct list_struct *new = malloc(sizeof(struct list_struct));

    new->next = next;
    new->tid = site->tid;
    new->pos = (site->strand) ? site->start : -1 * (site->start);
    if(site->type) {
        new->n_methylated = 1;
        new->n_unmethylated = 0;
    } else {
        new->n_methylated = 0;
        new->n_unmethylated = 1;
    }
    current->next = new;
    return new;
}

/*******************************************************************************
*
*  Initialize a Sites structure
*
*******************************************************************************/
Sites* init_sites() {
    Sites *output = malloc(sizeof(Sites));
    output->CpG = malloc(sizeof(Site)*1000000);
    output->CHG = malloc(sizeof(Site)*1000000);
    output->CHH = malloc(sizeof(Site)*1000000);
    output->num_CpG = 0;
    output->max_CpG = 1000000;
    output->num_CHG = 0;
    output->max_CHG = 1000000;
    output->num_CHH = 0;
    output->max_CHH = 1000000;
    output->only_CpG = 0;
    output->only_CHG = 0;
    output->only_CHH = 0;
    return output;
}

/*******************************************************************************
*
*  Free space used by a Sites structure
*
*******************************************************************************/
void destroy_sites(Sites *p) {
    free(p->CpG);
    free(p->CHG);
    free(p->CHH);
    free(p);
}

/*******************************************************************************
*
*  Site sorting comparison function used by qsort in sort calls
*
*******************************************************************************/
int site_comparison(const void *p1, const void *p2) {
    Site *site1 = (Site *) p1;
    Site *site2 = (Site *) p2;
    int output = 0;

    if(site1->tid == site2->tid) {
        if(site1->start == site2->start) {
            output = 0;
        } else {
            output = site1->start - site2->start;
        }
    } else {
        output = strcmp(global_header->target_name[site1->tid],global_header->target_name[site2->tid]);
    }
    return output;
}

/*******************************************************************************
*
*  Sort methylation sites according to chromosome and start position
*
*******************************************************************************/
void sort_sites(Sites *sites, int which) {
    if(which == 1) qsort((void *) sites->CpG, (size_t) sites->num_CpG, sizeof(Site), site_comparison);
    else if(which == 2) qsort((void *) sites->CHG, (size_t) sites->num_CHG, sizeof(Site), site_comparison);
    else if(which == 3) qsort((void *) sites->CHH, (size_t) sites->num_CHH, sizeof(Site), site_comparison);
}

void merge_calls(Sites *sites, int which) {
    Site *type;
    int nsites=0, i=0;
    struct list_struct *olist=NULL, *current=NULL;

    if(which == 1) {
        type = sites->CpG;
        olist = CpGlist;
        nsites = sites->num_CpG;
    } else if(which == 2) {
        type = sites->CHG;
        olist = CHGlist;
        nsites = sites->num_CHG;
    } else if(which == 3) {
        type = sites->CHH;
        olist = CHHlist;
        nsites = sites->num_CHH;
    }

    //Take care of the first call
    current = olist;
    while(i<nsites) {
        if(current->tid == type[i].tid) {
            if(abs(current->pos) == abs(type[i].start)) {
                if(type[i].type == 1) current->n_methylated++;
                else current->n_unmethylated++;
                i++;
            } else if(abs(current->next->pos) > type[i].start) {
                current = insert_call(current, type+i);
                i++;
            } else {
                if(current->next->tid == type[i].tid) {
                    current = current->next;
                } else {
                    current = insert_call(current, type+i);
                    i++;
                }
            }
        } else {
            if(current->next->tid == INT_MAX) {
                current = insert_call(current, type+i);
                i++;
            } else if(current->next->tid == type[i].tid) { //Changing chromosomes
                if(abs(current->next->pos) > type[i].start) {
                    current = insert_call(current, type+i);
                    i++;
                } else {
                    current = current->next;
                }
            } else if(strcmp(global_header->target_name[type[i].tid], global_header->target_name[current->next->tid]) < 0) {
                current = insert_call(current, type+i);
                i++;
            } else {
                current = current->next;
            }
        }
    }

    //Reset the appropriate counter
    if(which == 1) sites->num_CpG = 0;
    else if(which == 2) sites->num_CHG = 0;
    else if(which == 3) sites->num_CHH = 0;
}


/*******************************************************************************
*
*  This will write the actual output, return 1 on success and 0 on error.
*
*******************************************************************************/
int process_call(int32_t tid, unsigned int position, char call, Sites *sites, char strand) {
    Site *site;

    if(call == 'Z' || call == 'z') { //CpG (methylated == Z, unmethylated = z)
        if(sites->only_CHG || sites->only_CHH) return 1;
        site = sites->CpG+sites->num_CpG;
        site->tid = tid;
        site->strand = (strand == '+') ? 1 : 0;
        site->start = position;
        site->type = (call == 'Z') ? 1 : 0;
        (sites->num_CpG)++;
    } else if(call == 'H' || call == 'h') { //CHH (methylated == H, unmethylated == h)
        if(sites->only_CpG || sites->only_CHG) return 1;
        site = sites->CHH+sites->num_CHH;
        site->tid = tid;
        site->strand = (strand == '+') ? 1 : 0;
        site->start = position;
        site->type = (call == 'H') ? 1 : 0;
        (sites->num_CHH)++;
    } else if(call == 'X' || call == 'x') { //CHG (methylated == X, unmethylated == x)
        if(sites->only_CpG || sites->only_CHH) return 1;
        site = sites->CHG+sites->num_CHG;
        site->tid = tid;
        site->strand = (strand == '+') ? 1 : 0;
        site->start = position;
        site->start = position;
        site->type = (call == 'X') ? 1 : 0;
        (sites->num_CHG)++;
    } else {
        printf("(1) Got an unknown character in the XM string of a read: %c\n",call);
        return 0;
    }
    return 1;
}

/*******************************************************************************
*
*  Process either a single-end read or a non-overlapping paired-end read.
*
*  Return 1 on success and 0 on error.
*
*******************************************************************************/
int extractor_process_single(bam1_t *read, Sites *sites) {
    unsigned long long *positions = NULL;
    char *XM = bam_aux2Z(bam_aux_get(read,"XM"));
    char *XR = bam_aux2Z(bam_aux_get(read, "XR"));
    char *XG = bam_aux2Z(bam_aux_get(read, "XG"));
    uint8_t *QUAL = bam1_qual(read);
    char strand = (strcmp(bam_aux2Z(bam_aux_get(read,"XG")), "CT") == 0) ? '+' : '-';
    char call;
    int i, start = 0, end = strlen(XM); //These may be overridden

    /***************************************************************************
    *
    *  Do we need to increase the size of anything pointed to by sites?
    *
    ***************************************************************************/
    if(storeCpG) {
        if(sites->num_CpG + 100000 > sites->max_CpG) {
            sites->CpG = realloc(sites->CpG, (sites->max_CpG+100000)*sizeof(Site));
            sites->max_CpG += 100000;
        }
    }
    if(storeCHG) {
        if(sites->num_CHG + 100000 > sites->max_CHG) {
            sites->CHG = realloc(sites->CHG, (sites->max_CHG+100000)*sizeof(Site));
            sites->max_CHG += 100000;
        }
    }
    if(storeCHH) {
        if(sites->num_CHH + 100000 > sites->max_CHH) {
            sites->CHH = realloc(sites->CHH, (sites->max_CHH+100000)*sizeof(Site));
            sites->max_CHH += 100000;
        }
    }

    positions = calculate_positions(read);

    //Should we override "start" and "end"?
    if(read->core.flag & BAM_FREAD2) { //#2
        if(strcmp(XR, "GA") == 0 && strcmp(XG, "CT") == 0) { //OT
            if(OT[2] != 0) start = OT[2];
            if(OT[3] != 0) {
                if(end > OT[3]) end = OT[3];
            }
        } else if(strcmp(XR, "GA") == 0 && strcmp(XG, "GA") == 0) { //OB
            if(OB[2] != 0) start = OB[2];
            if(OB[3] != 0) {
                if(end > OB[3]) end = OB[3];
            }
        } else if(strcmp(XR, "CT") == 0 && strcmp(XG, "CT") == 0) { //CTOT
            if(CTOT[2] != 0) start = CTOT[2];
            if(CTOT[3] != 0) {
                if(end > CTOT[3]) end = CTOT[3];
            }
        } else if(strcmp(XR, "CT") == 0 && strcmp(XG, "CT") == 0) { //CTOT
            if(CTOB[2] != 0) start = CTOB[2];
            if(CTOB[3] != 0) {
                if(end > CTOB[3]) end = CTOB[3];
            }
        }
    } else { //#1
        if(strcmp(XR, "CT") == 0 && strcmp(XG, "CT") == 0) { //OT
            if(OT[0] != 0) start = OT[0];
            if(OT[1] != 0) {
                if(end > OT[1]) end = OT[1];
            }
        } else if(strcmp(XR, "CT") == 0 && strcmp(XG, "GA") == 0) { //OB
            if(OB[0] != 0) start = OB[0];
            if(OB[1] != 0) {
                if(end > OB[1]) end = OB[1];
            }
        } else if(strcmp(XR, "GA") == 0 && strcmp(XG, "CT") == 0) { //CTOT
            if(CTOT[0] != 0) start = CTOT[0];
            if(CTOT[1] != 0) {
                if(end > CTOT[1]) end = CTOT[1];
            }
        } else if(strcmp(XR, "GA") == 0 && strcmp(XG, "GA") == 0) { //CTOB
            if(CTOB[0] != 0) start = CTOB[0];
            if(CTOB[1] != 0) {
                if(end > CTOB[1]) end = CTOB[1];
            }
        }
    }

    for(i=start; i<end; i++) {
        if(*(positions+i) != ULLONG_MAX) {
            if(*(XM+i) != '.') {
                if(*(QUAL+i) < min_Phred) continue;
                call = *(XM+i);
                if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                    if(!process_call(read->core.tid, *(positions+i), *(XM+i), sites, strand)) {
                        printf("(2) Got an unknown character (%i) in the XM string of a single-ended read: %s\n", i, XM);
                        free(positions);
                        return 0;
                    }
                }
            }
        }
    }
    free(positions);
    return 1;
}

/*******************************************************************************
*
*  Process either overlapping paired-end reads
*
*  Return 1 on success and 0 on error.
*
*******************************************************************************/
int extractor_process_overlapping(bam1_t *read1, bam1_t *read2, Sites *sites) {
    unsigned long long *positions1 = calculate_positions(read1), *positions2 = calculate_positions(read2);
    char strand = (strcmp(bam_aux2Z(bam_aux_get(read1, "XG")), "CT") == 0) ? '+' : '-';
    char call;
    char *XR = bam_aux2Z(bam_aux_get(read1,"XR"));
    char *XG = bam_aux2Z(bam_aux_get(read1,"XG"));
    char *XM1 = bam_aux2Z(bam_aux_get(read1,"XM"));
    char *XM2 = bam_aux2Z(bam_aux_get(read2,"XM"));
    int i, j, end1 = (int) read1->core.l_qseq, end2 = (int) read2->core.l_qseq;
    int start1 = 0, start2 = 0;

    /***************************************************************************
    *
    *  Do we need to increase the size of anything pointed to by sites?
    *
    ***************************************************************************/
    if(sites->num_CpG + 100000 > sites->max_CpG) {
        sites->CpG = realloc(sites->CpG, (sites->max_CpG+100000)*sizeof(Site));
        sites->max_CpG += 100000;
    }
    if(sites->num_CHG + 100000 > sites->max_CHG) {
        sites->CHG = realloc(sites->CHG, (sites->max_CHG+100000)*sizeof(Site));
        sites->max_CHG += 100000;
    }
    if(sites->num_CHH + 100000 > sites->max_CHH) {
        sites->CHH = realloc(sites->CHH, (sites->max_CHH+100000)*sizeof(Site));
        sites->max_CHH += 100000;
    }

    //Should we override start1,start2, end1 and end2?
    if(strcmp(XR, "CT") == 0 && strcmp(XG, "CT") == 0) { //OT
        if(OT[0] != 0) start1 = OT[0];
        if(OT[1] != 0) {
            if(end1 > OT[1]) end1 = OT[1];
        }
        if(OT[2] != 0) start2 = OT[2];
        if(OT[3] != 0) {
            if(end2 > OT[3]) end2 = OT[3];
        }
    } else if(strcmp(XR, "CT") == 0 && strcmp(XG, "GA") == 0) { //OB
        if(OB[0] != 0) start1 = OB[0];
        if(OB[1] != 0) {
            if(end1 > OB[1]) end1 = OB[1];
        }
        if(OB[2] != 0) start2 = OB[2];
        if(OB[3] != 0) {
            if(end2 > OB[3]) end2 = OB[3];
        }
    } else if(strcmp(XR, "GA") == 0 && strcmp(XG, "CT") == 0) { //CTOT
        if(CTOT[0] != 0) start1 = CTOT[0];
        if(CTOT[1] != 0) {
            if(end1 > CTOT[1]) end1 = CTOT[1];
        }
        if(CTOT[2] != 0) start2 = CTOT[2];
        if(CTOT[3] != 0) {
            if(end2 > CTOT[3]) end2 = CTOT[3];
        }
    } else if(strcmp(XR, "GA") == 0 && strcmp(XG, "GA") == 0) { //CTOB
        if(CTOB[0] != 0) start1 = CTOB[0];
        if(CTOB[1] != 0) {
            if(end1 > CTOB[1]) end1 = CTOB[1];
        }
        if(CTOB[2] != 0) start2 = CTOB[2];
        if(CTOB[3] != 0) {
            if(end2 > CTOB[3]) end2 = CTOB[3];
        }
    }
    i = start1;
    XM1 += start1;
    j = start2;
    XM2 += start2;
    while(*(positions1+i) == ULLONG_MAX) {
        i++;
        start1++;
        XM1++;
    }
    while(*(positions2+j) == ULLONG_MAX) {
        j++;
        start2++;
        XM2++;
    }
    while(*(positions1+end1-1) == ULLONG_MAX) end1--;
    while(*(positions2+end2-1) == ULLONG_MAX) end2--;

    /***************************************************************************
    *
    *  If there is a 5' overhang when comparing the two sequences, then we
    *  should process that first before dealing with the overlap.
    *
    ***************************************************************************/
    if(*positions1 < *positions2) {
        while(*(positions1+i) < *positions2) {
            if(*(positions1+i) != ULLONG_MAX) {
                if(*XM1 != '.') {
                    call = *XM1;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read1->core.tid, *(positions1+i), *XM1, sites, strand)) {
                            printf("(3) Got an unknown character in the XM string: %s\n", XM1);
                            return 0;
                        }
                    }
                }
            }
            i++;
            XM1++;
            if(i == end1) break;
        }
    } else if(*positions2 < *positions1) {
        while(*(positions2+j) < *positions1) {
            if(*(positions2+j) != ULLONG_MAX) {
                if(*XM2 != '.') {
                    call = *XM2;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read2->core.tid, *(positions2+j), *XM2, sites, strand)) {
                            printf("(4) Got an unknown character in the XM string: %s\n", XM2);
                            return 0;
                        }
                    }
                }
            }
            j++;
            XM2++;
            if(j == end2) break;
        }
    }

    //We are now up to the overlapping section
    while((i<end1) && (j<end2)) {
        //Ensure we're on the same position
        if(*(positions1+i) != ULLONG_MAX && *(positions2+j) != ULLONG_MAX) {
            //Deal with InDels at the beginning of homopolymer repeats which may be different in the reads
            if(*(positions1+i) != *(positions2+j)) {
                if(*(positions1+i) < *(positions2+j)) {
                    if(*XM1 != '.') {
                        call = *XM1;
                        if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                            if(!process_call(read1->core.tid, *(positions1+i), *XM1, sites, strand)) {
                                printf("(5a) Got an unknown character in the XM string: %s\n", XM1);
                                return 0;
                            }
                        }
                    }
                    XM1++;
                    i++;
                } else {
                    if(*XM2 != '.') {
                        call = *XM2;
                        if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                            if(!process_call(read2->core.tid, *(positions2+j), *XM2, sites, strand)) {
                                printf("(5b) Got an unknown character in the XM string: %s\n", XM2);
                                return 0;
                            }
                        }
                    }
                    XM2++;
                    j++;
                }
                continue;
            }

            if(*XM1 == *XM2) {
                if(*XM1 != '.') {
                    call = *XM1;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read1->core.tid, *(positions1+i), *XM1, sites, strand)) {
                            printf("(6) Got an unknown character in the XM string: %s\n", XM1);
                            return 0;
                        }
                    }
                }
            } else { //bison will call '.' if there is an N in a read or an impossible conversion, so whichever read has a call is correct (the call becomes '.' if the reads have different calls)
                if(*XM2 != '.' && *XM1 == '.') {
                    call = *XM2;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read2->core.tid, *(positions2+j), *XM2, sites, strand)) {
                            printf("(7) Got an unknown character in the XM string: %s\n", XM2);
                            return 0;
                        }
                    }
                } else if(*XM1 != '.' && *XM2 == '.') {
                    call = *XM1;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read1->core.tid, *(positions1+i), *XM1, sites, strand)) {
                            printf("(8) Got an unknown character in the XM string: %s\n", XM1);
                            return 0;
                        }
                    }
                }
            }
            XM1++;
            XM2++;
            i++;
            j++;
        } else {
            if(*(positions1+i) == ULLONG_MAX) { 
                XM1++;
                i++;
            }
            if(*(positions2+j) == ULLONG_MAX) { 
                XM2++;
                j++;
            }
        }
    }

    if(i >= end1 && j >= end2) {
        free(positions1);
        free(positions2);
        return 1;
    }
    if(i >= end1) {
        while(j<end2) {
            if(*(positions2+j) != ULLONG_MAX) {
                if(*XM2 != '.') {
                    call = *XM2;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read2->core.tid, *(positions2+j), *XM2, sites, strand)) {
                            printf("(12) Got an unknown character in the XM string: %s\n", XM2);
                            free(positions1);
                            free(positions2);
                            return 0;
                        }
                    }
                }
            }
            XM2++;
            j++;
        }
    } else {
        while(i<end1) {
            if(*(positions1+i) != ULLONG_MAX) {
                if(*XM1 != '.') {
                    call = *XM1;
                    if(((call == 'Z' || call == 'z') && storeCpG) || ((call == 'X' || call == 'x') && storeCHG) || ((call == 'H' || call == 'h') && storeCHH)) {
                        if(!process_call(read1->core.tid, *(positions1+i), *XM1, sites, strand)) {
                            printf("(13) Got an unknown character in the XM string: %s\n", XM1);
                            free(positions1);
                            free(positions2);
                            return 0;
                        }
                    }
                }
            }
            XM1++;
            i++;
        }
    }
    free(positions1);
    free(positions2);
    return 1;
}

//Remove methylation calls on 5' ends of reads
void trim5(bam1_t *read, int digest_types) {
    int MspI = digest_types & 1;
    int TaqI = digest_types & 2;
    unsigned long long offset = genome_offset(lookup_chrom(read), read->core.pos);
    char *sequence;
    char *XM = bam_aux2Z(bam_aux_get(read, "XM"));
    int i;

    for(i=0; i<2; i++) {
        sequence = chromosomes.genome+offset+read->core.pos-i;
        if(MspI) {
            if(strcmp(sequence, "CCGG")) {
                *(XM+2-i) = '.';
                break;
            }
        }
        if(TaqI) {
            if(strcmp(sequence, "TCGA")) {
                *(XM+2-i) = '.';
                break;
            }
        }
    }
}

//Remove methylation calls on 3' ends of reads
void trim3(bam1_t *read, int digest_types) {
    int MspI = digest_types & 1;
    int TaqI = digest_types & 2;
    unsigned long long offset = genome_offset(lookup_chrom(read), read->core.pos);
    char *sequence;
    char *XM = bam_aux2Z(bam_aux_get(read, "XM"));
    uint32_t end = bam_calend(&(read->core), bam1_cigar(read));
    int i, len = strlen(XM);

    for(i=0; i<2; i++) {
        sequence = chromosomes.genome+offset+end-2-i;
        if(MspI) {
            if(strcmp(sequence, "CCGG")) {
                *(XM+len-2-i) = '.';
                break;
            }
        }
        if(TaqI) {
            if(strcmp(sequence, "TCGA")) {
                *(XM+len-2-i) = '.';
                break;
            }
        }
    }
}

void process_RRBS_read(bam1_t *read1, bam1_t *read2, int digest_types) {
    if(strncmp(bam_aux2Z(bam_aux_get(read1, "XG")), "CT", 2) == 0) { //OT or CTOT
        trim3(read1, digest_types);
        if(read1->core.flag & BAM_FPAIRED) trim3(read2, digest_types);
    } else { //OB or CTOB
        trim5(read1, digest_types);
        if(read1->core.flag & BAM_FPAIRED) trim5(read2, digest_types);
    }
}

void write_sites(struct of_struct *of, int which) {
    struct list_struct *list = NULL;
    int mpercent;
    FILE *f = NULL;

    if(which == 1) {
        list = CpGlist->next;
        f = of->CpG;
    } else if(which == 2) {
        list = CHGlist->next;
        f = of->CHG;
    } else if(which == 3) {
        list = CHHlist->next;
        f = of->CHH;
    }

    while(list->next != NULL) {
        mpercent = (int) (1000 * ((float) list->n_methylated)/(float)(list->n_methylated + list->n_unmethylated));
        fprintf(f, "%s\t%u\t%u\t%i\t%u\t%u\n", global_header->target_name[list->tid], \
            abs(list->pos), abs(list->pos)+1, mpercent, list->n_methylated, list->n_unmethylated);
        list = list->next;
    }
}

//Generate output file names and open them for writing
void generate_output_names(char *ifile, struct of_struct *of) {
    char *p, *tmp = strdup(ifile);
    char *oname = NULL;

    //Generate the basename by stripping off .sam or .bam
    p = strrchr(tmp, '.');
    if(strcmp(p, ".sam") == 0 || strcmp(p, ".bam") == 0) *p = '\0';
    oname = malloc(sizeof(char) * (strlen(tmp) + strlen("_CpG.bedGraph ")));

    if(storeCpG) {
        sprintf(oname, "%s_CpG.bedGraph", tmp);
        printf("CpG counts will be written to %s\n", oname);
        of->CpG = fopen(oname, "w");
    }
    if(storeCHG) {
        sprintf(oname, "%s_CHG.bedGraph", tmp);
        printf("CHG counts will be written to %s\n", oname);
        of->CHG = fopen(oname, "w");
    }
    if(storeCHH) {
        sprintf(oname, "%s_CHH.bedGraph", tmp);
        printf("CHH counts will be written to %s\n", oname);
        of->CHH = fopen(oname, "w");
    }

    free(tmp);
    free(oname);
}

//Fill the inclusion bounds
void fill_bounds(char *str, int bounds[4]) {
    int i;
    char *p;

    for(i=0; i<4; i++) {
        if(i==0) {
            p = strtok(str, ",");
        } else {
            p = strtok(NULL, ",");
        }
        if(p == NULL) break;
        bounds[i] = atoi(p);
    }
}

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] genome_directory/ input.(sam|bam)\n", prog);
    printf("\n\
    Extract methylation information into a bedGraph file or files. By default,\n\
    only CpG metrics are output\n\
\n\
    -h            Print this message.\n\
\n\
    -q            Read MAPQ value must at least this for inclusion (default 10).\n\
                  Specify 0 to include everything.\n\
\n\
    -phred        Minimum Phred score that a base must have for its methylation\n\
                  state to be included in the output. The default is 5.\n\
\n\
    --MspI        Library was MspI digested.\n\
\n\
    --TaqI        Library was TaqI digested (this can be in addition to\n\
                  MspI digestion).\n\
\n\
    -no_CpG       Don't output CpG sites (they're output by default).\n\
\n\
    -CHH          Output CHH statistics.\n\
\n\
    -CHG          Output CHG statistics.\n\
\n\
    -OT           Bounds for the region of reads mapped to the original top\n\
                  strand to include. It is highly recommended that bison_mbias\n\
                  and/or bison_mbias2pdf be run so that approximate bounds can\n\
                  be generated for this. The format is \"-OT A,B,C,D\", where \"A\"\n\
                  is the 5'-most and \"B\" the 3'-most bound of the included\n\
                  region for read #1. \"C\" and \"D\" are the equivalent bounds for\n\
                  read #2. A value of 0 means to leave that portion of the read\n\
                  unbound (e.g., \"-OT 0,90,20,0\" will not include methylation\n\
                  calls after the 90th base on read #1 or before the 20th base\n\
                  on read #2). The default is \"-OT 0,0,0,0\", meaning that all\n\
                  methylation calls are included.\n\
\n\
    -OB           Like -OT, but for reads mapping to the original bottom strand.\n\
\n\
    -CTOT         Like -OT, but for reads mapping to the complementary to\n\
                  original top strand.\n\
\n\
    -CTOB         Like -OT, but for reads mapping to the complementary to\n\
                  original bottom strand.\n\
    --genome-size Many of the bison tools need to read the genome into memory.\n\
                  By default, they allocate 3000000000 bases worth of memory for\n\
                  this and increase that as needed. However, this can sometimes\n\
                  be far more than is needed (meaning wasted memory) or far too\n\
                  little (in which case the process can become quite slow). If\n\
                  you input the approximate size of your genome here (in bases),\n\
                  then you can maximize performance and minimize wasted space.\n\
                  It's convenient to round up a little.\n\
\n\
    -max-sites-size N    This option can increase or decrease memory\n\
                  requirements by changing the number of methylation calls\n\
                  stored in memory prior to sorting and merging. The default is\n\
                  50,000.\n");
}

int main(int argc, char *argv[]) {
    int i, max_sites_size = 50000, MspI = 0, TaqI = 0;
    int min_MAPQ = 10;
    samfile_t *fp = NULL;
    bam1_t *read1 = bam_init1(), *read2 = bam_init1();
    struct of_struct *of = calloc(1, sizeof(struct of_struct));
    unsigned int r1_pos = 0, r2_pos = 0, total_reads = 0;
    Sites *sites = init_sites();

    CpGlist = init_list();
    CHGlist = init_list();
    CHHlist = init_list();
    config.genome_dir = NULL;
    chromosomes.max_genome = 3000000000;
    chromosomes.nchromosomes = 0;
    storeCpG = storeCHG = storeCHH = 0;
    storeCpG = 1;
    min_Phred = 10;
    for(i=0; i<4; i++) OT[i] = OB[i] = CTOT[i] = CTOB[i] = 0;

    /* read in the file names */
    if(argc < 3) {
        usage(argv[0]);
        return 1;
    };
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-no_CpG") == 0) {
            storeCpG = 0;
        } else if(strcmp(argv[i], "-CHG") == 0) {
            storeCHG = 1;
        } else if(strcmp(argv[i], "-CHH") == 0) {
            storeCHH = 1;
        } else if(strcmp(argv[i], "--MspI") == 0) {
            MspI = 1;
        } else if(strcmp(argv[i], "--TaqI") == 0) {
            TaqI = 1;
        } else if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp(argv[i], "-q") == 0) {
            i++;
            min_MAPQ = atoi(argv[i]);
        } else if(strcmp(argv[i], "-phred") == 0) {
            i++;
            min_Phred = atoi(argv[i]);
        } else if(strcmp(argv[i], "-max-sites-size") == 0) {
            max_sites_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-OT") == 0) {
            fill_bounds(argv[++i], OT);
        } else if (strcmp(argv[i], "-OB") == 0) {
            fill_bounds(argv[++i], OT);
        } else if (strcmp(argv[i], "-CTOT") == 0) {
            fill_bounds(argv[++i], CTOT);
        } else if (strcmp(argv[i], "-CTOB") == 0) {
            fill_bounds(argv[++i], CTOB);
        } else if(config.genome_dir == NULL) {
            config.genome_dir = strdup(argv[i]);
            if(*(config.genome_dir+strlen(config.genome_dir)-1) != '/') {
                config.genome_dir = realloc(config.genome_dir, sizeof(char) * (strlen(config.genome_dir)+2));
                sprintf(config.genome_dir, "%s/", config.genome_dir);
            }
        } else if(strcmp(argv[i], "--genome-size") == 0) {
            i++;
            chromosomes.max_genome = strtoull(argv[i], NULL, 10);
        } else if(fp == NULL) {
            if(argv[i][strlen(argv[i])-3] == 'b') {
                fp = samopen(argv[i], "rb", NULL);
            } else {
                fp = samopen(argv[i], "r", NULL);
            }
            global_header = fp->header;
        } else {
            printf("Unknown parameter %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if(config.genome_dir == NULL || fp == NULL) {
        printf("Genome directory or SAM/BAM input file not specified!\n");
        usage(argv[0]);
    }

    //Generate the output names and open the output files
    generate_output_names(argv[argc-1], of);

    //Read in the genome
    printf("Allocating space for %llu characters\n", chromosomes.max_genome); fflush(stdout);
    chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
    *chromosomes.genome = '\0';
    if(chromosomes.genome == NULL) {
        printf("Could not allocate enough room to hold the genome!\n");
        return -1;
    }
    read_genome();

    //Process the reads
    while(samread(fp, read1) > 1) {
        if(read1->core.flag & BAM_FPAIRED) {
            samread(fp, read2);
            r1_pos = read1->core.pos+1;
            r2_pos = read2->core.pos+1;
        } else {
            r1_pos = read1->core.pos+1;
            r2_pos = INT_MAX;
        }
        if(read1->core.flag & BAM_FDUP) continue;
        if(read1->core.qual < min_MAPQ) continue;
        if(TaqI+MspI) process_RRBS_read(read1, read2, MspI+2*TaqI);

        //Are the reads even overlapping? If not, this is easy.
        if(r2_pos == INT_MAX) { //Unpaired read
            if(!extractor_process_single(read1, sites)) { printf("Error!\n"); break; }
        } else if(r1_pos < r2_pos && r1_pos + read1->core.l_qseq - 1 < r2_pos) { //No Overlap
            if(!extractor_process_single(read1, sites)) { printf("Error!\n"); break; }
            if(!extractor_process_single(read2, sites)) { printf("Error!\n"); break; }
        } else if(r2_pos < r1_pos && r2_pos + read2->core.l_qseq - 1 < r1_pos) { //No Overlap
            if(!extractor_process_single(read1, sites)) { printf("Error!\n"); break; }
            if(!extractor_process_single(read2, sites)) { printf("Error!\n"); break; }
        } else { //Overlap
            if(!extractor_process_overlapping(read1, read2, sites)) { printf("Error!\n"); break; }
        }

        if(sites->num_CpG >= max_sites_size) {
            sort_sites(sites, 1);
            merge_calls(sites, 1);
        }
        if(sites->num_CHG >= max_sites_size) {
            sort_sites(sites, 2);
            merge_calls(sites, 2);
        }
        if(sites->num_CHH >= max_sites_size) {
            sort_sites(sites, 3);
            merge_calls(sites, 3);
        }

        total_reads++;
        if(total_reads % 1000000 == 0) {
            printf("Processed %u reads\n", total_reads);
            fflush(stdout);
        }
    }
    if(samread(fp, read1) > 1) {
        printf("We must have exited on an error as there are still reads left\n");
        fflush(stdout);
    }

    //Do the final sort and merge
    if(sites->num_CpG) {
        sort_sites(sites, 1);
        merge_calls(sites, 1);
    }
    if(sites->num_CHG) {
        sort_sites(sites, 2);
        merge_calls(sites, 2);
    }
    if(sites->num_CHH) {
        sort_sites(sites, 3);
        merge_calls(sites, 3);
    }

    //Write output
    if(storeCpG) write_sites(of, 1);
    if(storeCHG) write_sites(of, 2);
    if(storeCHH) write_sites(of, 3);

    //Close things up
    if(of->CpG != NULL) fclose(of->CpG);
    if(of->CHG != NULL) fclose(of->CHG);
    if(of->CHH != NULL) fclose(of->CHH);
    free(of);
    free(chromosomes.genome);
    for(i=0; i<chromosomes.nchromosomes; i++) {
        free((chromosomes.chromosome[i])->chrom);
        free(*(chromosomes.chromosome+i));
    }
    free(chromosomes.chromosome);
    if(config.genome_dir != NULL) free(config.genome_dir);
    destroy_methyl_list(CpGlist);
    destroy_methyl_list(CHGlist);
    destroy_methyl_list(CHHlist);
    destroy_sites(sites);
    bam_destroy1(read1);
    bam_destroy1(read2);
    samclose(fp);

    return 0;
};
