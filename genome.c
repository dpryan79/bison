#include "bison.h"

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
void get_seq(char *seq, FILE *fp) {
    char *line = malloc(MAXREAD*sizeof(char));
    assert(fgets(line, MAXREAD, fp) != NULL);
    assert(fgets(line, MAXREAD, fp) != NULL);
    *(line+strlen(line)-1) = '\0'; //remove the \n
    strcpy(seq, line);
    assert(fgets(line, MAXREAD, fp) != NULL);
    assert(fgets(line, MAXREAD, fp) != NULL);
    free(line);
}

/******************************************************************************
*
*   Return the length of a given chromosome.
*
*   char *chrom: the chromosome of interest
*
*******************************************************************************/
unsigned long long genome_chrom_length(char *chrom) {
    int i;
    unsigned long long output = 0;

    for(i=0; i<chromosomes.nchromosomes; i++) {
        if(strcmp(chromosomes.chromosome[i]->chrom, chrom) == 0) {
            output = chromosomes.chromosome[i]->length;
            break;
        }
    }
    return output;
}

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
*******************************************************************************/
char * get_genomic_context(unsigned long long offset, unsigned long long position, int change, unsigned long long chrom_length) {
    int i;
    char *output = calloc(4, sizeof(char));

    if(change > 0) {
        for(i=0; i<3; i++) {
            if(position+i < chrom_length) {
                *(output+i) = toupper(*(chromosomes.genome+offset+position+i));
            } else {
                *(output+i) = 'N';
            }
        }
    } else {
        for(i=0; i<3; i++) {
            if(position-2+i >= 0) {
                *(output+i) = toupper(*(chromosomes.genome+offset+position-2+i));
            } else {
                *(output+i) = 'N';
            }
        }
    }
    return output;
}
