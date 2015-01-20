#include "bison.h"

/*******************************************************************************
*
*  Create a position array to account for any InDels
*  This function assumes that the first base is not marked as an InDel or
*  clipped in any way. If that occurs then things will break.
*
*  The output needs to be free()d
*
*******************************************************************************/
unsigned long long *calculate_positions(bam1_t *read) {
    unsigned long long *positions = malloc(sizeof(unsigned long long) * (size_t)read->core.l_qseq);
    int i, j, offset = 0, op, op_len;
    uint32_t *CIGAR = bam_get_cigar(read);
    unsigned int previous_position = (unsigned int) read->core.pos;

    for(i=0; i<read->core.n_cigar; i++) {
        op = *(CIGAR+i) & 15;
        op_len = (*(CIGAR+i)) >> 4;
        for(j=0; j<op_len; j++) {
            if(op == 0 || op == 7 || op == 8) { //M, =, X
                *(positions+offset) = previous_position++;
                offset++;
            } else if(op == 1 || op == 4 || op == 5) { //I, S, H
                *(positions+offset) = ULLONG_MAX; //This sets a practical limit on a contig's length (though this should never occur in reality
                offset++;
            } else if(op == 2 || op == 3) { //D, N
                previous_position++;
            } else { //P
                fprintf(stderr,"We encountered a CIGAR operation that we're not ready to deal with in %s\n", bam_get_qname(read));
            }
        }
    }
    return positions;
}

/******************************************************************************
*
*   Given a struct dirent, return 1 if it's a fasta file. Otherwise, return 0.
*
*******************************************************************************/
int filter(const struct dirent *file) {
    char *p = strrchr(file->d_name, '.');

    if(p == NULL) return 0; //No file extension!
    if(strcmp(p, ".fa") == 0 || strcmp(p, ".fasta") == 0) return 1; //A fasta file
    return 0;
}

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
void read_genome() {
//    DIR *dir = opendir(config.genome_dir);
    FILE *fp;
    char *p, *line = malloc(sizeof(char)*MAXREAD), *fullpath = NULL;
    char *g = chromosomes.genome;
    struct dirent **files;
    unsigned long long offset = 0;
    unsigned long long length = 0;
    int end, nchromosomes, i, j, nfiles;
    chromosome_struct *chromosome = NULL;

    nfiles = scandir(config.genome_dir, &files, filter, alphasort);
    for(j=0; j<nfiles; j++) {
        //This is a fasta file that we need to read into the genome array and append a chromosome_struct onto chromosomes_struct
        fullpath = realloc(fullpath, sizeof(char)*(strlen(config.genome_dir)+strlen(files[j]->d_name)+1));
        sprintf(fullpath, "%s%s",config.genome_dir,files[j]->d_name);
        fp = fopen(fullpath, "r");
        if(!config.quiet) fprintf(stderr, "Reading in %s\n", fullpath);
        fflush(stderr);
        while(fgets(line, MAXREAD, fp) != NULL) {
            end=strlen(line);
            if(line[end-1] == '\n') line[end-1] = '\0';
            if(line[0] == '>') {
                //Store the length of the previous contig, if there was one
                if(chromosome != NULL) {
                    chromosome->length = length;
                }

                //Initialize a new chromosome_struct and lengthen the global chromosomes struct
                nchromosomes = ++chromosomes.nchromosomes;
                chromosomes.chromosome = realloc(chromosomes.chromosome, sizeof(chromosome_struct*) * nchromosomes);
                chromosomes.chromosome[nchromosomes-1] = malloc(sizeof(chromosome_struct));
                chromosome = chromosomes.chromosome[nchromosomes-1];
                chromosome->offset = offset;
                p = strchr(line, ' ');
                if(p != NULL) *p = '\0'; //If there's anything after the name, ignore it
                chromosome->chrom = malloc(sizeof(char)*strlen(line));
                strcpy(chromosome->chrom, (line+1)); //ignore the ">"
                length = 0;
                chromosome->offset = offset;
            } else {
                //Ensure that we have enough space in chromosomes.genome
                if(offset + 10000 >= chromosomes.max_genome) {
                    chromosomes.max_genome += 1000000;
                    chromosomes.genome = realloc(chromosomes.genome, sizeof(char) * chromosomes.max_genome);
                    g = chromosomes.genome + offset;
                 }
                 offset += end-1;
                 length += end-1;
                 for(i=0; i<strlen(line); i++) *(line+i) = toupper(*(line+i)); //Make everything upper case
                 strncpy(g, line, end);
                 g += end-1;
            }
        }
        //Store the last contig's length
        chromosome->length = length;
        if(!config.quiet) fprintf(stderr, "Finished %s\n", fullpath);
        fflush(stderr);
        fclose(fp);
        free(files[j]);
    }
    free(line);
    free(fullpath);
    free(files);
}

/******************************************************************************
*
*   Reverse complement a sequence (in place)
*
*   char *seq: the sequence
*
*******************************************************************************/
void reverse_complement(char *seq) {
    char *tmp = strdup(seq);
    char current, new;
    int i, j;

    for(i=0, j=strlen(tmp)-1; j>=0; i++, j--) {
        current = *(tmp+j);
        new = 'N';
        if(current == 'A' || current == 'a') new = 'T';
        if(current == 'T' || current == 't') new = 'A';
        if(current == 'C' || current == 'c') new = 'G';
        if(current == 'G' || current == 'g') new = 'C';
        *(seq+i) = new;
    }
    free(tmp);
}

/******************************************************************************
*
*   Determine the appropriate offset in chromosomes.genome
*
*   char *chrom: Chromosome name
*   int32_t pos: 0-based position on Chromosome. This is read->core.pos
*
*******************************************************************************/
unsigned long long genome_offset(char *chrom, int32_t pos) {
    int i;
    unsigned long long chrom_offset = 0;

    for(i=0; i<chromosomes.nchromosomes; i++) {
        if(strcmp(chromosomes.chromosome[i]->chrom, chrom) == 0) {
            chrom_offset = chromosomes.chromosome[i]->offset;
            chrom_offset += pos;
            break;
        }
    }

    if(chrom_offset == 0 && pos != 0) fprintf(stderr, "Unable to calculate the genomic offset for %s:%i!\n", chrom, (int) pos);
    return chrom_offset;
}

/******************************************************************************
*
*   Return a pointer to the chromosome name onto which a read maps.
*
*   bam1_t *read: The read in question
*
*******************************************************************************/
inline char *lookup_chrom(bam1_t *read) {
    int32_t tid = read->core.tid;
    return global_header->target_name[tid];
}
