#include "../bison.h"
#include "sam.h"

struct CpG {
    int tid;
    int start;
    int end;
    unsigned int n_methylated;
    unsigned int n_unmethylated;
};

//Copy CpG call information from one struct to another
//Just changing the pointer wouldn't work because we aren't malloc()ing things
inline void copy_CpG(struct CpG *target, struct CpG *source) {
    target->tid = source->tid;
    target->start = source->start;
    target->end = source->end;
    target->n_methylated = source->n_methylated;
    target->n_unmethylated = source->n_unmethylated;
}

//Print a single CpG struct that contains information on a full CpG to file
void print_CpG_line(FILE *of, struct CpG *CpG) {
    //#chrom, start, end
    fprintf(of, "%s\t%i\t%i", chromosomes.chromosome[CpG->tid]->chrom, CpG->start, CpG->end);
    //ratio, totalC, methC, strand, next
    fprintf(of, "\t%4.2f\t%u\t%u\t+\tG", ((double) CpG->n_methylated)/(CpG->n_methylated+CpG->n_unmethylated), \
        CpG->n_methylated+CpG->n_unmethylated, CpG->n_methylated);
    //Plus, totalC, methC, Minus, totalC, methC, localSeq
    fprintf(of, "\t+\t%u\t%u\t-\t0\t0\t\n", CpG->n_methylated+CpG->n_unmethylated, CpG->n_methylated);
}

//Given two CpGs (one of which may be NULL), print a line in MOAB format
void print_line(FILE *of, struct CpG *first, struct CpG *second) {
    char strand;
    int tid, start, end;
    unsigned int n_total, n_meth;

    //Determine full-CpG metrics
    if(first == NULL) {
        strand = '-';
        tid = second->tid;
        start = second->start-1;
        end = second->end;
        n_total = second->n_methylated + second->n_unmethylated;
        n_meth = second->n_methylated;
    } else if(second == NULL) {
        strand = '+';
        tid = first->tid;
        start = first->start;
        end = first->end+1;
        n_total = first->n_methylated + first->n_unmethylated;
        n_meth = first->n_methylated;
    } else {
        strand = 'B';
        tid = first->tid;
        start = first->start;
        end = second->end;
        n_total = first->n_methylated + first->n_unmethylated;
        n_total += second->n_methylated + second->n_unmethylated;
        n_meth = first->n_methylated + second->n_methylated;
    }

    //#chrom, start, end
    fprintf(of, "%s\t%i\t%i", chromosomes.chromosome[tid]->chrom, start, end);
    //ratio, totalC, methC, strand, next
    fprintf(of, "\t%4.2f\t%u\t%u\t%c\tG", ((double) n_meth)/n_total, n_total, n_meth, strand);
    //Plus, totalC, methC
    if(first != NULL) {
        fprintf(of, "\t+\t%u\t%u", first->n_methylated + first->n_unmethylated, first->n_methylated);
    } else {
        fprintf(of, "\t+\t0\t0");
    }
    //Minus, totalC, methC, localSeq
    if(second != NULL) {
        fprintf(of, "\t-\t%u\t%u\t\n", second->n_methylated + second->n_unmethylated, second->n_methylated);
    } else {
        fprintf(of, "\t-\t0\t0\t\n");
    }
}

FILE * generate_output_name(char *iname) {
    FILE *of = NULL;
    char *p;
    char *oname = malloc(sizeof(char) * (strlen(iname) + 8));
    strcpy(oname, iname);
    p = strrchr(oname, '.');
    if(strcmp(p, ".bedGraph") == 0 || strcmp(p, ".bedgraph") == 0) {
        *p = '\0';
    } else {
        oname = realloc(oname, sizeof(char) * (strlen(oname) + strlen(".moabs")));
    }
    sprintf(oname, "%s.moabs", oname);

    printf("Output will be written to %s\n", oname);
    of = fopen(oname, "w");
    free(oname);
    return of;
}

//Given a chromosome name, return the numeric index of its placement in the chromsomes struct
//It would make sense to memoize the result to prevent continuous lookup
inline int char2tid(char *chrom) {
    int i;
    for(i=0; i<chromosomes.nchromosomes; i++) {
        if(strcmp(chromosomes.chromosome[i]->chrom, chrom) == 0) return i;
    }
    return chromosomes.nchromosomes;
}

inline void process_line(char *line, struct CpG *current_line) {
    char *col;

    //chrom
    col = strtok(line, "\t");
    current_line->tid = char2tid(col);

    //start
    col = strtok(NULL, "\t");
    current_line->start = (int32_t) atoi(col);

    //end
    col = strtok(NULL, "\t");
    current_line->end = (int32_t) atoi(col);

    //1000*methylation percentage
    col = strtok(NULL, "\t");

    //n_methylated
    col = strtok(NULL, "\t");
    current_line->n_methylated = (int32_t) atoi(col);

    //n_unmethylated
    col = strtok(NULL, "\t");
    current_line->n_unmethylated = (int32_t) atoi(col);
}

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] genome_directory/ file.bedGraph [file2.beGraph file3.bedGraph ...]\n", prog);
    printf("\n\
    Convert a bedGraph file to the format required by MOABS. Ideally, the\n\
    metrics in the file should be single-C, not merged into CpGs, but both are\n\
    supported. If metrics have been merged to be per-CpG, then the output\n\
    metrics will always be on the + strand.\n\
\n\
    The output file(s) will be named file.moabs if the input is named\n\
    file.bedGraph. The format of this file is tab-seperated table with the\n\
    following columns:\n\
    #chrom	chromosome name\n\
    start	start coordinate of the CpG\n\
    end		end coordinate of the CpG\n\
    ratio	100*percent methylation\n\
    totalC	Total coverage of a CpG\n\
    methC	Total number methylated bases\n\
    strand	'B': calls on both strands\n\
    		'+': calls only on plus strand\n\
    		'-': calls only on minus strand\n\
    next	'G', denoting a CpG\n\
    Plus	'+', this seems to serve no purpose\n\
    totalC	As with previous totalC, but only for plus strand\n\
    methC	As with previous methC, but only for plus strand\n\
    Minus	'-'\n\
    totalC	As above but for minus strand\n\
    methC	As above but for minus strand\n\
    localSeq	Blank.\n\
\n\
    -h      Print this message.\n\
\n\
    --genome-size Many of the bison tools need to read the genome into memory. By\n\
            default, they allocate 3000000000 bases worth of memory for this and\n\
            increase that as needed. However, this can sometimes be far more\n\
            than is needed (meaning wasted memory) or far too little (in which\n\
            case the process can become quite slow). If you input the\n\
            approximate size of your genome here (in bases), then you can\n\
            maximize performance and minimize wasted space. It's convenient to\n\
            round up a little.\n\
\n");
}

int main(int argc, char *argv[]) {
    int i, nfiles = 0;
    char *fname = NULL, *line = malloc(sizeof(char) * MAXREAD);
    char **fnames = NULL;
    char base;
    unsigned long long offset;
    FILE *of, *ifile;
    struct CpG current_line, last_line;

    config.genome_dir = NULL;
    chromosomes.max_genome = 3000000000;
    chromosomes.nchromosomes = 0;
    last_line.tid = -1; //This will denote "not set"

    /* read in the file names */
    if(argc < 3) {
        usage(argv[0]);
        return 0;
    };
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(config.genome_dir == NULL) {
            config.genome_dir = strdup(argv[i]);
            if(*(config.genome_dir+strlen(config.genome_dir)-1) != '/') {
                config.genome_dir = realloc(config.genome_dir, sizeof(char) * (strlen(config.genome_dir)+2));
                sprintf(config.genome_dir, "%s/", config.genome_dir);
            }
        } else if(strcmp(argv[i], "--genome-size") == 0) {
            i++;
            chromosomes.max_genome = strtoull(argv[i], NULL, 10);
        } else {
            nfiles++;
            fnames = realloc(fnames, sizeof(char*) * nfiles);
            fnames[nfiles-1] = argv[i];
        }
    }

    if(config.genome_dir == NULL || fnames == NULL) {
        printf("Genome directory or bedGraph input file(s) not specified!\n");
        usage(argv[0]);
        return -1;
    }

    //Read in the genome
    printf("Allocating space for %llu characters\n", chromosomes.max_genome); fflush(stdout);
    chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
    *chromosomes.genome = '\0';
    if(chromosomes.genome == NULL) {
        printf("Could not allocate enough room to hold the genome!\n");
        return -1;
    }
    read_genome();

    //per-file
    for(i=0; i<nfiles; i++) {
        fname = fnames[i];

        //Generate the output names and open the output files
        of = generate_output_name(fname);

        ifile = fopen(fname, "r");
        fprintf(of, "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC\tlocalSeq\n");
        while(fgets(line, MAXREAD, ifile) != NULL) {
            if(strncmp("track", line, 5) == 0) continue; //Skip the track line

            process_line(line, &current_line);

            //If this describe a CpG rather than a single C, just process it
            if(current_line.end-current_line.start == 2) {
                print_CpG_line(of, &current_line);
                continue;
            }

            //Compare the current and the last line
            if(last_line.tid != -1) {
                if(last_line.tid == current_line.tid && last_line.end == current_line.start) {
                    //Ensure that last_line describe a 'C'
                    offset = chromosomes.chromosome[last_line.tid]->offset;
                    base = toupper(*(chromosomes.genome+offset+last_line.start));
                    if(base == 'C') {
                        print_line(of, &last_line, &current_line);
                        last_line.tid = -1;
                    } else {
                        print_line(of, NULL, &last_line);
                        copy_CpG(&last_line, &current_line);
                    }
                } else {
                    offset = chromosomes.chromosome[last_line.tid]->offset;
                    base = toupper(*(chromosomes.genome+offset+last_line.start));
                    if(base == 'C') {
                        print_line(of, &last_line, NULL);
                    } else {
                        print_line(of, NULL, &last_line);
                    }
                    copy_CpG(&last_line, &current_line);
                }
            } else {
                copy_CpG(&last_line, &current_line);
            }
        }
        //Deal with any remaining calls
        if(last_line.tid != -1) {
            offset = chromosomes.chromosome[last_line.tid]->offset;
            base = toupper(*(chromosomes.genome+offset+last_line.start));
            if(base == 'C') {
                print_line(of, &last_line, NULL);
            } else {
                print_line(of, NULL, &last_line);
            }
        }

        //Close things for this file
        fclose(of);
        fclose(ifile);
    }

    //Close things up
    if(config.genome_dir != NULL) free(config.genome_dir);
    free(line);
    free(chromosomes.genome);
    for(i=0; i<chromosomes.nchromosomes; i++) {
        free((chromosomes.chromosome[i])->chrom);
        free(*(chromosomes.chromosome+i));
    }
    free(chromosomes.chromosome);
    free(fnames);

    return 0;
};
