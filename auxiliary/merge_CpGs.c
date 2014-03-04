#include "../bison.h"
#include "sam.h"

struct CpG {
    int tid;
    int start;
    int end;
    unsigned int n_methylated;
    unsigned int n_unmethylated;
};

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] genome_directory/ file.bedGraph [file2.bedGraph file3.bedGraph]\n", prog);
    printf("\n\
    Merge strand metrics for individual CpG calls (i.e. if there are separate\n\
    methylation metrics for the C's on the + and - strand of a CpG site, combine\n\
    them). If you specify more than one bedGraph file, they will each be merged\n\
    independently and written to separate files.\n\
\n\
    -h      Print this message.\n\
\n\
    --genome-size Many of the bison tools need to read the genome into memory. By\n \
            default, they allocate 3000000000 bases worth of memory for this and\n \
            increase that as needed. However, this can sometimes be far more\n \
            than is needed (meaning wasted memory) or far too little (in which\n \
            case the process can become quite slow). If you input the\n \
            approximate size of your genome here (in bases), then you can\n \
            maximize performance and minimize wasted space. It's convenient to\n \
            round up a little.\n \
\n");
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
        oname = realloc(oname, sizeof(char) * (strlen(oname) + strlen(".merged.bedGraph ")));
    }
    sprintf(oname, "%s.merged.bedGraph", oname);

    printf("Writting %s\n", oname);
    of = fopen(oname, "w");
    fprintf(of, "track type=bedGraph\n");
    free(oname);
    return of;
}

//Given a chromosome name, return the numeric index of it's placement in the chromsomes struct
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

int main(int argc, char *argv[]) {
    int i, fstart = -1, last_tid = 0, mpercent;
    char *line = malloc(sizeof(char) * MAXREAD);
    char *chrom, *last_chrom = NULL;
    char base;
    unsigned long long offset;
    FILE *of, *ifile;
    struct CpG current_line, last_line;

    last_line.tid = -1; //This will mean that the last line has been written
    last_line.start = 0;
    last_line.end = 0;
    last_line.n_methylated = 0;
    last_line.n_unmethylated = 0;
    config.genome_dir = NULL;
    chromosomes.max_genome = 3000000000;
    chromosomes.nchromosomes = 0;

    /* read in the file names */
    if(argc < 3) {
        usage(argv[0]);
        return 0;
    };
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp(argv[i], "--genome-size") == 0) {
            chromosomes.max_genome = strtoull(argv[++i], NULL, 10);
        } else if(config.genome_dir == NULL) {
            config.genome_dir = strdup(argv[i]);
            if(*(config.genome_dir+strlen(config.genome_dir)-1) != '/') {
                config.genome_dir = realloc(config.genome_dir, sizeof(char) * (strlen(config.genome_dir)+2));
                sprintf(config.genome_dir, "%s/", config.genome_dir);
            }
        } else if(fstart == -1) {
            fstart = i;
            break;
        } else {
            printf("Got an unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if(config.genome_dir == NULL || fstart == -1) {
        printf("Genome directory or SAM/BAM input file not specified!\n");
        usage(argv[0]);
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

    for(i=fstart; i<argc; i++) {
        //Generate the output names and open the output files
        ifile = fopen(argv[i], "r");

        of = generate_output_name(argv[i]);
        while(fgets(line, MAXREAD, ifile) != NULL) {
            if(strncmp(line, "track", 5) == 0) if(fgets(line, MAXREAD, ifile) == NULL) break; //Skip the track line
            chrom = strtok(line, "\t");
            if(last_chrom == NULL || strcmp(chrom, last_chrom) != 0) {
                last_tid = char2tid(chrom);
                if(last_chrom != NULL) free(last_chrom);
                last_chrom = strdup(chrom);
            }
            current_line.tid = last_tid;
            process_line(line, &current_line);

            //Compare the current and last calls
            if(current_line.tid == last_line.tid && current_line.start == last_line.end) {
                //Are these different strands of a single CpG?
                offset = chromosomes.chromosome[last_line.tid]->offset;
                base = toupper(*(chromosomes.genome+offset+last_line.start));
                if(base == 'C') { //Yes
                    last_line.end++;
                    last_line.n_methylated += current_line.n_methylated;
                    last_line.n_unmethylated += current_line.n_unmethylated;
                    mpercent = (int) (1000 * ((float) last_line.n_methylated)/(float)(last_line.n_methylated + last_line.n_unmethylated));
                    fprintf(of, "%s\t%i\t%i\t%i\t%i\t%i\n", chromosomes.chromosome[last_tid]->chrom, last_line.start, \
                        last_line.end, mpercent, last_line.n_methylated, last_line.n_unmethylated);
                    last_line.tid = -1;
                } else { //No
                    last_line.start--;
                    mpercent = (int) (1000 * ((float) last_line.n_methylated)/(float)(last_line.n_methylated + last_line.n_unmethylated));
                    fprintf(of, "%s\t%i\t%i\t%i\t%i\t%i\n", chromosomes.chromosome[last_tid]->chrom, last_line.start, \
                        last_line.end, mpercent, last_line.n_methylated, last_line.n_unmethylated);
                    last_line.tid = current_line.tid;
                    last_line.start = current_line.start;
                    last_line.end = current_line.end;
                    last_line.n_methylated = current_line.n_methylated;
                    last_line.n_unmethylated = current_line.n_unmethylated;
                }
            } else {
                if(last_line.tid != -1) {
                    offset = chromosomes.chromosome[last_line.tid]->offset;
                    base = toupper(*(chromosomes.genome+offset+last_line.start));
                    if(base == 'C') { //Yes
                        last_line.end++;
                    } else {
                        last_line.start--;
                    }
                    mpercent = (int) (1000 * ((float) last_line.n_methylated)/(float)(last_line.n_methylated + last_line.n_unmethylated));
                    fprintf(of, "%s\t%i\t%i\t%i\t%i\t%i\n", chromosomes.chromosome[last_line.tid]->chrom, last_line.start, \
                        last_line.end, mpercent, last_line.n_methylated, last_line.n_unmethylated);
                }
                last_line.tid = current_line.tid;
                last_line.start = current_line.start;
                last_line.end = current_line.end;
                last_line.n_methylated = current_line.n_methylated;
                last_line.n_unmethylated = current_line.n_unmethylated;
            }
        }
        //Attend to a possible remnant line
        if(last_line.tid != -1) {
            offset = chromosomes.chromosome[last_tid]->offset;
            base = toupper(*(chromosomes.genome+offset+last_line.start));
            if(base == 'C') { //Yes
                last_line.end++;
            } else {
                last_line.start--;
            }
            mpercent = (int) (1000 * ((float) last_line.n_methylated)/(float)(last_line.n_methylated + last_line.n_unmethylated));
            fprintf(of, "%s\t%i\t%i\t%i\t%i\t%i\n", chromosomes.chromosome[last_tid]->chrom, last_line.start, \
                last_line.end, mpercent, last_line.n_methylated, last_line.n_unmethylated);
        }
        last_line.tid = -1; //Otherwise, if this will be written as the first line of the next file!
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

    return 0;
};
