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
    printf("Usage: %s genome_directory file.bedGraph\n", prog);
    printf("\n\
    Convert a CpG bedGraph file to the format required for methylKit.\n\
    The CpGs in the file should not be merged (i.e., they should represent\n\
    individual strand)! See the methylKit documentation for the file format.\n\
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

FILE * generate_output_name(char *iname) {
    FILE *of = NULL;
    char *p;
    char *oname = malloc(sizeof(char) * (strlen(iname) + 8));
    strcpy(oname, iname);
    p = strrchr(oname, '.');
    if(strcmp(p, ".bedGraph") == 0 || strcmp(p, ".bedgraph") == 0) {
        *p = '\0';
    } else {
        oname = realloc(oname, sizeof(char) * (strlen(oname) + strlen(".methylKit")));
    }
    sprintf(oname, "%s.methylKit", oname);

    printf("Output will be written to %s\n", oname);
    of = fopen(oname, "w");
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
    int i, last_tid = 0;
    char *fname = NULL, *line = malloc(sizeof(char) * MAXREAD);
    char *chrom, *last_chrom = NULL;
    char base, strand;
    unsigned long long offset;
    FILE *of, *ifile;
    struct CpG current_line;

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
        } else if(config.genome_dir == NULL) {
            config.genome_dir = argv[i];
        } else if(fname == NULL) {
            fname = argv[i];
        } else if(strcmp(argv[i], "--genome-size") == 0) {
            i++;
            chromosomes.max_genome = strtoull(argv[i], NULL, 10);
        } else {
            printf("Got an unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if(config.genome_dir == NULL || fname == NULL) {
        printf("Genome directory or SAM/BAM input file not specified!\n");
        usage(argv[0]);
    }

    //Generate the output names and open the output files
    of = generate_output_name(fname);

    //Read in the genome
    printf("Allocating space for %llu characters\n", chromosomes.max_genome); fflush(stdout);
    chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
    *chromosomes.genome = '\0';
    if(chromosomes.genome == NULL) {
        printf("Could not allocate enough room to hold the genome!\n");
        return -1;
    }
    read_genome();

    ifile = fopen(fname, "r");
    fprintf(of, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    while(fgets(line, MAXREAD, ifile) != NULL) {
        chrom = strtok(line, "\t");
        if(last_chrom == NULL || strcmp(chrom, last_chrom) != 0) {
            last_tid = char2tid(chrom);
            if(last_chrom != NULL) free(last_chrom);
            last_chrom = strdup(chrom);
        }
        current_line.tid = last_tid;
        process_line(line, &current_line);

        //Determine the strand
        offset = chromosomes.chromosome[last_tid]->offset;
        base = toupper(*(chromosomes.genome+offset+current_line.start));
        strand='R';
        if(base=='C') strand='F';
        fprintf(of, "%s.%i\t%s\t%i\t%c\t%i\t%5.2f\t%5.2f\n", chrom, current_line.start+1, chrom, current_line.start+1, strand, \
            current_line.n_methylated+current_line.n_unmethylated, \
            100*((float) current_line.n_methylated)/(float)(current_line.n_methylated + current_line.n_unmethylated), \
            100*((float) current_line.n_unmethylated)/(float)(current_line.n_methylated + current_line.n_unmethylated));
    }

    //Close things up
    free(line);
    fclose(of);
    fclose(ifile);
    free(chromosomes.genome);
    for(i=0; i<chromosomes.nchromosomes; i++) {
        free((chromosomes.chromosome[i])->chrom);
        free(*(chromosomes.chromosome+i));
    }
    free(chromosomes.chromosome);

    return 0;
};
