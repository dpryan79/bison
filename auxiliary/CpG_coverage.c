#include "../bison.h"
#include "sam.h"

//This will hold the coverage. The last bin is actually for everything >250
unsigned long long coverage[252];
struct {
    char *chrom;
    unsigned long long position;
    unsigned long long end;
    unsigned long coverage;
} cur_line;

void next_line(FILE *fp, char *buffer) {
    if(fgets(buffer, 1024, fp) != NULL) {
        cur_line.chrom = strtok(buffer, "\t");
        cur_line.position = strtoull(strtok(NULL, "\t"), NULL, 10);
        cur_line.end = strtoull(strtok(NULL, "\t"), NULL, 10);
        strtok(NULL, "\t");
        cur_line.coverage = strtoul(strtok(NULL, "\t"), NULL, 10);
        cur_line.coverage += strtoul(strtok(NULL, "\n"), NULL, 10);
    }
}

void usage(char *prog) {
    printf("Usage: %s genome_directory input.bedGraph output.txt\n", prog);
    printf("\n\
    Calculate a histogram of per-CpG coverage. N.B., the genome and bedGraph\n\
    file need to be in the same order (they will be if the bedGraph file was\n\
    produced with bison and the same genome is used).\n\
\n\
    -h            Print this message.\n\
\n");
}

int main(int argc, char *argv[]) {
    FILE *fp = NULL;
    FILE *ofile = NULL;
    int32_t i = 0;
    uint32_t j = 0;
    char *GenomeChrom = NULL;
    unsigned long temp_coverage = 0;
    char *line = malloc(sizeof(char) * 1024);
    unsigned long long k;
    unsigned long long nCpGs = 0;

    config.genome_dir = NULL;
    chromosomes.nchromosomes = 0;

    /* read in the file names */
    if(argc < 4 || strcmp(argv[1], "-h") == 0) {
        usage(argv[0]);
        return 1;
    };
    config.genome_dir = argv[1];
    fp = fopen(argv[2], "r");
    ofile = fopen(argv[3], "w");

    for(i=0; i<252; i++) coverage[i] = 0;

    //Read in the genome
    chromosomes.max_genome = 3000000000;
    printf("Allocating space for %llu characters\n", chromosomes.max_genome); fflush(stdout);
    chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
    *chromosomes.genome = '\0';
    if(chromosomes.genome == NULL) {
        printf("Could not allocate enough room to hold the genome!\n");
        return -1;
    }
    read_genome();

    //Start reading in the file
    next_line(fp, line);

    //Iterate through the genome
    for(i=0; i<chromosomes.nchromosomes; i++) {
        GenomeChrom = chromosomes.chromosome[i]->chrom;
        j = chromosomes.chromosome[i]->offset;
        k = 0; //0-based chromosome position
        while(j < chromosomes.chromosome[i]->length - 1) {
            if(*(chromosomes.genome+j) == 'C' && *(chromosomes.genome+j+1) == 'G') {
                nCpGs++;
                while(strcmp(cur_line.chrom, GenomeChrom) == 0 && k > cur_line.position) next_line(fp, line); //We should never go beyond 1 line...
                if(strcmp(cur_line.chrom, GenomeChrom) == 0 && (k == cur_line.position || k == cur_line.position-1)) {
                    temp_coverage = cur_line.coverage;
                    if(cur_line.end-cur_line.position == 1) { //Single-C resolution rather than merged as CpGs
                        next_line(fp, line);
                        if(strcmp(cur_line.chrom, GenomeChrom) == 0 && k == cur_line.position-1) {
                            temp_coverage += cur_line.coverage;
                        }
                    }
                    if(temp_coverage > 250) temp_coverage = 251;
                    coverage[temp_coverage]++;
                } else {
                    coverage[0]++;
                }
            }
            j++;
            k++;
        }
    }

    //Print some output
    for(i=0; i<251; i++) fprintf(ofile, "%i\t%llu\n", i, coverage[i]);
    fprintf(ofile, "251+\t%llu\n", coverage[251]);
    printf("There were %llu CpGs\n", nCpGs);

    //Close things up
    free(line);
    free(chromosomes.genome);
    for(i=0; i<chromosomes.nchromosomes; i++) {
        free((chromosomes.chromosome[i])->chrom);
        free(*(chromosomes.chromosome+i));
    }
    free(chromosomes.chromosome);
    fclose(fp);
    fclose(ofile);

    return 0;
};
