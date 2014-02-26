#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAXLINE 512
#define MAXChromosome 400000000

void usage(char *prog_name) {
    printf("Usage: %s [options] GENOME.FA OUTPUT.FA\n",prog_name);
    printf("\t-n X	Maximum number of bases in each read (prior to CG/TG/etc. trimming)\n\t\tDefault is 36. N.B. 10%% more is used to the output.\n");
    printf("\t-TaqI	Create a reduced representation genome that was cut by TaqI as well as MspI.\n");
    printf("\t-h	Print this message\n");
    return;
}

void output_fragment(FILE *of, char *fragment) {
    fprintf(of,"%s",fragment);
    return;
}

unsigned long get_left_mask(char *f, int read_size) {
    unsigned long output = 0, i = 0, max_len = strlen(f)-1;

    while(i < read_size) {
        if(*(f+output) != '\n') {
            i++;
        }
        output++;
        if(output >= max_len) {
            output = max_len;
            break;
        }
    }

    return output;
}

unsigned long get_right_mask(char *f, int read_size) {
    unsigned long output = strlen(f)-1, i = 0;

    while(i < read_size) {
        if(*(f+output) != '\n') {
            i++;
        }
        output--;
        if(output <= 0) {
            output = 0;
            break;
        }
    }

    return output;
}

void process_fragment(FILE *of, char *fragment, int read_size) {
    char *fp = fragment;
    unsigned long i = 0, left_mask, right_mask;

    //Determine the masking coordinates
    if(strlen(fragment)-1 <= read_size) {
        output_fragment(of, fragment);
        return;
    }
    left_mask = get_left_mask(fragment, read_size);
    right_mask = get_right_mask(fragment, read_size);

    if(left_mask < right_mask) {
        for(i=0; i <= right_mask; i++) {
            if(i>=left_mask) {
                if(*fp != '\n') {
                    *fp = 'N';
                }
            }
            fp++;
        }
    }
    output_fragment(of, fragment);

    return;
}

void process_chromosome(FILE *of, char *chrom, int read_size, int Taq) {
    char *fragment = malloc(MAXChromosome * sizeof(char));
    char *cp = chrom, *fp = fragment;

    while(*cp != '\0') {
        //Are we at a TaqI site? Do we even care?
        if(*cp == 'T' && Taq) {
            if(strncmp(cp,"TCGA", 4) == 0) {
                //Add on the last T and a null
                *fp = 'T';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'A';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 4;
                continue;
            } else if(strncmp(cp,"T\nCGA", 5) == 0) {
                //Add on the last T and a null
                *fp = 'T';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = '\n';
                fp++;
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'A';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            } else if(strncmp(cp,"TC\nGA", 5) == 0) {
                //Add on the last T and a null
                *fp = 'T';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = '\n';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'A';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            } else if(strncmp(cp,"TCG\nA", 5) == 0) {
                //Add on the last T and a null
                *fp = 'T';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = '\n';
                fp++;
                *fp = 'A';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            }
        } else if(*cp == 'C') { //MspI site
            if(strncmp(cp,"CCGG", 4) == 0) {
                //Add on the last T and a null
                *fp = 'C';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'G';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 4;
                continue;
            } else if(strncmp(cp,"C\nCGG", 5) == 0) {
                //Add on the last T and a null
                *fp = 'C';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = '\n';
                fp++;
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'G';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            } else if(strncmp(cp,"CC\nGG", 5) == 0) {
                //Add on the last T and a null
                *fp = 'C';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = '\n';
                fp++;
                *fp = 'G';
                fp++;
                *fp = 'G';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            } else if(strncmp(cp,"CCG\nG", 5) == 0) {
                //Add on the last T and a null
                *fp = 'C';
                *(++fp) = '\0';
                process_fragment(of,fragment,read_size);

                fp = fragment; //Move the pointer back to the front of the fragment
                //Start the next fragment
                *fp = 'C';
                fp++;
                *fp = 'G';
                fp++;
                *fp = '\n';
                fp++;
                *fp = 'G';
                fp++;

                //Move the chromosome pointer past the cut site and loop
                cp += 5;
                continue;
            }
        }

        //We are not at a cut site
        *fp = *cp;
        fp++;
        cp++;
    }

    //Don't forget the last fragment!
    *fp = '\0';
    process_fragment(of,fragment,read_size);

    return;
}

int main(int argc, char *argv[]) {
    int Taq = 0;
    int i = 1;
    int read_size = 36; //Set by -n, bp maximum in each read and, therefore number of bases on each end of a fragment to print.
    char *infile = NULL, *outfile = NULL;
    FILE *f = NULL, *of = NULL;
    char *chrom_sequence = malloc(MAXChromosome * sizeof(char)), *line = malloc(MAXLINE * sizeof(char));
    char *p = chrom_sequence;

    if(argc < 3) {
        usage(argv[0]);
        return 1;
    }

    //Parse the input
    while(i<argc) {
        if(strcmp(argv[i],"-h") == 0) {
            usage(argv[0]);
            return 1;
        } else if(strcmp(argv[i],"-n") == 0) {
            read_size = atoi(argv[i+1]);
            i++;
            printf("Changing from default read size to %ibp.\n", read_size);
        } else if(strcmp(argv[i], "-TaqI") == 0) {
            Taq = 1;
            printf("TaqI sites will also be cut.\n");
        } else {
            if(infile == NULL) {
                infile = argv[i];
            } else if(outfile == NULL) {
                outfile = argv[i];
            } else {
                usage(argv[0]);
                return 1;
            }
        };
        i++;
    }
    read_size *= 1.1;
    printf("The first and last %ibp of each fragment will not be masked\n",read_size);

    //Open files for I/O, we should really check if they exist
    f = fopen(infile,"r");
    of = fopen(outfile,"w");

    //Read in the chromomes
    while(fgets(line, MAXLINE, f) != NULL) {
        //Have we switched chromosomes?
        if(*line == '>') {
            if(chrom_sequence != p) {
                *p = '\0'; //Ensure that we end in a null
                process_chromosome(of, chrom_sequence, read_size, Taq);
            }
            fprintf(of,"%s",line);
            p = chrom_sequence;
        } else {
            for(i = 0; i<strlen(line); i++) {
                if(line[i] == 'T' || line[i] == 't') {
                    *p = 'T';
                } else if(line[i] == 'G' || line[i] == 'g') {
                    *p = 'G';
                } else if(line[i] == 'C' || line[i] == 'c') {
                    *p = 'C';
                } else if(line[i] == 'A' || line[i] == 'a') {
                    *p = 'A';
                } else if(line[i] == 'N' || line[i] == 'n') {
                    *p = 'N';
                } else if(line[i] == '\n') {
                    *p = '\n';
                } else {
                    printf("Uhoh, found %c\n",*p);
                }
                p++;
            }
        }
    }

    //Deal with the last chromosome
    *p = '\0'; //Ensure that we end in a null
    process_chromosome(of, chrom_sequence, read_size, Taq);

    free(chrom_sequence);
    fclose(f);
    fclose(of);
    return 0;
}
