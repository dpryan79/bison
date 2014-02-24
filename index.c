#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>

#define MAXLINE 1024

typedef struct {
    char *options;
    char *odir;
    char conversion;
} btoptions_struct;

void usage(char *prog) {
    printf("Usage: %s [options] reference(s)\n", prog);
    printf("\n \
Note, references is a comma separated list of FASTA files. A \"bisulfite_genome\"\n \
directory with CT_conversion and GA_conversion subdirectories will be created.\n \
While the directory structure and indexing method are identical to bismark, the\n \
resulting indexes are not compatible, owing to bismark's changing of\n \
chromosome/contig names.\n \
\n \
Options are currently identical to those for bowtie2-build\n \
(http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer),\n \
as this program is effectively just a wrapper.\n");
}

void * bt2_build(void *a) {
    btoptions_struct *options = (btoptions_struct *) a;
    int rv;
    char *cmd;

    //Create the command
    cmd = malloc(sizeof(char) * (strlen(options->options) + 2*strlen(options->odir) + 200));
    if(options->conversion == 'C') {
        sprintf(cmd, "bowtie2-build %s %s/genome.fa %s/BS_CT", options->options, options->odir, options->odir);
    } else {
        sprintf(cmd, "bowtie2-build %s %s/genome.fa %s/BS_GA", options->options, options->odir, options->odir);
    }
    printf("Now executing: %s\n", cmd);
    rv = system(cmd);
    if(rv) printf("%s returned with status %i!\n", cmd, rv);
    return NULL;
}

int main(int argc, char *argv[]) {
    char *odir, *p, *CT_dir, *GA_dir;
    char *CT_line, *GA_line;
    char *options;
    FILE *fp, *CT, *GA;
    btoptions_struct CT_data, GA_data;
    pthread_t threads[2];
    int i;

    if(argc == 1) {
        usage(argv[0]);
        return 0;
    } else if(strcmp(argv[1], "-h") == 0) {
        usage(argv[0]);
        return 0;
    }

    //Store bowtie2-build options
    options = (char *) calloc(1, sizeof(char));
    for(i=1; i<argc-1; i++) {
        options = realloc(options, sizeof(char) * (strlen(options) + strlen(argv[i]) + 2));
        sprintf(options, "%s %s", options, argv[i]);
    }

    //Create the basename for dir
    odir = strdup(argv[argc-1]);
    p = strchr(odir, ',');
    if(p != NULL) *p = '\0';
    p = strrchr(odir, '/');
    if(p != NULL) {
        p++;
        *p = '\0';
    } else {
        *odir = '\0';
    }

    //Make the output directories
    odir = realloc(odir, sizeof(char) * (strlen(odir) + 1 + strlen("bisulfite_genome")));
    odir = strcat(odir, "bisulfite_genome");
    mkdir(odir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    printf("Output will be placed under %s\n", odir);
    CT_dir = malloc(sizeof(char) * (strlen(odir) + strlen("/CT_conversion/genome.fa") +1));
    sprintf(CT_dir, "%s/CT_conversion", odir);
    mkdir(CT_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    GA_dir = malloc(sizeof(char) * (strlen(odir) + strlen("/CT_conversion/genome.fa") +1));
    sprintf(GA_dir, "%s/GA_conversion", odir);
    mkdir(GA_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    //Iterate through the files, converting them
    CT_dir = strcat(CT_dir, "/genome.fa");
    GA_dir = strcat(GA_dir, "/genome.fa");
    CT = fopen(CT_dir, "w");
    GA = fopen(GA_dir, "w");
    CT_line = malloc(sizeof(char)*MAXLINE);
    GA_line = malloc(sizeof(char)*MAXLINE);
    p = strtok(argv[argc-1], ",");
    while(p != NULL) {
        printf("Reading in and converting %s...\n", p);
        fp = fopen(p, "r");
        while(fgets(CT_line, MAXLINE, fp)) {
            if(*CT_line != '>') {
                for(i=0; i<strlen(CT_line); i++) *(CT_line+i) = toupper(*(CT_line+i));
                strcpy(GA_line, CT_line);
                for(i=0; i<strlen(CT_line); i++) {
                    if(*(CT_line+i) == 'C') *(CT_line+i) = 'T';
                    if(*(GA_line+i) == 'G') *(GA_line+i) = 'A';
                }
            } else {
                strcpy(GA_line, CT_line);
            }
            fputs(CT_line, CT);
            fputs(GA_line, GA);
        }
        fclose(fp);
        p = strtok(NULL, ",");
    }
    fclose(CT);
    fclose(GA);

    //Invoke bowtie2-build in 2 threads
    p = strrchr(CT_dir, '/');
    *p = '\0';
    p = strrchr(GA_dir, '/');
    *p = '\0';
    CT_data.odir = CT_dir;
    CT_data.options = options;
    CT_data.conversion = 'C';
    GA_data.odir = GA_dir;
    GA_data.options = options;
    GA_data.conversion = 'G';
    pthread_create(&threads[0], NULL, &bt2_build, (void *) &CT_data);
    pthread_create(&threads[1], NULL, &bt2_build, (void *) &GA_data);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);

    //We don't actually need these anymore
    *(CT_dir+strlen(CT_dir)) = '/';
    *(GA_dir+strlen(GA_dir)) = '/';
    printf("Removing %s\n", CT_dir);
    printf("Removing %s\n", GA_dir);
    unlink(CT_dir);
    unlink(GA_dir);

    //Cleaning up
    free(options);
    free(odir);
    free(CT_dir);
    free(GA_dir);
    return 0;
}
