#include "bison.h"

char * reverse_qual(char *qual) {
    char *output = malloc(sizeof(char)*(1+strlen(qual)));
    int i, j;
    for(i=0, j=strlen(qual)-1; i<strlen(qual); i++, j--) *(output+i) = *(qual+j);
    *(output+i) = '\0';
    free(qual);
    return output;
}

/******************************************************************************
*
*   Write an unmapped read to a gzipped fastq file.
*
*   FILE *fp: gzipped fastq file
*   bam1_t *read: read to write in fastq format
*
*******************************************************************************/
void write_unmapped(FILE *fp, bam1_t *read) {
    char *seq = calloc(1+read->core.l_qseq, sizeof(char));
    char *qual = calloc(1+read->core.l_qseq, sizeof(char));
    uint8_t b, *seqp = bam1_seq(read), *qualp = bam1_qual(read);
    int i;

    for(i=0; i<read->core.l_qseq; i++) {
        b = bam1_seqi(seqp, i);
        if(b == 1) *(seq+i) = 'A';
        else if(b == 2) *(seq+i) = 'C';
        else if(b == 4) *(seq+i) = 'G';
        else if(b == 8) *(seq+i) = 'T';
        else if(b == 15) *(seq+i) = 'N';
        *(qual+i) = qualp[i] + 33;
    }
    if(read->core.flag & BAM_FREVERSE) {
        reverse_complement(seq);
        qual = reverse_qual(qual);
    }

    fprintf(fp, "@%s\n", bam1_qname(read));
    fprintf(fp, "%s\n", seq);
    fprintf(fp, "+\n");
    fprintf(fp, "%s\n", qual);

    free(seq);
    free(qual);
}

/******************************************************************************
*
*   Construct the output directory name, putting it in config.odir
*
*******************************************************************************/
void update_odir() {
    char *p, *tmp;

    if(config.odir == NULL) {
        tmp = strdup(config.FASTQ1);
        p = strrchr(tmp, '/');
        if(p != NULL) {
            *(p+1) = '\0';
            config.odir = tmp;
        } else {
            config.odir = malloc(sizeof(char) * 3);
            sprintf(config.odir, "./");
        }
    } else {
        if(config.odir[strlen(config.odir)-1] != '/') {
            config.odir = realloc(config.odir, (strlen(config.odir)+2) * sizeof(char));
            strcat(config.odir, "/");
        }
    }
}

/******************************************************************************
*
*   Given the name of a (possibly gzipped) fastq file, return the file name
*   with the .fastq.gz, .fq.gz, .fastq, or .fq extension removed.
*
*   CAUTION, THE OUTPUT MUST BE free()d!
*
*******************************************************************************/
char * get_basename(char *file) {
    char *output = malloc(sizeof(char) * (strlen(file) + 1));
    char *p = NULL;

    //Create the basename of the input
    strcpy(output, file);
    p = strrchr(output, '.');
    if(p != NULL) {
        if(strcmp(p, ".gz") == 0) {
            *p = '\0';
            p = strrchr(output, '.');
            if(p != NULL) {
                if(strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) *p = '\0';
            }
        } else if(strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) {
            *p = '\0';
        }
    }

    //Remove any preceding path
    p = strrchr(output, '/');
    if(p != NULL) {
        p++;
        memmove(output, p, strlen(p)+1);
    }
    return output;
}

/******************************************************************************
*
*   These functions are executed via pthreads to convert the fastq sequences.
*
*******************************************************************************/
void * convert1(void *a) {
    char *cmd = malloc(sizeof(char) * (strlen(config.FASTQ1) + 6));
    char *line1 = malloc(MAXREAD*sizeof(char));
    char *line2 = malloc(MAXREAD*sizeof(char));
    FILE *f, *of1, *of2 = NULL;
    unsigned long long total = 0;
    unsigned int limit = *((unsigned int *) a);
    int i;
    char *p;

    //Determine how we should read in the file
    p = strrchr(config.FASTQ1, '.');
    if(strcmp(p, ".gz") == 0 || strcmp(p, ".GZ") == 0) {
        sprintf(cmd, "zcat %s", config.FASTQ1);
    } else if(strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0) {
        sprintf(cmd, "bzcat %s", config.FASTQ1);
    } else {
        sprintf(cmd, "cat %s", config.FASTQ1);
    }
    f = popen(cmd, "r");

    //CT
    cmd = realloc(cmd,sizeof(char) * (strlen(config.FASTQ1CT) + 8));
    sprintf(cmd, "gzip > %s", config.FASTQ1CT);
    of1 = popen(cmd, "w");

    //GA
    if(!config.directional) {
        sprintf(cmd, "gzip > %s", config.FASTQ1GA);
        of2 = popen(cmd, "w");
    }

    //Iterate through
    while(1) {
        //Read name
        if(fgets(line1, MAXREAD, f) == NULL) break;
        total++;
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);
        //Sequence
        assert(fgets(line1, MAXREAD, f) != NULL);
        if(!config.directional) strcpy(line2, line1);
        for(i=0; i<strlen(line1); i++) {
            if(*(line1+i) == 'C' || *(line1+i) == 'c') *(line1+i) = 'T';
            if(!config.directional) if(*(line2+i) == 'G' || *(line2+i) == 'g') *(line2+i) = 'A';
        }
        fputs(line1, of1);
        if(!config.directional) fputs(line2, of2);
        //QUAL header
        assert(fgets(line1, MAXREAD, f) != NULL);
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);
        //QUAL
        assert(fgets(line1, MAXREAD, f) != NULL);
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);

        if(limit) if(total >= limit) break;
    }

    if(!config.quiet) printf("%s contained %llu reads\n", config.FASTQ1, total);
    pclose(f);
    pclose(of1);
    if(!config.directional) pclose(of2);
    free(cmd);
    free(line1);
    free(line2);

    return NULL;
}
void * convert2(void *a) {
    char *cmd = malloc(sizeof(char) * (strlen(config.FASTQ2) + 6));
    char *line1 = malloc(MAXREAD*sizeof(char));
    char *line2 = malloc(MAXREAD*sizeof(char));
    FILE *f, *of1, *of2 = NULL;
    unsigned long long total = 0;
    unsigned int limit = *((unsigned int *) a);
    int i;
    char *p;

    //Determine how we should read in the file
    p = strrchr(config.FASTQ2, '.');
    if(strcmp(p, ".gz") == 0 || strcmp(p, ".GZ") == 0) {
        sprintf(cmd, "zcat %s", config.FASTQ2);
    } else if(strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0) {
        sprintf(cmd, "bzcat %s", config.FASTQ2);
    } else {
        sprintf(cmd, "cat %s", config.FASTQ2);
    }
    f = popen(cmd, "r");

    //GA
    cmd = realloc(cmd, sizeof(char) * (strlen(config.FASTQ2GA) + 8));
    sprintf(cmd, "gzip > %s", config.FASTQ2GA);
    of1 = popen(cmd, "w");

    //CT
    if(!config.directional) {
        sprintf(cmd, "gzip > %s", config.FASTQ2CT);
        of2 = popen(cmd, "w");
    }

    //Iterate through
    while(1) {
        //Read name
        if(fgets(line1, MAXREAD, f) == NULL) break;
        total++;
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);
        //Sequence
        assert(fgets(line1, MAXREAD, f) != NULL);
        if(!config.directional) strcpy(line2, line1);
        for(i=0; i<strlen(line1); i++) {
            if(*(line1+i) == 'G' || *(line1+i) == 'g') *(line1+i) = 'A';
            if(!config.directional) if(*(line2+i) == 'C' || *(line2+i) == 'c') *(line2+i) = 'T';
        }
        fputs(line1, of1);
        if(!config.directional) fputs(line2, of2);
        //QUAL header
        assert(fgets(line1, MAXREAD, f) != NULL);
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);
        //QUAL
        assert(fgets(line1, MAXREAD, f) != NULL);
        fputs(line1, of1);
        if(!config.directional) fputs(line1, of2);
        if(limit) if(total >= limit) break;
    }

    if(!config.quiet) printf("%s contained %llu reads\n", config.FASTQ2, total);
    pclose(f);
    pclose(of1);
    if(!config.directional) pclose(of2);
    free(cmd);
    free(line1);
    free(line2);

    return NULL;
}

/******************************************************************************
*
*   Invoke the C->T and G->A conversion threads of the fastq files (located in
*   the global config structure).
*
*   FLAGS: integer bit field denoting the conversions to make
*       0x8 fastq #1 C->T
*       0x4 fastq #1 G->A
*       0x2 fastq #2 C->T
*       0x1 fastq #2 G->A
*
*******************************************************************************/
void convert_fastq(int FLAGS, unsigned int limit) {
    pthread_t *threads;
    int rc;

    if(!config.quiet) {
        if(FLAGS & 8) printf("Will C->T convert %s and store the results in %s.\n", config.FASTQ1, config.FASTQ1CT);
        if(FLAGS & 4) printf("Will G->A convert %s and store the results in %s.\n", config.FASTQ1, config.FASTQ1GA);
        if(FLAGS & 2) printf("Will C->T convert %s and store the results in %s.\n", config.FASTQ2, config.FASTQ2CT);
        if(FLAGS & 1) printf("Will G->A convert %s and store the results in %s.\n", config.FASTQ2, config.FASTQ2GA);
    }

    if(config.paired) {
        threads = calloc(2, sizeof(pthread_t));
        rc = pthread_create(&(threads[0]), NULL, convert1, (void *) &limit);
        if(rc) {
            printf("An error occured with invoking pthread_create; %d\n", rc);
            exit(-1);
        }
        rc = pthread_create(&(threads[1]), NULL, convert2, (void *) &limit);
        if(rc) {
            printf("An error occured with invoking pthread_create; %d\n", rc);
            exit(-1);
        }
    } else {
        threads = calloc(1, sizeof(pthread_t));
        rc = pthread_create(&(threads[0]), NULL, convert1, (void *) &limit);
        if(rc) {
            printf("An error occured with invoking pthread_create; %d\n", rc);
            exit(-1);
        }
    }
    pthread_join(threads[0], NULL);
    if(config.paired) pthread_join(threads[1], NULL);

    free(threads);
}

/******************************************************************************
*
*   Take the config.FASTQ1 and config.FASTQ2 filenames and use them to generate
*   the config.FASTQ1CT... filenames. These must subsequently be free()d, which
*   is done in the quit() function.
*
*******************************************************************************/
void create_fastq_names(char *f1, char *f2) {
    char *basename1 = malloc(sizeof(char) * (strlen(f1) + 20));
    char *basename2 = NULL;
    char *p;

    basename1 = strcpy(basename1, f1);
    if(config.paired) {
        basename2 = malloc(sizeof(char) * (strlen(f2) + 20));
        basename2 = strcpy(basename2, f2);
    }

    //Create the basename of FASTQ1, trim off [.fastq/.fq].(bz/gz/bz2/fastq/fq)
    p = strrchr(basename1, '.');
    if(p != NULL) {
        if(strcmp(p, ".gz") == 0 || strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0 || strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) {
            *p = '\0';
            p = strrchr(basename1, '.');
            if(p != NULL) {
                if(strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) *p = '\0';
            }
        }
    }
    config.FASTQ1CT = malloc(sizeof(char) * (strlen(basename1) + 10));
    config.FASTQ1GA = malloc(sizeof(char) * (strlen(basename1) + 10));
    if(config.odir != NULL) {
        p = strrchr(basename1, '/');
        if(p!=NULL) {
            p++;
        } else {
            p = basename1;
        }
        config.unmapped1 = malloc(sizeof(char) * (strlen(config.odir) + strlen(p) + strlen(".unmapped.fq.gz") + 1));
        sprintf(config.unmapped1, "%s%s.unmapped.fq.gz", config.odir, p);
    } else {
        config.unmapped1 = malloc(sizeof(char) * (strlen(basename1) + strlen(".unmapped.fq.gz") + 1));
        sprintf(config.unmapped1, "%s.unmapped.fq.gz", basename1);
    }
    sprintf(config.FASTQ1CT, "%s.CT.fq.gz", basename1);
    sprintf(config.FASTQ1GA, "%s.GA.fq.gz", basename1);

    //Create the basename of FASTQ2
    if(config.paired) {
        p = strrchr(basename2, '.');
        if(p != NULL) {
            if(strcmp(p, ".gz") == 0 || strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0 || strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) {
                *p = '\0';
                p = strrchr(basename2, '.');
                if(p != NULL) {
                    if(strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) *p = '\0';
                }
            }
        }
        config.FASTQ2CT = malloc(sizeof(char) * (strlen(basename2) + 10));
        config.FASTQ2GA = malloc(sizeof(char) * (strlen(basename2) + 10));
        if(config.odir != NULL) {
            p = strrchr(basename2, '/');
            if(p!=NULL) {
                p++;
            } else {
                p = basename2;
            }
            config.unmapped2 = malloc(sizeof(char) * (strlen(config.odir) + strlen(p) + strlen(".unmapped.fq.gz") + 1));
            sprintf(config.unmapped2, "%s%s.unmapped.fq.gz", config.odir, p);
        } else {
            config.unmapped2 = malloc(sizeof(char) * (strlen(basename2) + strlen(".unmapped.fq.gz") + 1));
            sprintf(config.unmapped2, "%s.unmapped.fq.gz", basename2);
        }
        sprintf(config.FASTQ2CT, "%s.CT.fq.gz", basename2);
        sprintf(config.FASTQ2GA, "%s.GA.fq.gz", basename2);
        free(basename2);
    }

    free(basename1);
}
