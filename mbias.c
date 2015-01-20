#include "bison.h"
#include <math.h>

unsigned long long *r1_m[4];
unsigned long long *r1_um[4];
unsigned long long *r2_m[4];
unsigned long long *r2_um[4];
int min_phred = 5;

void store_calls(unsigned long long *m, unsigned long long *um, bam1_t *read, int reversed) {
    char *meth = bam_aux2Z(bam_aux_get(read, "XM"));
    uint8_t *qual = bam_get_qual(read);
    int i;

    if(!reversed) {
        for(i=0; i<strlen(meth); i++) {
            if(*(qual+i) < min_phred) continue;
            if(*(meth+i) == 'Z') *(m+i) += 1;
            if(*(meth+i) == 'z') *(um+i) += 1;
        }
    } else {
        for(i=strlen(meth)-1; i>=0; i--) {
            if(*(qual+i) < min_phred) continue;
            if(*(meth+i) == 'Z') *(m+i) += 1;
            if(*(meth+i) == 'z') *(um+i) += 1;
        }
    }
}

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] file.bam\n", prog);
    printf("\n\
    Compute the methylation percentage as a function of read position for a BAM\n\
    file. The output can be conveniently plotted with the accompanying\n\
    plot_mbias.R script\n\
\n\
    -phred Minimum Phred score that a base must have for inclusion in the\n\
           metrics (default 5).\n\
\n\
    -q     Read MAPQ value must at least this for inclusion (default 10).\n\
           Specify 0 to include everything.\n\
\n\
    -pdf   Run the R script to convert the output to pdf format, including\n\
           recommended inclusion bounds. R must be installed and in your PATH.\n");
}

int main(int argc, char *argv[]) {
    htsFile *ifile = NULL;
    FILE *ofile = NULL;
    char *prefix = NULL;
    char *p, *XR, *XG;
    bam1_t *read = bam_init1();
    bam_hdr_t *header = NULL;
    int max_length = 50;
    int paired = 0, reversed = 0, hasComp = 0;
    int i, j, min_mapq = 10, pdf = 0;
    unsigned long long treads = 0;

    if(argc < 2) {
        usage(argv[0]);
        return 1;
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-phred") == 0) {
            min_phred = atoi(argv[++i]);
        } else if(strcmp(argv[i], "-q") == 0) {
            min_mapq = atoi(argv[++i]);
        } else if(strcmp(argv[i], "-pdf") == 0) {
            pdf = 1;
        } else if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(prefix == NULL) {
            ifile = sam_open(argv[i], "rb");
            header = sam_hdr_read(ifile);
            prefix = strdup(argv[i]);
        } else {
            fprintf(stderr, "Unknown option: %s", argv[i]);
            free(prefix);
            bam_hdr_destroy(header);
            sam_close(ifile);
            usage(argv[0]);
            return -1;
        }
    }

    if(prefix == NULL) {
        fprintf(stderr, "No BAM file specified!\n");
        usage(argv[0]);
        return -1;
    }

    //Create enough space to hold the metrics
    for(i=0; i<4; i++) {
        r1_m[i] = calloc(max_length, sizeof(unsigned long long));
        r1_um[i] = calloc(max_length, sizeof(unsigned long long));
        r2_m[i] = calloc(max_length, sizeof(unsigned long long));
        r2_um[i] = calloc(max_length, sizeof(unsigned long long));
    }

    //Create the output name
    p = strrchr(prefix, '.');
    *p = '\0';
    prefix = realloc(prefix, sizeof(char) * (strlen(prefix) + strlen("_mbias.txt ")));
    sprintf(prefix, "%s_mbias.txt", prefix);
    ofile = fopen(prefix, "w");

    while(sam_read1(ifile, header, read) > 1) {
        if(++treads % 10000000 == 0) fprintf(stderr, "Processed %llu reads\n", treads);
        if(read->core.qual < min_mapq) continue;
        if(read->core.flag & BAM_FUNMAP) continue;

        //Lengthen the output arrays if needed
        if(read->core.l_qseq > max_length) {
            for(i=0; i<4; i++) {
                r1_m[i] = realloc(r1_m[i], read->core.l_qseq * sizeof(unsigned long long));
                r1_um[i] = realloc(r1_um[i], read->core.l_qseq * sizeof(unsigned long long));
                r2_m[i] = realloc(r2_m[i], read->core.l_qseq * sizeof(unsigned long long));
                r2_um[i] = realloc(r2_um[i], read->core.l_qseq * sizeof(unsigned long long));
                for(j=max_length; j<read->core.l_qseq; j++) {
                    *(r1_m[i]+j) = 0;
                    *(r1_um[i]+j) = 0;
                    *(r2_m[i]+j) = 0;
                    *(r2_um[i]+j) = 0;
                }
            }
            max_length = read->core.l_qseq;
        }
        reversed = (read->core.flag & BAM_FREVERSE) ? 1 : 0;

        if(bam_aux_get(read, "XR") == NULL || bam_aux_get(read, "XG") == NULL) { 
            kstring_t str;
            str.l = str.m = 0; str.s = NULL;
            sam_format1(header, read, &str);
            fprintf(stderr, "Error: Couldn't retrieve strand information for %s\n", str.s);
            free(str.s);
        }
        XR = bam_aux2Z(bam_aux_get(read, "XR"));
        XG = bam_aux2Z(bam_aux_get(read, "XG"));
        if(!(read->core.flag & BAM_FREAD2)) {
            if(strcmp(XG, "CT") == 0) { //OT or CTOT
                if(strcmp(XR, "CT") == 0) { //OT
                    store_calls(r1_m[0], r1_um[0], read, reversed);
                } else { //CTOT
                    hasComp = 1;
                    store_calls(r1_m[1], r1_um[1], read, reversed);
                }
            } else {
                if(strcmp(XR, "CT") == 0) { //OB
                    store_calls(r1_m[2], r1_um[2], read, reversed);
                } else { //CTOB
                    hasComp = 1;
                    store_calls(r1_m[3], r1_um[3], read, reversed);
                }
            }
        } else {
            paired = 1;
            if(strcmp(XG, "CT") == 0) { //OT or CTOT
                if(strcmp(XR, "GA") == 0) { //OT
                    store_calls(r2_m[0], r2_um[0], read, reversed);
                } else { //CTOT
                    hasComp = 1;
                    store_calls(r2_m[1], r2_um[1], read, reversed);
                }
            } else {
                if(strcmp(XR, "GA") == 0) { //OB
                    store_calls(r2_m[2], r2_um[2], read, reversed);
                } else { //CTOB
                    hasComp = 1;
                    store_calls(r2_m[3], r2_um[3], read, reversed);
                }
            }
        }
    }

    //Output the calls
    fprintf(ofile, "Strand\tRead\tPosition\tnMethylated\tnUnmethylated\n");
    for(i=0; i<max_length; i++) { //OT
        if(r1_m[0][i] > 0 || r1_um[0][i] > 0) fprintf(ofile, "OT\t1\t%i\t%llu\t%llu\n", i+1, r1_m[0][i], r1_um[0][i]);
        if(paired) {
            if(r2_m[0][i] > 0 || r2_um[0][i] > 0) fprintf(ofile, "OT\t2\t%i\t%llu\t%llu\n", i+1, r2_m[0][i], r2_um[0][i]);
        }
    }
    for(i=0; i<max_length; i++) { //OB
        if(r1_m[2][i] > 0 || r1_um[2][i] > 0) fprintf(ofile, "OB\t1\t%i\t%llu\t%llu\n", i+1, r1_m[2][i], r1_um[2][i]);
        if(paired) {
            if(r2_m[2][i] > 0 || r2_um[2][i] > 0) fprintf(ofile, "OB\t2\t%i\t%llu\t%llu\n", i+1, r2_m[2][i], r2_um[2][i]);
        }
    }
    if(hasComp) {
        for(i=0; i<max_length; i++) { //CTOT
            if(r1_m[1][i] > 0 || r1_um[1][i] > 0) fprintf(ofile, "CTOT\t1\t%i\t%llu\t%llu\n", i+1, r1_m[1][i], r1_um[1][i]);
            if(paired) {
                if(r2_m[1][i] > 0 || r2_um[1][i] > 0) fprintf(ofile, "CTOT\t2\t%i\t%llu\t%llu\n", i+1, r2_m[1][i], r2_um[1][i]);
            }
        }
        for(i=0; i<max_length; i++) { //CTOB
            if(r1_m[3][i] > 0 || r1_um[3][i] > 0) fprintf(ofile, "CTOB\t1\t%i\t%llu\t%llu\n", i+1, r1_m[3][i], r1_um[3][i]);
            if(paired) {
                if(r2_m[3][i] > 0 || r2_um[3][i] > 0) fprintf(ofile, "CTOB\t2\t%i\t%llu\t%llu\n", i+1, r2_m[3][i], r2_um[3][i]);
            }
        }
    }

    fprintf(stderr, "Processed %llu reads\n", treads);
    fclose(ofile);

    if(pdf) {
        char *cmd = malloc(sizeof(char) * (strlen("bison_mbias2pdf ") + strlen(prefix) + 1));
        sprintf(cmd, "bison_mbias2pdf %s", prefix);
        fprintf(stderr, "Executing %s\n", cmd);
        if(system(cmd) == -1) fprintf(stderr, "N.B. an error occured while running bison_mbias2pdf!\n");
        free(cmd);
    } else {
        fprintf(stderr, "The output may be converted to PDF with recommended inclusion bounds by running bison_mbias2pdf %s\n", prefix);
    }

    //Cleanup
    free(prefix);
    for(i=0; i<4; i++) {
        free(r1_m[i]);
        free(r1_um[i]);
        free(r2_m[i]);
        free(r2_um[i]);
    }
    bam_hdr_destroy(header);
    bam_destroy1(read);
    sam_close(ifile);
    return(0);
}
