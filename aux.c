#include "bison.h"

KSTREAM_INIT(gzFile, gzread, 16384)
KHASH_MAP_INIT_STR(ref, uint64_t)
FILE *popen_fd;

struct __tamFile_t {
        gzFile fp;
        kstream_t *ks;
        kstring_t *str;
        uint64_t n_lines;
        int is_first;
};

/******************************************************************************
*
*   Return the number of nodes what will actually be run (as opposed to the
*   number allocated)
*
*******************************************************************************/
#ifdef DEBUG
int effective_nodes() {
    return(8);
}
#else
int effective_nodes() {
    int output, remainder;

    MPI_Comm_size(MPI_COMM_WORLD, &output);
    --output; //Ignore the master node

    if(config.directional) {
        remainder = output % 2;
    } else {
        remainder = output % 4;
    }
    output -= remainder;
    return(output);
}
#endif

/******************************************************************************
*
*   quit, while performing some cleanup
*
*   int FLAG: What to free/close/etc.
*             0x1 things created by create_fastq_names()
*             0x2 things pthreads are closed and bam headers destroyed
*             In addition, the master node will free chromosomes.genome, close
*             the BAM file, and free everything in the chromosomes struct.
*
*   int rv: return value
*
*******************************************************************************/
void quit(int FLAG, int rv) {
    int taskid, i;

    free(config.bowtie2_options);

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if(FLAG & 1) { //FASTQ filenames set
#ifndef DEBUG
        if(taskid == MASTER) {
            if(config.FASTQ1CT != NULL) remove(config.FASTQ1CT);
            if(config.paired && (config.FASTQ2GA != NULL)) remove(config.FASTQ2GA);
            if(!config.directional) {
                if(config.FASTQ1GA != NULL) remove(config.FASTQ1GA);
                if(config.paired && (config.FASTQ2CT != NULL)) remove(config.FASTQ2CT);
            }
        }
#endif
        if(config.FASTQ1CT != NULL) free(config.FASTQ1CT);
        if(config.FASTQ1GA != NULL) free(config.FASTQ1GA);
        if(config.unmapped1 != NULL) free(config.unmapped1);
        if(config.paired) {
            if(config.FASTQ2CT != NULL) free(config.FASTQ2CT);
            if(config.FASTQ2GA != NULL) free(config.FASTQ2GA);
            if(config.unmapped2 != NULL) free(config.unmapped2);
        }
        free(config.basename);
        free(config.outname);
    }

    if(taskid == MASTER) {
        free(chromosomes.genome);
        for(i=0; i<chromosomes.nchromosomes; i++) {
            free((chromosomes.chromosome[i])->chrom);
            free(*(chromosomes.chromosome+i));
        }
        free(chromosomes.chromosome);
        if(FLAG && OUTPUT_BAM) bam_close(OUTPUT_BAM);
    }
    MPI_Finalize();
    if(taskid == MASTER && FLAG > 0) {
#ifdef DEBUG
        if(fp1) bam_close(fp1);
        if(fp2) bam_close(fp2);
        if(!config.directional) {
            if(fp3) bam_close(fp3);
            if(fp4) bam_close(fp4);
        }
#else
        if(config.unmapped) {
            pclose(unmapped1);
            if(config.paired) pclose(unmapped2);
        }
#endif
    }
    exit(rv);
}

void print_metrics() {
    char *of = malloc(sizeof(char) * (strlen(config.odir)+5+strlen(config.basename)));
    FILE *fp;
    unsigned long long m_reads = m_reads_OT + m_reads_OB + m_reads_CTOT + m_reads_CTOB;
    sprintf(of, "%s%s.txt", config.odir, config.basename);
    fp = fopen(of, "w");

    if(!config.quiet) printf("Alignment:\n");
    fprintf(fp,"Alignment:\n");
    if(config.paired) {
        if(!config.quiet) {
            printf("\t%llu total paired-end reads analysed\n", t_reads);
            printf("\t%llu paired-end reads mapped (%6.2f%%).\n", m_reads, ((float) (100*m_reads))/((float) t_reads));
            printf("\n");
        }
        fprintf(fp, "\t%llu total paired-end reads analysed\n", t_reads);
        fprintf(fp, "\t%llu paired-end reads mapped (%6.2f%%).\n", m_reads, ((float) (100*m_reads))/((float) t_reads));
        fprintf(fp, "\n");
    } else {
        if(!config.quiet) {
            printf("\t%llu total reads analysed\n", t_reads);
            printf("\t%llu reads mapped (%6.2f%%).\n", m_reads, ((float) (100*m_reads))/((float) t_reads));
            printf("\n");
        }
        fprintf(fp,"\t%llu total reads analysed\n", t_reads);
        fprintf(fp,"\t%llu reads mapped (%6.2f%%).\n", m_reads, ((float) (100*m_reads))/((float) t_reads));
        fprintf(fp,"\n");
    }
    if(!config.quiet) {
        printf("Number of hits aligning to each of the orientations:\n");
        printf("\t%llu\t%6.2f%%\tOT (original top strand)\n", m_reads_OT, ((float) (100*m_reads_OT))/((float) t_reads));
        printf("\t%llu\t%6.2f%%\tOB (original bottom strand)\n", m_reads_OB, ((float) (100*m_reads_OB))/((float) t_reads));
        if(!config.directional) printf("\t%llu\t%6.2f%%\tCTOT (complementary to the original top strand)\n", m_reads_CTOT, ((float) (100*m_reads_CTOT))/((float) t_reads));
        if(!config.directional) printf("\t%llu\t%6.2f%%\tCTOB (complementary to the original bottom strand)\n", m_reads_CTOB, ((float) (100*m_reads_CTOB))/((float) t_reads));
        printf("\n");
        printf("Cytosine Methylation (N.B., statistics from overlapping mates are added together!):\n");
        printf("\tNumber of C's in a CpG context: %llu\n", t_CpG);
        printf("\tPercentage of methylated C's in a CpG context: %6.2f%%\n", ((float) (100*m_CpG))/((float) t_CpG));
        printf("\tNumber of C's in a CHG context: %llu\n", t_CHG);
        printf("\tPercentage of methylated C's in a CHG context: %6.2f%%\n", ((float) (100*m_CHG))/((float) t_CHG));
        printf("\tNumber of C's in a CHH context: %llu\n", t_CHH);
        printf("\tPercentage of methylated C's in a CHH context: %6.2f%%\n", ((float) (100*m_CHH))/((float) t_CHH));
    }
    fprintf(fp,"Number of hits aligning to each of the orientations:\n");
    fprintf(fp,"\t%llu\t%6.2f%%\tOT (original top strand)\n", m_reads_OT, ((float) (100*m_reads_OT))/((float) t_reads));
    fprintf(fp,"\t%llu\t%6.2f%%\tOB (original bottom strand)\n", m_reads_OB, ((float) (100*m_reads_OB))/((float) t_reads));
    if(!config.directional) fprintf(fp,"\t%llu\t%6.2f%%\tCTOT (complementary to the original top strand)\n", m_reads_CTOT, ((float) (100*m_reads_CTOT))/((float) t_reads));
    if(!config.directional) fprintf(fp,"\t%llu\t%6.2f%%\tCTOB (complementary to the original bottom strand)\n", m_reads_CTOB, ((float) (100*m_reads_CTOB))/((float) t_reads));
    fprintf(fp,"\n");
    fprintf(fp,"Cytosine Methylation (N.B., statistics from overlapping mates are added together!):\n");
    fprintf(fp,"\tNumber of C's in a CpG context: %llu\n", t_CpG);
    fprintf(fp,"\tPercentage of methylated C's in a CpG context: %6.2f%%\n", ((float) (100*m_CpG))/((float) t_CpG));
    fprintf(fp,"\tNumber of C's in a CHG context: %llu\n", t_CHG);
    fprintf(fp,"\tPercentage of methylated C's in a CHG context: %6.2f%%\n", ((float) (100*m_CHG))/((float) t_CHG));
    fprintf(fp,"\tNumber of C's in a CHH context: %llu\n", t_CHH);
    fprintf(fp,"\tPercentage of methylated C's in a CHH context: %6.2f%%\n", ((float) (100*m_CHH))/((float) t_CHH));

    fclose(fp);
    free(of);
}

tamFile sam_popen(char *cmd) {
    tamFile fp = calloc(1, sizeof(struct __tamFile_t));
    gzFile gzfp;
    int fid, fid2;
    popen_fd = popen(cmd, "r"); //Global

    if(popen_fd == NULL) return 0;
    fid = fileno(popen_fd);
    fid2 = dup(fid); //otherwise, the file descriptor is closed by zlib and pclose() won't work!!
    gzfp = gzdopen(fid2, "r");
    fp->str = (kstring_t*) calloc(1, sizeof(kstring_t));
    fp->fp = gzfp;
    fp->ks = ks_init(fp->fp);
    fp->n_lines = 0;
    fp->is_first = 1;
    return fp;
}

void sam_pclose(tamFile fp) {
    if(fp) {
        ks_destroy(fp->ks);
        gzclose(fp->fp);
        pclose(popen_fd); //global
        free(fp->str->s);
        free(fp->str);
        free(fp);
    }
}
