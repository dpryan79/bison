#include "bison.h"
#include <sys/time.h>

/******************************************************************************
*
*   The main worker node function.
*
*   int thread_id: the thread_id
*
*******************************************************************************/
void worker_node(int thread_id) {
    int cmd_length = 1, max_qname = 0, status;
    char *cmd, *last_qname = calloc(1, sizeof(char));;
    MPI_Header *packed_header;
    MPI_read *packed_read = calloc(1, sizeof(MPI_read));
    bam_hdr_t *header;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    samFile *fp;
    MPI_Status stat;
#ifdef DEBUG
    int current_p_size = 100;
    htsFile *of;
    bam_hdr_t *debug_header = bam_hdr_init();
    bam1_t *debug_read = bam_init1();
    global_header = bam_hdr_init();
    void *p = calloc(100,1);
    int NODE_ID = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &NODE_ID);
    if(!config.quiet) fprintf(stderr, "NODE_ID: %i\n",NODE_ID); fflush(stderr);
    char *oname;
#else
    int start = 0, i = 0;
#endif
    time_t t0, t1;
    int swapped = 0; //This is only used to indicate paired-end reads in the wrong order!!

    packed_read->size = 0;
    packed_read->packed = NULL;

    //construct the bowtie2 command
    cmd_length += (int) strlen("bowtie2 -q --reorder") + 1;
    cmd_length += (int) strlen(config.bowtie2_options) + 1;
    cmd_length += (int) strlen("--norc -x") + 1;
    cmd_length += (int) strlen(config.genome_dir) + strlen("bisulfite_genome/CT_conversion/BS_CT") + 1;
    cmd_length += (int) 2*(strlen("-1 ") + strlen(config.FASTQ1CT)) + 3;

    cmd = (char *) malloc(sizeof(char) * cmd_length);
    if(thread_id == 1) { //OT Read#1 C->T, Read#2 G->A, Genome C->T only the + strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder %s --norc -x %sbisulfite_genome/CT_conversion/BS_CT -1 %s -2 %s", config.bowtie2_options, config.genome_dir, config.FASTQ1CT, config.FASTQ2GA);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder %s --norc -x %sbisulfite_genome/CT_conversion/BS_CT -U %s", config.bowtie2_options, config.genome_dir, config.FASTQ1CT);
        }
#ifdef DEBUG
        oname = malloc(sizeof(char) *(1+strlen(config.odir)+strlen(config.basename)+strlen("_OT.bam")));
        sprintf(oname, "%s%s_OT.bam", config.odir, config.basename);
        of = sam_open(oname, "wb");
        free(oname);
#endif
    } else if(thread_id == 2) { //OB Read#1 C->T, Read#2 G->A, Genome G->A only the - strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder %s --nofw -x %sbisulfite_genome/GA_conversion/BS_GA -1 %s -2 %s", config.bowtie2_options, config.genome_dir, config.FASTQ1CT, config.FASTQ2GA);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder %s --nofw -x %sbisulfite_genome/GA_conversion/BS_GA -U %s", config.bowtie2_options, config.genome_dir, config.FASTQ1CT);
        }
#ifdef DEBUG
        oname = malloc(sizeof(char) *(1+strlen(config.odir)+strlen(config.basename)+strlen("_OB.bam")));
        sprintf(oname, "%s%s_OB.bam", config.odir, config.basename);
        of = sam_open(oname, "wb");
        free(oname);
#endif
    } else if(thread_id == 3) { //CTOT Read#1 G->A, Read#2 C->T, Genome C->T, only the - strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder %s --nofw -x %sbisulfite_genome/CT_conversion/BS_CT -1 %s -2 %s", config.bowtie2_options, config.genome_dir, config.FASTQ1GA, config.FASTQ2CT);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder %s --nofw -x %sbisulfite_genome/CT_conversion/BS_CT -U %s", config.bowtie2_options, config.genome_dir, config.FASTQ1GA);
        }
#ifdef DEBUG
        oname = malloc(sizeof(char) *(1+strlen(config.odir)+strlen(config.basename)+strlen("_CTOT.bam")));
        sprintf(oname, "%s%s_CTOT.bam", config.odir, config.basename);
        of = sam_open(oname, "wb");
        free(oname);
#endif
    } else if(thread_id == 4) { //CTOB Read#1 G->A, Read#2 C->T, Genome G->A, only the + strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder %s --norc -x %sbisulfite_genome/GA_conversion/BS_GA -1 %s -2 %s", config.bowtie2_options, config.genome_dir, config.FASTQ1GA, config.FASTQ2CT);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder %s --norc -x %sbisulfite_genome/GA_conversion/BS_GA -U %s", config.bowtie2_options, config.genome_dir, config.FASTQ1GA);
        }
#ifdef DEBUG
        oname = malloc(sizeof(char) *(1+strlen(config.odir)+strlen(config.basename)+strlen("_CTOB.bam")));
        sprintf(oname, "%s%s_CTOB.bam", config.odir, config.basename);
        of = sam_open(oname, "wb");
        free(oname);
#endif
    } else {
        fprintf(stderr, "Oh shit, got thread_id %i!\n", thread_id);
        return;
    }

    //Wait for the signal to start
#ifndef DEBUG
    while(start == 0) MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
#endif

    //Start the process
    if(!config.quiet) fprintf(stderr, "Node %i executing: %s\n", thread_id, cmd); fflush(stderr);
    fp = sam_popen(cmd);
    header = sam_hdr_read(fp);
#ifdef DEBUG
    sam_hdr_write(of, header);
#endif

#ifndef DEBUG
    if(thread_id == 1) header = modifyHeader(header, config.argc, config.argv);
    packed_header = pack_header(header);
    if(thread_id == 1) {
        //Send the header
        MPI_Send((void *) &(packed_header->size), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        status = MPI_Send((void *) packed_header->packed, packed_header->size, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
        if(status != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Send returned %i\n", status);
            fflush(stderr);
        }
    }
#else
    packed_header = pack_header(header);
    void *tmp_pointer = malloc(packed_header->size);
    MPI_Request request;
    MPI_Isend((void *) packed_header->packed, packed_header->size, MPI_BYTE, NODE_ID, 2, MPI_COMM_WORLD, &request);
    status = MPI_Recv(tmp_pointer, packed_header->size, MPI_BYTE, NODE_ID, 2, MPI_COMM_WORLD, &stat);
    if(status != MPI_SUCCESS) fprintf(stderr, "We seem to have received an error when sending to ourselves!\n");
    MPI_Wait(&request, &stat);
    unpack_header(debug_header, tmp_pointer);
    global_header = debug_header;
    free(tmp_pointer);
#endif

    t0 = time(NULL);
    if(!config.quiet) fprintf(stderr, "Node %i began sending reads @%s", thread_id, ctime(&t0)); fflush(stderr);
    while(sam_read1(fp, header, read1) == 0) {
#ifdef DEBUG
        sam_write1(of, global_header, read1);
#endif
        if(strcmp(bam_get_qname(read1), last_qname) == 0) { //Multimapper
            if(config.paired) {
                sam_read1(fp, header, read2);
#ifdef DEBUG
                sam_write1(of, global_header, read2);
#endif
            }
            continue;
        } else {
            if(read1->core.l_qname > max_qname) {
                max_qname = read1->core.l_qname + 10;
                last_qname = realloc(last_qname, sizeof(char) * max_qname);
            }
            strcpy(last_qname, bam_get_qname(read1));
        }

        //Are paired-end reads in the wrong order?
        swapped = 0;
        if(config.paired) {
            if(read1->core.flag & BAM_FREAD2) {
                swapped = 1;
                sam_read1(fp, header, read2);
                packed_read = pack_read(read2, packed_read);
#ifndef DEBUG
                MPI_Send((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
#else
                sam_write1(of, global_header, read2);
                if(packed_read->size > current_p_size) p = realloc(p, packed_read->size);
                MPI_Isend((void *) packed_read->packed, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &request);
                status = MPI_Recv(p, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &stat);
                MPI_Wait(&request, &stat);
                debug_read = unpack_read(debug_read, p);
#endif
            }
        }

        //Send the read
        packed_read = pack_read(read1, packed_read);
#ifndef DEBUG
        MPI_Send((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
#else
        if(packed_read->size > current_p_size) p = realloc(p, packed_read->size);
        MPI_Isend(packed_read->packed, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &request);
        status = MPI_Recv(p, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &stat);
        MPI_Wait(&request, &stat);
#endif
        //Deal with paired-end reads
        if(config.paired && !swapped) {
            sam_read1(fp, header, read2);
            packed_read = pack_read(read2, packed_read);
#ifndef DEBUG
            MPI_Send((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
#else
            sam_write1(of, global_header, read2);
            if(packed_read->size > current_p_size) p = realloc(p, packed_read->size);
            MPI_Isend((void *) packed_read->packed, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &request);
            status = MPI_Recv(p, packed_read->size, MPI_BYTE, NODE_ID, 5, MPI_COMM_WORLD, &stat);
            MPI_Wait(&request, &stat);
            debug_read = unpack_read(debug_read, p);
#endif
        }
#ifndef DEBUG
        i++;
#endif
    }
    t1 = time(NULL);
    if(!config.quiet) fprintf(stderr, "Node %i finished sending reads @%s\t(%f sec elapsed)\n", thread_id, ctime(&t1), difftime(t1, t0)); fflush(stderr);

    //Notify the master node
    packed_read->size = 0;
#ifndef DEBUG
    void *A = malloc(1);
    MPI_Send(A, 1, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
    free(A);
#endif

    //Close things up
    bam_hdr_destroy(header);
    bam_destroy1(read1);
    bam_destroy1(read2);
    free(cmd);
    if(packed_read->packed != NULL) free(packed_read->packed);
    free(packed_read);
    if(packed_header->packed != NULL) free(packed_header->packed);
    free(packed_header);
    free(last_qname);
    sam_pclose(fp);
#ifdef DEBUG
    sam_close(of);
    bam_hdr_destroy(debug_header);
    bam_destroy1(debug_read);
    free(p);
#endif
    if(!config.quiet) fprintf(stderr, "Exiting worker node %i\n", thread_id); fflush(stderr);
};
