#include "../bison.h"
#include <sys/time.h>

struct packed_struct *first_writer, *first_writer_sentinel;
struct packed_struct *second_writer, *second_writer_sentinel;

//write read #1
void * first_writer_func(void *a) {
    int thread_id = ((slurp_fastq_struct *) a)->thread_id;
    char *fastq1 = ((slurp_fastq_struct *) a)->fastq1;
    FILE *f1 = fopen(fastq1, "w");
    fastq *read = malloc(sizeof(fastq));
    int strand;

    //Determine the conversions to make
    if(config.directional) {
        strand = (thread_id-1) % 2;
    } else {
        strand = (thread_id-1) % 4;
    }

    //Initialize the fastq struct
    read->max_name1 = 10;
    read->max_name2 = 10;
    read->max_seq1 = 10;
    read->max_seq2 = 10;
    read->max_qual1 = 10;
    read->max_qual2 = 10;
    read->name1 = malloc(sizeof(char) * 10);
    read->seq1 = malloc(sizeof(char) * 10);
    read->qual1 = malloc(sizeof(char) * 10);
    read->name2 = malloc(sizeof(char) * 10);
    read->seq2 = malloc(sizeof(char) * 10);
    read->qual2 = malloc(sizeof(char) * 10);

    while(1) {
        while(!is_ready(first_writer, 0)); //Sleeping slows things down too much
        if(is_finished(first_writer)) break;

        //Unpack
        read = unpack_fastq(read, first_writer->next->packed);
        //Remove from the linked list
        remove_raw_element(first_writer);
        //Convert
        switch(strand) {
            case 0 :
            case 1 :
                convertCT(read, 0);
                break;
            case 2 :
            case 3 :
                convertGA(read, 0);
                break;
        }
        fprintf(f1, "%s%s+\n%s", read->name1, read->seq1, read->qual1);
    }

    //Free things up
    fclose(f1);
    free(read->name1);
    free(read->seq1);
    free(read->qual1);
    free(read->name2);
    free(read->seq2);
    free(read->qual2);
    free(read);
    destroy_list(first_writer);
    return NULL;
}

//write read #2
void * second_writer_func(void *a) {
    int thread_id = ((slurp_fastq_struct *) a)->thread_id;
    char *fastq2 = ((slurp_fastq_struct *) a)->fastq2;
    FILE *f2 = fopen(fastq2, "w");
    fastq *read = malloc(sizeof(fastq));
    int strand;

    //Determine the conversions to make
    if(config.directional) {
        strand = (thread_id-1) % 2;
    } else {
        strand = (thread_id-1) % 4;
    }

    //Initialize the fastq struct
    read->max_name1 = 10;
    read->max_name2 = 10;
    read->max_seq1 = 10;
    read->max_seq2 = 10;
    read->max_qual1 = 10;
    read->max_qual2 = 10;
    read->name1 = malloc(sizeof(char) * 10);
    read->seq1 = malloc(sizeof(char) * 10);
    read->qual1 = malloc(sizeof(char) * 10);
    read->name2 = malloc(sizeof(char) * 10);
    read->seq2 = malloc(sizeof(char) * 10);
    read->qual2 = malloc(sizeof(char) * 10);

    while(1) {
        while(!is_ready(second_writer, 0)); //Sleeping slows things down too much
        if(is_finished(second_writer)) break;

        //Unpack
        read = unpack_fastq(read, second_writer->next->packed);
        //Remove from the linked list
        remove_raw_element(second_writer);
        //Convert
        switch(strand) {
            case 0 :
            case 1 :
                convertGA(read, 1);
                break;
            case 2 :
            case 3 :
                convertCT(read, 1);
                break;
        }
        fprintf(f2, "%s%s+\n%s", read->name2, read->seq2, read->qual2);
    }

    //Free things up
    fclose(f2);
    free(read->name1);
    free(read->seq1);
    free(read->qual1);
    free(read->name2);
    free(read->seq2);
    free(read->qual2);
    free(read);
    destroy_list(second_writer);
    return NULL;
}

/******************************************************************************
*
*   This receives the reads, converts them, and writes them to the FIFO(s)
*   
*   void *a: a pointer to a struct with the following components:
*
*   int thread_id: the thread_id
*   char *fastq1: FIFO from which bowtie2 can get read1
*   char *fastq2: FIFO from which bowtie2 can get read2 (if it exists)
*
*******************************************************************************/
void * slurp_fastq(void *a) {
    pthread_t threads[2];
    void *p = NULL, *p2 = NULL, *p3 = NULL;
    int size = 0, current_p_size = 0;
    MPI_Status status;
    fastq *read = malloc(sizeof(fastq));

    first_writer = malloc(sizeof(struct packed_struct));
    first_writer_sentinel = malloc(sizeof(struct packed_struct));
    first_writer = initialize_list(first_writer);
    first_writer_sentinel = first_writer->next;
    pthread_create(&(threads[0]), NULL, &first_writer_func, a);
    if(config.paired) {
        //If we have pairs, then writing simultaneuosly to two fifos (that will be read sequentially by bowtie2) won't work, since bowtie2 will read from a single fifo multiple times!!!
        second_writer = malloc(sizeof(struct packed_struct));
        second_writer_sentinel = malloc(sizeof(struct packed_struct));
        second_writer = initialize_list(second_writer);
        second_writer_sentinel = second_writer->next;
        pthread_create(&(threads[1]), NULL, &second_writer_func, a);
    }

    //Initialize the fastq struct
    read->max_name1 = 10;
    read->max_name2 = 10;
    read->max_seq1 = 10;
    read->max_seq2 = 10;
    read->max_qual1 = 10;
    read->max_qual2 = 10;
    read->name1 = malloc(sizeof(char) * 10);
    read->seq1 = malloc(sizeof(char) * 10);
    read->qual1 = malloc(sizeof(char) * 10);
    read->name2 = malloc(sizeof(char) * 10);
    read->seq2 = malloc(sizeof(char) * 10);
    read->qual2 = malloc(sizeof(char) * 10);

    //Receive and process the raw reads
    while(1) {
        MPI_Probe(0, 3, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_BYTE, &size);
        if(size > current_p_size) {
            p = realloc(p, (size_t) size);
        }
        MPI_Recv(p, size, MPI_BYTE, 0, 3, MPI_COMM_WORLD, &status);
        //Are we finished receiving?
        if(size <= 1) break;

        //Copy if needed
        if(config.paired) {
            p2 = malloc(size);
            memcpy(p2,p,size);
            add_element(second_writer_sentinel, p2);
        }
        p3 = malloc(size);
        memcpy(p3,p,size);
        add_element(first_writer_sentinel, p3);
    }
    add_finished(first_writer_sentinel);
    if(config.paired) add_finished(second_writer_sentinel);

    //Wait for the other thread
    pthread_join(threads[0], NULL);
    if(config.paired) {
        pthread_join(threads[1], NULL);
    }

    //Free things up
    free(p);
    free(read->name1);
    free(read->seq1);
    free(read->qual1);
    free(read->name2);
    free(read->seq2);
    free(read->qual2);
    free(read);

    return NULL;
}

/******************************************************************************
*
*   The main worker node function.
*
*   int thread_id: the thread_id
*   char *fastq1: FIFO from which bowtie2 can get read1
*   char *fastq2: FIFO from which bowtie2 can get read2 (if it exists)
*
*******************************************************************************/
void herd_worker_node(int thread_id, char *fastq1, char *fastq2) {
    int cmd_length = 1, max_qname = 0, status, strand;
    char *cmd, *last_qname = calloc(1, sizeof(char));
    MPI_Header *packed_header;
    MPI_read *packed_read = calloc(1, sizeof(MPI_read));
    bam_header_t *header;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    tamFile fp;
#ifdef DEBUG
    MPI_Status stat;
    int current_p_size = 100;
    bamFile of;
    bam_header_t *debug_header = bam_header_init();
    bam1_t *debug_read = bam_init1();
    global_header = bam_header_init();
    void *p = calloc(100,1);
    char *oname = NULL;
#else
    int i = 0;
#endif
    time_t t0, t1;

    //Which strand should we be aligning to?
    if(config.directional) {
        strand = (thread_id-1) % 2;
    } else {
        strand = (thread_id-1) % 4;
    }

    packed_read->size = 0;
    packed_read->packed = NULL;

    //construct the bowtie2 command
    cmd_length += (int) strlen("bowtie2 -q --reorder --no-mixed --no-discordant") + 1;
    cmd_length += (int) strlen(config.bowtie2_options) + 1;
    cmd_length += (int) strlen("--norc -x") + 1;
    cmd_length += (int) strlen(config.genome_dir) + strlen("bisulfite_genome/CT_conversion/BS_CT") + 1;
    cmd_length += (int) 2*(strlen("-1 ") + strlen(fastq1)) + 3;
    if(config.paired) cmd_length += (int) strlen(fastq2); //This is likely unneeded.

#ifdef DEBUG
    oname = malloc(sizeof(char) *(1+strlen(config.odir)+strlen(config.basename)+strlen("_X.bam")));
    sprintf(oname, "%s%s_%i.bam", config.odir, config.basename, thread_id);
    if(!config.quiet) printf("Writing output to %s\n", oname);
    of = bam_open(oname, "w");
    free(oname);
#endif

    cmd = (char *) malloc(sizeof(char) * cmd_length);
    if(strand == 0) { //OT Read#1 C->T, Read#2 G->A, Genome C->T only the + strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --norc -x %sbisulfite_genome/CT_conversion/BS_CT -1 %s -2 %s", config.bowtie2_options, config.genome_dir, fastq1, fastq2);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --norc -x %sbisulfite_genome/CT_conversion/BS_CT -U %s", config.bowtie2_options, config.genome_dir, fastq1);
        }
    } else if(strand == 1) { //OB Read#1 C->T, Read#2 G->A, Genome G->A only the - strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --nofw -x %sbisulfite_genome/GA_conversion/BS_GA -1 %s -2 %s", config.bowtie2_options, config.genome_dir, fastq1, fastq2);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --nofw -x %sbisulfite_genome/GA_conversion/BS_GA -U %s", config.bowtie2_options, config.genome_dir, fastq1);
        }
    } else if(strand == 2) { //CTOT Read#1 G->A, Read#2 C->T, Genome C->T, only the - strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --nofw -x %sbisulfite_genome/CT_conversion/BS_CT -1 %s -2 %s", config.bowtie2_options, config.genome_dir, fastq1, fastq2);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --nofw -x %sbisulfite_genome/CT_conversion/BS_CT -U %s", config.bowtie2_options, config.genome_dir, fastq1);
        }
    } else if(strand == 3) { //CTOB Read#1 G->A, Read#2 C->T, Genome G->A, only the + strand
        if(config.paired) {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --norc -x %sbisulfite_genome/GA_conversion/BS_GA -1 %s -2 %s", config.bowtie2_options, config.genome_dir, fastq1, fastq2);
        } else {
            sprintf(cmd, "bowtie2 -q --reorder --no-mixed --no-discordant %s --norc -x %sbisulfite_genome/GA_conversion/BS_GA -U %s", config.bowtie2_options, config.genome_dir, fastq1);
        }
    } else {
        printf("Oh shit, got strand %i!\n", strand);
        return;
    }

    //Start the process
    if(!config.quiet) printf("Node %i executing: %s\n", thread_id, cmd); fflush(stdout);
    fp = sam_popen(cmd);
    header = sam_header_read(fp);
#ifdef DEBUG
    bam_header_write(of, header);
#endif

#ifndef DEBUG
    packed_header = pack_header(header);
    if(thread_id == 1) {
        //Send the header
        MPI_Send((void *) &(packed_header->size), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        status = MPI_Send((void *) packed_header->packed, packed_header->size, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
        if(status != MPI_SUCCESS) {
            printf("MPI_Send returned %i\n", status);
            fflush(stdout);
        }
    }
#else
    packed_header = pack_header(header);
    void *tmp_pointer = malloc(packed_header->size);
    MPI_Request request;
    MPI_Isend((void *) packed_header->packed, packed_header->size, MPI_BYTE, 0, 2, MPI_COMM_WORLD, &request);
    status = MPI_Recv(tmp_pointer, packed_header->size, MPI_BYTE, 0, 2, MPI_COMM_WORLD, &stat);
    if(status != MPI_SUCCESS) printf("We seem to have not been able to send the message to ourselves!\n");
    MPI_Wait(&request, &stat);
    unpack_header(debug_header, tmp_pointer);
    global_header = debug_header;
    free(tmp_pointer);
#endif

    t0 = time(NULL);
    if(!config.quiet) printf("Node %i began sending reads @%s", thread_id, ctime(&t0)); fflush(stdout);
    while(sam_read1(fp, header, read1) > 1) {
#ifdef DEBUG
        bam_write1(of, read1);
#endif
        if(strcmp(bam1_qname(read1), last_qname) == 0) { //Multimapper
            if(config.paired) {
                sam_read1(fp, header, read2);
#ifdef DEBUG
                bam_write1(of, read2);
#endif
            }
            continue;
        } else {
            if(read1->core.l_qname > max_qname) {
                max_qname = read1->core.l_qname + 10;
                last_qname = realloc(last_qname, sizeof(char) * max_qname);
            }
            strcpy(last_qname, bam1_qname(read1));
        }

        //Send the read
        packed_read = pack_read(read1, packed_read);
#ifndef DEBUG
        MPI_Send((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
#else
        if(packed_read->size > current_p_size) p = realloc(p, packed_read->size);
        MPI_Isend(packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD, &request);
        status = MPI_Recv(p, packed_header->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD, &stat);
        MPI_Wait(&request, &stat);
#endif
        //Deal with paired-end reads
        if(config.paired) {
            sam_read1(fp, header, read2);
            packed_read = pack_read(read2, packed_read);
#ifndef DEBUG
            MPI_Send((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
#else
            bam_write1(of, read2);
            if(packed_read->size > current_p_size) p = realloc(p, packed_read->size);
            MPI_Isend((void *) packed_read->packed, packed_read->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD, &request);
            status = MPI_Recv(p, packed_header->size, MPI_BYTE, 0, 5, MPI_COMM_WORLD, &stat);
            MPI_Wait(&request, &stat);
            debug_read = unpack_read(debug_read, p);
#endif
        }
#ifndef DEBUG
        i++;
#endif
    }
    t1 = time(NULL);
    if(!config.quiet) printf("Node %i finished sending reads @%s\t(%f sec elapsed)\n", thread_id, ctime(&t1), difftime(t1, t0)); fflush(stdout);

    //Notify the master node
    packed_read->size = 0;
#ifndef DEBUG
    void *A = malloc(1);
    MPI_Send(A, 1, MPI_BYTE, 0, 5, MPI_COMM_WORLD);
    free(A);
#endif

    //Close things up
    bam_header_destroy(header);
    bam_destroy1(read1);
    bam_destroy1(read2);
    free(cmd);
    if(packed_read->packed != NULL) free(packed_read->packed);
    free(packed_read);
    if(packed_header->packed != NULL) free(packed_header->packed);
    free(packed_header);
    free(last_qname);
    sam_pclose(fp);
    //Remove the FIFO(s)
    unlink(fastq1);
    if(config.paired) unlink(fastq2);
#ifdef DEBUG
    bam_close(of);
    bam_header_destroy(debug_header);
    bam_destroy1(debug_read);
    free(p);
#endif
    if(!config.quiet) printf("Exiting worker node %i\n", thread_id); fflush(stdout);
};
