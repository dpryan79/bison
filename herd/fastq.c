#include "../bison.h"
#include <bzlib.h>
#include <zlib.h>
#include <wordexp.h>

//This serve as the buffer for reading from compressed files
struct local_buffer {
    char *buf;
    unsigned long pos;
    int finished; //0 no, 1 yes
    int type; //0: txt, 1: gz, 2:bz2
    union {
        FILE *fptxt;
        gzFile fpgz;
        BZFILE *fpbz2;
    } x;
};

/******************************************************************************
*
*   Take a fastq struct and convert it C->T, the conversion is in place
*
*   fastq *read, input struct
*   int which, which of the reads to convert
*
*******************************************************************************/
void convertCT(fastq *read, int which) {
    char *p;
    if(which == 0) {
        p = read->seq1;
    } else {
        p = read->seq2;
    }
    while(*p != '\n') {
        if(*p == 'C' || *p == 'c') *p = 'T';
        p++;
    }
}

/******************************************************************************
*
*   Take a fastq struct and convert it G->A, the conversion is in place
*
*   fastq *read, input struct
*   int which, which of the reads to convert
*
*******************************************************************************/
void convertGA(fastq *read, int which) {
    char *p;
    if(which == 0) {
        p = read->seq1;
    } else {
        p = read->seq2;
    }
    while(*p != '\n') {
        if(*p == 'G' || *p == 'g') *p = 'A';
        p++;
    }
}

/******************************************************************************
*
*   Read a full line into a buffer, increasing its size as needed and returning
*   its max size.
*
*   FILE *fp, input file stream
*   char *cur_buf, the buffer to expand and insert into
*   int size, current maximum malloc()ed size of cur_buf
*   char *buf, a buffer of length sizeof(char)*MAXREAD to use, this simply saves
*              us from constantly malloc()ing one.
*   int ignore, if 1, read in the line until the end but don't store it
*
*   size is updated on success and set to -1 on error or EOF
*
*******************************************************************************/
char * read_line(struct local_buffer *fp, char *cur_buf, int *size, int ignore) {

    if(fp->type == 0 || fp->type == 2) { //plain text input
        while(1) {
            if(fp->finished == 1) {
                //We hit the end of the file in the last go around
                *size = -1;
                break;
            }
            if(ignore) {
                if(fgets(fp->buf, BT2BUF_SZ, fp->x.fptxt) == NULL) {
                    fp->finished = 1;
                    *size = -1;
                    break;
                }
                while(fp->buf[strlen(fp->buf)-1] != '\n') {
                    if(fgets(fp->buf, BT2BUF_SZ, fp->x.fptxt) == NULL) fp->finished = 1;
                    if(fp->finished == 1) break; //Broken input
                }
                break;
            } else {
                if(fgets(cur_buf, *size, fp->x.fptxt) == NULL) {
                    fp->finished = 1;
                    *size = -1;
                    break;
                }
                while(cur_buf[strlen(cur_buf)-1] != '\n') {
                    cur_buf = realloc(cur_buf, sizeof(char) * (*size + BT2BUF_SZ));
                    *size += BT2BUF_SZ;
                    if(fgets(fp->buf, BT2BUF_SZ, fp->x.fptxt) == NULL) fp->finished = 1;
                    if(fp->finished == 1) break; //Broken input
                    cur_buf = strcat(cur_buf, fp->buf);
                }
                break;
            }
        }
    } else if(fp->type == 1) { //gzipped input
        while(1) {
            if(fp->finished == 1) {
                //We hit the end of the file in the last go around
                *size = -1;
                break;
            }
            if(ignore) {
                if(gzgets(fp->x.fpgz, fp->buf, BT2BUF_SZ) == NULL) {
                    fp->finished = 1;
                    *size = -1;
                    break;
                }
                while(fp->buf[strlen(fp->buf)-1] != '\n') {
                    if(gzgets(fp->x.fpgz, fp->buf, BT2BUF_SZ) == NULL) fp->finished = 1;
                    if(fp->finished == 1) break; //Broken input
                }
                break;
            } else {
                if(gzgets(fp->x.fpgz, cur_buf, *size) == NULL) {
                    fp->finished = 1;
                    *size = -1;
                    break;
                }
                while(cur_buf[strlen(cur_buf)-1] != '\n') {
                    cur_buf = realloc(cur_buf, sizeof(char) * (*size + BT2BUF_SZ));
                    *size += BT2BUF_SZ;
                    if(gzgets(fp->x.fpgz, fp->buf, BT2BUF_SZ) == NULL) fp->finished = 1;
                    if(fp->finished == 1) break; //Broken input
                    cur_buf = strcat(cur_buf, fp->buf);
                }
                break;
            }
        }
    }

    return cur_buf;
}

/******************************************************************************
*
*   Read in an actual fastq read into a fastq struct, resizing as needed
*
*   FILE *fp, input file pointer
*   fastq *read, input struct
*   int which, 0 for read1 and 1 for read2
*
*   returns an int, which is -1 on EOF or error
*
*******************************************************************************/
int read_fastq(struct local_buffer *fp, fastq *read, int which) {
    int *max_name = NULL, *max_seq = NULL, *max_qual = NULL;
    int orig_maxname;
    char *name = NULL, *seq = NULL, *qual = NULL;

    //Point everything to the correct read
    if(which == 0) { //read1
        name = read->name1;
        seq = read->seq1;
        qual = read->qual1;
        max_name = &(read->max_name1);
        max_seq = &(read->max_seq1);
        max_qual = &(read->max_qual1);
    } else {
        name = read->name2;
        seq = read->seq2;
        qual = read->qual2;
        max_name = &(read->max_name2);
        max_seq = &(read->max_seq2);
        max_qual = &(read->max_qual2);
    }

    //name
    orig_maxname = *max_name;
    name = read_line(fp, name, max_name, 0);
    if(*max_name == -1) {
        *max_name = orig_maxname;
        return -1;
    }
    //Seq
    seq = read_line(fp, seq, max_seq, 0);
    //+
    read_line(fp, NULL, 0, 1);
    //Qual
    qual = read_line(fp, qual, max_qual, 0);

    //Reset the pointers if they've moved
    if(which == 0) { //Read1
        read->name1 = name;
        read->seq1 = seq;
        read->qual1 = qual;
    } else {
        read->name2 = name;
        read->seq2 = seq;
        read->qual2 = qual;
    }

    return 0;
}

/******************************************************************************
*
*   Read in the fastq file(s) sending the reads to the appropriate nodes and
*   also storing the unconverted reads in a linked list on the master node.
*
*   This will act as its own thread on the master node.
*
*   void *a is an unsigned long
*
*******************************************************************************/
void * send_store_fastq(void *a) {
    char *line = malloc(MAXREAD*sizeof(char));
    struct local_buffer *f1 = NULL, *f2 = NULL;
    int i=0, nnodes = effective_nodes(), status;
    int nnode_groups = nnodes/((config.directional) ? 2 : 4);
    int j, max_j = 4, multiplier = 4;
    int current_file = 0;
    unsigned long upto = *((unsigned long *) a);
    unsigned long total = 0;
    char *cmd = NULL;
    char *p, *fname1 = NULL, *fname2 = NULL, *save_ptr1=NULL, *save_ptr2=NULL;
    char *finished_signal = NULL;
    fastq *read = malloc(sizeof(fastq));
    MPI_Fastq *packed = NULL;
    int rv1 = 0, rv2 = 0, wordexp_offset=0;
    wordexp_t fnames1_wordexp, fnames2_wordexp;
    void *A = malloc(1);
#ifdef DEBUG
    int taskid = global_debug_taskid;
#endif
    f1 = calloc(1, sizeof(struct local_buffer));
    f1->buf = malloc(BT2BUF_SZ*sizeof(char));
    f1->buf[0] = '\0'; //Just so that we know that we're at the start of a buffer
    if(config.paired) {
        f2 = calloc(1, sizeof(struct local_buffer));
        f2->buf = malloc(BT2BUF_SZ*sizeof(char));
        f2->buf[0] = '\0'; //Just so that we know that we're at the start of a buffer
    }

    //Initialize the read struct
    read->max_name1 = 10;
    read->max_seq1 = 10;
    read->max_qual1 = 10;
    read->max_name2 = 10;
    read->max_seq2 = 10;
    read->max_qual2 = 10;
    read->name1 = malloc(sizeof(char)*10);
    read->seq1 = malloc(sizeof(char)*10);
    read->qual1 = malloc(sizeof(char)*10);
    read->name2 = malloc(sizeof(char)*10);
    read->seq2 = malloc(sizeof(char)*10);
    read->qual2 = malloc(sizeof(char)*10);

    //These will be used later
    if(config.directional) {
        max_j = 2;
        multiplier = 2;
    }

    fname1 = strtok_r(config.FASTQ1,",", &save_ptr1);
    rv1 = wordexp(fname1, &fnames1_wordexp, WRDE_SHOWERR | WRDE_UNDEF);
    fnames1[current_file] = strdup(fnames1_wordexp.we_wordv[wordexp_offset]); //This will need to be free()d
    if(config.paired) {
        fname2 = strtok_r(config.FASTQ2,",", &save_ptr2);
        rv2 = wordexp(fname2, &fnames2_wordexp, WRDE_SHOWERR | WRDE_UNDEF);
        fnames2[current_file] = strdup(fnames2_wordexp.we_wordv[wordexp_offset]);
    }
    while(fname1 != NULL) {
        if(rv1 != 0 || rv2 != 0) {
            printf("An error ocurred when trying to expand the first filename.\n");
            if(rv1 == WRDE_BADCHAR) {
                printf("%s contains an illegal character\n", fname1);
            } else if(rv1 == WRDE_BADVAL) {
                printf("%s contains an undefined shell variable\n", fname1);
            } else if(rv1 == WRDE_NOSPACE) {
                printf("Out of memory when processing %s\n", fname1);
            } else if(rv1 == WRDE_SYNTAX) {
                printf("%s had a syntax error\n", fname1);
            }
            if(config.paired) {
                if(rv2 == WRDE_BADCHAR) {
                    printf("%s contains an illegal character\n", fname2);
                } else if(rv2 == WRDE_BADVAL) {
                    printf("%s contains an undefined shell variable\n", fname2);
                } else if(rv2 == WRDE_NOSPACE) {
                    printf("Out of memory when processing %s\n", fname2);
                } else if(rv2 == WRDE_SYNTAX) {
                    printf("%s had a syntax error\n", fname2);
                }
            }
            goto finish; //Yeah yeah, an evil "goto"
        }
        //Determine how we should read in the file(s)
        p = strrchr(fnames1_wordexp.we_wordv[wordexp_offset], '.');
        if(strcmp(p, ".gz") == 0 || strcmp(p, ".GZ") == 0) {
            f1->type = 1;
            f1->x.fpgz = gzopen(fnames1_wordexp.we_wordv[wordexp_offset], "rb");
        } else if(strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0) {
            f1->type = 2;
            cmd = realloc(cmd, sizeof(char) * (strlen(fnames1_wordexp.we_wordv[wordexp_offset]) + strlen("bzcat  ")));
            sprintf(cmd, "bzcat %s", fnames1_wordexp.we_wordv[wordexp_offset]);
            f1->x.fptxt = popen(cmd, "r");
        } else {
            f1->type = 0;
            f1->x.fptxt = fopen(fnames1_wordexp.we_wordv[wordexp_offset], "r");
        }
        f1->finished = 0;
        f1->pos = 0;
        f1->buf[0] = '\0'; //Just so that we know that we're at the start of a buffer
        if(config.paired) {
            p = strrchr(fnames2_wordexp.we_wordv[wordexp_offset], '.');
            if(strcmp(p, ".gz") == 0 || strcmp(p, ".GZ") == 0) {
                f2->type = 1;
                f2->x.fpgz = gzopen(fnames2_wordexp.we_wordv[wordexp_offset], "rb");
            } else if(strcmp(p, ".bz") == 0 || strcmp(p, ".bz2") == 0) {
                f2->type = 2;
                cmd = realloc(cmd, sizeof(char) * (strlen(fnames2_wordexp.we_wordv[wordexp_offset]) + strlen("bzcat  ")));
                sprintf(cmd, "bzcat %s", fnames2_wordexp.we_wordv[wordexp_offset]);
                f2->x.fptxt = popen(cmd, "r");
            } else {
                f2->type = 0;
                f2->x.fptxt = fopen(fnames2_wordexp.we_wordv[wordexp_offset], "r");
            }
            f2->finished = 0;
            f2->pos = 0;
            f2->buf[0] = '\0'; //Just so that we know that we're at the start of a buffer
        }

        //read everything in
        total = 0;
        while(1) {
            if(upto) {
                if(total >= upto) break;
            }

            if(read_fastq(f1, read, 0) == -1) break;
            if(config.paired) read_fastq(f2, read, 1);
 
            //Pack the struct
            packed = pack_fastq(read);

            //Store this in the linked-list
#ifdef DEBUG
            if(global_debug_taskid == MASTER) {
#endif
            add_element(last_fastq_sentinel_node[i], packed->packed);
#ifdef DEBUG
            }
#endif

            //Send it to the appropriate nodes
            for(j=1; j<=max_j; j++) {
#ifdef DEBUG
                if(global_debug_taskid != MASTER) {
                    if(j+multiplier*i == taskid) {
                        status = MPI_Send((void *) packed->packed, packed->size, MPI_BYTE, 0, 3, MPI_COMM_WORLD);
                        if(status != MPI_SUCCESS) {
                            printf("MPI_Send returned %i\n", status);
                            fflush(stdout);
                        }
                    }
                }
#else
                //Send to j+multiplier*i
                status = MPI_Send((void *) packed->packed, packed->size, MPI_BYTE, j+multiplier*i, 3, MPI_COMM_WORLD);
                if(status != MPI_SUCCESS) {
                    printf("MPI_Send returned %i\n", status);
                    fflush(stdout);
                }
#endif
            }
            i++;
            if(i >= nnode_groups) i=0;

            //Free packed (packed->packed is in the linked list!)
#ifdef DEBUG
            if(global_debug_taskid != MASTER) free(packed->packed);
#endif
            free(packed);
            total++;

#ifndef NOTHROTTLE
            if(config.reads_in_queue > 0) {
                if(total % THROTTLE_CHECK_INTERVAL == 0) {
                    while(total - nwritten[current_file] > config.reads_in_queue) sleep(1);
                }
            }
#endif
        }
        flengths[current_file] = total; //Otherwise, the writer thread will keep waiting
        //Notify the master_processor_threads that they need to update the methylation metrics
        for(j=0; j<nnode_groups; j++) { //This is actually excessive, but we otherwise need to
            finished_signal = malloc(2*sizeof(char)); //We need to malloc() this or it won't be properly free()d after being added to the linked-list.
            sprintf(finished_signal, "\2");
            add_element(last_fastq_sentinel_node[j], (void *) finished_signal);
        }
        if(!config.quiet) printf("finished sending reads from %s (%lu reads)\n", fnames1[current_file], total); fflush(stdout);
        //Close the input files
        if(f1->type == 0) fclose(f1->x.fptxt);
        else if(f1->type == 1) { gzclearerr(f1->x.fpgz); gzclose(f1->x.fpgz); }
        else if(f1->type == 2) pclose(f1->x.fptxt);
        if(config.paired) {
            if(f2->type == 0) fclose(f2->x.fptxt);
            else if(f2->type == 1) { gzclearerr(f2->x.fpgz); gzclose(f2->x.fpgz); }
            else if(f2->type == 2) pclose(f2->x.fptxt);
        }

        current_file++;
        if(++wordexp_offset >= fnames1_wordexp.we_wordc) {
            //Ensure we move to the next file
            wordexp_offset = 0;
            fname1 = strtok_r(NULL,",", &save_ptr1);
            if(fname1 == NULL) break;
            rv1 = wordexp(fname1, &fnames1_wordexp, WRDE_SHOWERR | WRDE_UNDEF | WRDE_REUSE);
            if(config.paired) {
                fname2 = strtok_r(NULL,",", &save_ptr2);
                rv2 = wordexp(fname2, &fnames2_wordexp, WRDE_SHOWERR | WRDE_UNDEF | WRDE_REUSE);
                fnames2[current_file] = strdup(fnames2_wordexp.we_wordv[wordexp_offset]); //This will need to be free()d
            }
        } //Else we've incremented to the next file
        fnames1[current_file] = strdup(fnames1_wordexp.we_wordv[wordexp_offset]); //This will need to be free()d
        if(config.paired) {
            fnames2[current_file] = strdup(fnames2_wordexp.we_wordv[wordexp_offset]); //This will need to be free()d
        }
    }

finish: //We'll only ever "goto" here on an error, otherwise we'll get here normally
    //Send a 1-byte package to signal completion
#ifdef DEBUG
    if(global_debug_taskid != MASTER) {
        status = MPI_Send(A, 1, MPI_BYTE, 0, 3, MPI_COMM_WORLD);
    }
#else
    for(j=1; j<=effective_nodes(); j++) {
        status = MPI_Send(A, 1, MPI_BYTE, j, 3, MPI_COMM_WORLD);
        if(status != MPI_SUCCESS) printf("Couldn't send 'finished' message to worker %i!\n", j);
    }
#endif

    //Add the "finished" element
#ifdef DEBUG
    if(global_debug_taskid == MASTER) {
#endif
        for(i=0; i<nnode_groups; i++) {
            add_finished(last_fastq_sentinel_node[i]);
        }
#ifdef DEBUG
    }
#endif

    //Clean things up
    if(cmd != NULL) free(cmd);
    free(A);
    free(line);
    free(read->name1);
    free(read->seq1);
    free(read->qual1);
    free(read->name2);
    free(read->seq2);
    free(read->qual2);
    free(read);
    wordfree(&fnames1_wordexp);
    if(config.paired) wordfree(&fnames2_wordexp);
    if(!config.quiet) printf("Finished reading in fastq files!\n"); fflush(stdout);
    return NULL;
}
