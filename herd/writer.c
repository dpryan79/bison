#include "../bison.h"
#include <sys/time.h>

int nMergeThreads;
pthread_t *tids;

void *mergeTempWrapper(void *pa) {
    alignmentBuffer *abuf = (alignmentBuffer *) pa;

    //Perform the actual merging
    mergeTemp(abuf);

    sam_close(abuf->fp);
    free(abuf->buf);
    free(abuf);
    return NULL;
}

/******************************************************************************
*
*   Update the CpG/CHG/CHH metrics according to the methylation calls in a read
*
*******************************************************************************/
void herd_update_counts(bam1_t *read) {
    char *XM = bam_aux2Z(bam_aux_get(read, "XM"));
    char base;
    int i;

    for(i=0; i<read->core.l_qseq; i++) {
        base = *(XM+i);
        if(base != '.') {
            if(base == 'Z') {
                t_CpG++;
                m_CpG++;
            } else if(base == 'z') {
                t_CpG++;
            } else if(base == 'X') {
                t_CHG++;
                m_CHG++;
            } else if(base == 'x') {
                t_CHG++;
            } else if(base == 'H') {
                t_CHH++;
                m_CHH++;
            } else if(base == 'h') {
                t_CHH++;
            }
        }
    }
}

void herd_setup(char *fname1, char *fname2) {
    char *cmd = NULL;
    if(config.basename) free(config.basename);
    config.basename = get_basename(fname1);
    config.outname = realloc(config.outname, sizeof(char)*(strlen(config.odir)+ strlen(config.basename)+6));
    assert(config.outname);

    if(config.isCRAM) sprintf(config.outname, "%s%s.cram", config.odir, config.basename);
    else sprintf(config.outname, "%s%s.bam", config.odir, config.basename);
    //Open the output file handles
    if(config.unmapped) {
        create_fastq_names(fname1, fname2);
        cmd = malloc(sizeof(char) * (strlen(config.unmapped1) + 8));
        assert(cmd);
        if(!config.quiet) fprintf(stderr, "Unmapped reads will be written to %s\n", config.unmapped1);
        sprintf(cmd, "gzip > %s", config.unmapped1);
        unmapped1 = popen(cmd, "w");
        if(config.paired) {
            cmd = realloc(cmd, sizeof(char) * (strlen(config.unmapped2) + 8));
            assert(cmd);
            if(!config.quiet) fprintf(stderr, "Unmapped reads will be written to %s\n", config.unmapped2);
            sprintf(cmd, "gzip > %s", config.unmapped2);
            unmapped2 = popen(cmd, "w");
        }
        free(cmd);
    }

    //Open a file for output
    if(!config.isCRAM) OUTPUT_BAM = sam_open(config.outname, "wb");
    else {
        OUTPUT_BAM = sam_open(config.outname, "wc");
        hts_set_fai_filename(OUTPUT_BAM, config.fai);
    }

    if(OUTPUT_BAM == NULL) {
        fprintf(stderr, "Could not open %s for writing!\n", config.outname);
        quit(2,-1);
    }
    if(!config.quiet) fprintf(stderr, "Alignments will be written to %s\n",config.outname);
    //If we have multiple files AND we're sorting then the API only allows a single compression thread
    if(config.n_compression_threads > 1) hts_set_threads(OUTPUT_BAM, config.n_compression_threads);
    global_header = modifyHeader(global_header, config.argc, config.argv);
    sam_hdr_write(OUTPUT_BAM, global_header);
    if(!config.quiet) fprintf(stderr, "Alignment metrics will be printed to %s%s.txt\n",config.odir,config.basename);
    fflush(stderr);
}

/******************************************************************************
*
*   This function will run as its own thread and process the linked lists
*   output from the master processor threads, writing them in order to a BAM
*   file. This will also write all of the other output (aside from metrics).
*   Furthermore, this provides a readout of the current number of reads
*   processed.
*
*   Output is NULL, as is the input (needed by pthreads).
*
*******************************************************************************/
void * bam_writer(void *a) {
    int i, j, *times = malloc(sizeof(int)*config.nmthreads);
    int times_per_thread = effective_nodes();
    int nfinished = 0;
    int nlooped = 0, current_file = 0;
    bam1_t *best_read1 = NULL;
    time_t now;
    char ctime_buffer[26];
    alignmentBuffer *abuf = NULL;
    assert(times);

    nMergeThreads = 0;
    tids = NULL;

    //If we write output in the exact same order as the input, we need to know
    //how many times to write from each master_processor_thread before going to the next
    if(config.directional){
        times_per_thread /= 2;
    } else {
        times_per_thread /= 4;
    }
    for(i=0; i<config.nmthreads; i++) *(times+i) = times_per_thread/config.nmthreads;
    *(times+config.nmthreads-1) += times_per_thread % config.nmthreads;

    //Sleep until we've received the global header
    while(global_header == NULL) sleep(1);
    herd_setup(fnames1[current_file], fnames2[current_file]); //Setup the various names

    if(config.sort) {
        abuf = calloc(1, sizeof(alignmentBuffer));
        assert(abuf);
        abuf->maxMem = config.maxMem;
        abuf->opref = strdup(config.basename);
        abuf->fp = OUTPUT_BAM;
    }

    while(nfinished < config.nmthreads) {
        for(i=0; i<config.nmthreads; i++) {
            for(j=0; j<*(times+i); j++) {
                //Do we need to go to a new file?
                if(flengths[current_file] > 0 && flengths[current_file] == t_reads) {
                    print_metrics();
                    t_reads = 0;
                    t_concordant = 0;
                    t_discordant = 0;
                    t_singletons = 0;
                    m_reads_OT = 0;
                    m_reads_OB = 0;
                    m_reads_CTOT = 0;
                    m_reads_CTOB = 0;
                    t_CpG = 0;
                    m_CpG = 0;
                    t_CHG = 0;
                    m_CHG = 0;
                    t_CHH = 0;
                    m_CHH = 0;
                    //Are we finished?
                    if(is_finished(to_write_node[i])) goto finished;
                    if(unmapped1 != NULL) pclose(unmapped1);
                    if(unmapped2 != NULL) pclose(unmapped2);
                    if(config.sort) {
                        tids = realloc(tids, sizeof(pthread_t)*(++nMergeThreads));
                        assert(tids);
                        pthread_create(&tids[nMergeThreads-1], NULL, mergeTempWrapper, (void*) abuf);
                    }
                    current_file++;
                    herd_setup(fnames1[current_file], fnames2[current_file]);
                    if(config.sort) {
                        //Reinitialize the buffer
                        abuf = calloc(1, sizeof(alignmentBuffer));
                        assert(abuf);
                        abuf->maxMem = config.maxMem;
                        abuf->opref = strdup(config.basename);
                        abuf->fp = OUTPUT_BAM;
                    }
                    i=0;
                    j=0;
                }
                //Just poll every second if we haven't yet written anything or if we've already looped a few times
                if(i == 0 && j == 0) nfinished = 0;
                if(!config.reorder) {
                    if(!is_ready(to_write_node[i], 0)) {
                        if(config.nmthreads == 1) {
                            sleep(1); //This is the same as --reorder
                            break;
                        }
                        if(t_reads == 0) sleep(1);
                        if(++nlooped > 100) {
                            if(nlooped > 1000) nlooped = 1000;
                            sleep(1);
                        }
                        break;
                    }
                } else {
                    while(!is_ready(to_write_node[i], 0)) sleep(1);
                }
                if(is_finished(to_write_node[i])) {
                    nfinished += 1;
                    break;
                }
                best_read1 = to_write_node[i]->next->packed;
                if(!(best_read1->core.flag & BAM_FUNMAP)) {
                    if(config.sort) pushAlignmentBuffer(abuf, bam_dup1(best_read1));
                    else sam_write1(OUTPUT_BAM, global_header, best_read1);
                    herd_update_counts(best_read1);
                } else {
                    if(config.unmapped) {
                        if(best_read1->core.flag & BAM_FREAD1) write_unmapped(unmapped1, best_read1);
                        if(best_read1->core.flag & BAM_FREAD2) write_unmapped(unmapped2, best_read1);
                    }
                }

                remove_element(to_write_node[i]);
                nlooped = 0;
                t_reads++;
                nwritten[current_file]++; //Only keep track of this if we're throttling

                //Give some status
                if((t_reads % 100000) == 0) {
                    now = time(NULL);
                    if(!config.quiet) fprintf(stderr, "%llu reads written @ %s", t_reads, ctime_r(&now, ctime_buffer)); fflush(stderr);
                }
            }
        }
    }

//This isn't elegant, but...
finished:
    if(config.sort) {
        if(!config.quiet) {
            fprintf(stderr, "Waiting for merge threads to complete...");
            fflush(stderr);
        }
        mergeTemp(abuf); //Can't split off a new thread, since the close function will break!
        free(abuf->buf);
        free(abuf);
        for(i=0; i<nMergeThreads; i++) pthread_join(tids[i], 0);
        if(!config.quiet) fprintf(stderr, "\n");
    }

    if(t_reads != 0) print_metrics(); //There seems to be a race condition for the last sample in a list. This gets around that.
    now = time(NULL);
    if(!config.quiet) fprintf(stderr, "Finished writing output @%s", ctime_r(&now, ctime_buffer)); fflush(stderr);

    free(times);
    for(i=0; i<config.nmthreads; i++) destroy_list(to_write_node[i]);

    return NULL;
}
