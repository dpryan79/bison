#include "../bison.h"
#include <math.h>
#include <sys/time.h>

/*******************************************************************************
*
*   The master node function.
*
*   void *a: Actually an int*, the thread_id
*
*******************************************************************************/
void * herd_master_processer_thread(void *a) {
    int thread_id = *((int *) a), best_node, j, quit = 0, multiplier;
    int ngroups = effective_nodes();
    int node_base, node_final;
    int tmp_j = 0;
    char **seq = malloc(sizeof(char *) * 2);
    bam1_t **node1_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node2_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node3_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **node4_read = malloc(sizeof(bam1_t*) * 2);
    bam1_t **best_read = NULL;
    fastq *read = malloc(sizeof(fastq));
    time_t now;
    char ctime_buffer[26];
    unsigned long long local_m_reads_OT = 0, local_m_reads_OB = 0;
    unsigned long long local_m_reads_CTOT = 0, local_m_reads_CTOB = 0;
    unsigned long long local_total = 0;

    //Properly set the number of node groups and other small things
    if(config.directional) {
        ngroups /= 2;
        multiplier = 2;
    } else {
        ngroups /= 4;
        multiplier = 4;
    }
    read->max_name1 = 0;
    read->max_seq1 = 0;
    read->max_qual1 = 0;
    read->max_name2 = 0;
    read->max_seq2 = 0;
    read->max_qual2 = 0;
    read->name1 = NULL;
    read->seq1 = NULL;
    read->qual1 = NULL;
    read->name2 = NULL;
    read->seq2 = NULL;
    read->qual2 = NULL;

    //Get the minimum and maximum node group to work on
    node_base = thread_id*(ngroups/config.nmthreads);
    node_final = node_base+(ngroups/config.nmthreads)-1;
    if(thread_id+1-config.nmthreads == 0) node_final += ngroups % config.nmthreads;
    node_final++;

    //Process read i/o
    while(quit < node_final-node_base) {
        //Currently, we output everything in the same order as the original input
        //We could also invoke one thread per node-group (or multiple groups) and
        //then either output in order or randomly (easier to implement).
        //I'll have to benchmark things to see if this can keep up
        for(j=node_base; j<node_final; j++) {
            while(!is_ready(nodes[multiplier*j], 0));
            if(is_finished(nodes[multiplier*j])) {
                quit += 1;
                if(!config.quiet) printf("Thread %i received a finished signal from node group %i (%i of %i groups are finished)\n", thread_id, j, quit, node_final-node_base); fflush(stdout);
                continue;
            }
            *(node1_read) = update_read(nodes[multiplier*j], 0);
            if(config.paired) {
                while(!is_ready(nodes[multiplier*j], 1));
                *(node1_read+1) = update_read(nodes[multiplier*j], 1);
            }
            while(!is_ready(nodes[multiplier*j+1], 0));
            *(node2_read) = update_read(nodes[multiplier*j+1], 0);
            if(config.paired) {
                while(!is_ready(nodes[multiplier*j+1], 1));
                *(node2_read+1) = update_read(nodes[multiplier*j+1], 1);
            }
            if(!config.directional) {
                while(!is_ready(nodes[multiplier*j+2], 0));
                *node3_read = update_read(nodes[multiplier*j+2], 0);
                if(config.paired) {
                    while(!is_ready(nodes[multiplier*j+2], 1));
                    *(node3_read+1) = update_read(nodes[multiplier*j+2], 1);
                }
                while(!is_ready(nodes[multiplier*j+3], 0));
                *node4_read = update_read(nodes[multiplier*j+3], 0);
                if(config.paired) {
                    while(!is_ready(nodes[multiplier*j+3], 1));
                    *(node4_read+1) = update_read(nodes[multiplier*j+3], 1);
                }
            }

            //Give some output, it's a bit misleading as the count is actually only for this thread and it'll only display for thread 0.
            if(!config.quiet) {
                if(++local_total % 100000 == 0) {
                    now = time(NULL);
                    printf("%llu reads in thread %i @ %s", local_total, thread_id, ctime_r(&now, ctime_buffer)); fflush(stdout);
                }
            }

            //Get the appropriate read from the linked list
            while(!is_ready(fastq_nodes[j], 0)) sleep(1);
            //Do we need to flush our statistics?
            if(*((char *)(fastq_nodes[j]->next->packed)) == '\2') {
                //Update!
//lock
                pthread_mutex_lock(&metrics_mutex);
                m_reads_OT += local_m_reads_OT;
                m_reads_OB += local_m_reads_OB;
                m_reads_CTOT += local_m_reads_CTOT;
                m_reads_CTOB += local_m_reads_CTOB;
                pthread_mutex_unlock(&metrics_mutex);
//unlock
                local_m_reads_OT = 0;
                local_m_reads_OB = 0;
                local_m_reads_CTOT = 0;
                local_m_reads_CTOB = 0;
                local_total = 0;
                tmp_j = j;
                for(j=node_base; j<node_final; j++) {
                    while(!is_ready(fastq_nodes[j], 0)) sleep(1);
                    remove_raw_element(fastq_nodes[j]);
                }
                j = tmp_j;
                while(!is_ready(fastq_nodes[j], 0)) sleep(1);
            }
            read = unpack_fastq(read, (fastq_nodes[j])->next->packed);
            *seq = read->seq1;
            *(*seq + strlen(*seq) - 1) = '\0'; //remove the \n, 
            if(config.paired) {
                *(seq+1) = read->seq2;
                *(*(seq+1) + strlen(*(seq+1)) - 1) = '\0';
            }

            //Process the reads
            if(!config.paired) {
                best_node = process_single(*node1_read, *node2_read, *node3_read, *node4_read, *seq); //Output is stored in read1
            } else {
                best_node = process_paired(node1_read, node2_read, node3_read, node4_read, seq); //Output is stored in read
            }

            if(best_node == 1) {
                best_read = node1_read;
                if(!((*best_read)->core.flag & BAM_FUNMAP)) local_m_reads_OT++;
            } else if(best_node == 2) {
                best_read = node2_read;
                if(!((*best_read)->core.flag & BAM_FUNMAP)) local_m_reads_OB++;
            } else if(best_node == 3) {
                best_read = node3_read;
                if(!((*best_read)->core.flag & BAM_FUNMAP)) local_m_reads_CTOT++;
            } else if(best_node == 4) {
                best_read = node4_read;
                if(!((*best_read)->core.flag & BAM_FUNMAP)) local_m_reads_CTOB++;
            }

            //Store the reads and free up space (N.B., the writer thread will free up the space used by the best read)
            if(best_node != 1) {
                remove_element(nodes[multiplier*j]);
                if(config.paired) remove_element(nodes[multiplier*j]);
            } else {
                move_element(nodes[multiplier*j], to_write_sentinel_node[thread_id]);
            }
            if(best_node != 2) {
                remove_element(nodes[multiplier*j+1]);
                if(config.paired) remove_element(nodes[multiplier*j+1]);
            } else {
                move_element(nodes[multiplier*j+1], to_write_sentinel_node[thread_id]);
            }
            if(!config.directional) {
                if(best_node != 3) {
                    remove_element(nodes[multiplier*j+2]);
                    if(config.paired) remove_element(nodes[multiplier*j+2]);
                } else {
                    move_element(nodes[multiplier*j+2], to_write_sentinel_node[thread_id]);
                }
                if(best_node != 4) {
                    remove_element(nodes[multiplier*j+3]);
                    if(config.paired) remove_element(nodes[multiplier*j+3]);
                } else {
                    move_element(nodes[multiplier*j+3], to_write_sentinel_node[thread_id]);
                }
            }
            remove_raw_element(fastq_nodes[j]);
        }
    }

    //Tell the writer thread that we're finished
    add_finished(to_write_sentinel_node[thread_id]);

    //Update the global metrics
//lock
    pthread_mutex_lock(&metrics_mutex);
    m_reads_OT += local_m_reads_OT;
    m_reads_OB += local_m_reads_OB;
    m_reads_CTOT += local_m_reads_CTOT;
    m_reads_CTOB += local_m_reads_CTOB;
    pthread_mutex_unlock(&metrics_mutex);
//unlock

    //Clean up
    free(seq);
    free(node1_read);
    free(node2_read);
    free(node3_read);
    free(node4_read);
    free(read->name1);
    free(read->seq1);
    free(read->qual1);
    if(config.paired) {
        free(read->name2);
        free(read->seq2);
        free(read->qual2);
    }
    free(read);
    if(!config.quiet) printf("Thread %i finishing!\n", thread_id); fflush(stdout);

    return NULL;
}