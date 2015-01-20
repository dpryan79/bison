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
    fastq *read = malloc(sizeof(fastq));
    time_t now;
    char ctime_buffer[26];
    unsigned long long local_m_reads_OT = 0, local_m_reads_OB = 0;
    unsigned long long local_m_reads_CTOT = 0, local_m_reads_CTOB = 0;
    unsigned long long local_total = 0;
    unsigned long long local_concordant = 0, local_discordant = 0, local_singletons = 0;

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
                if(!config.quiet) fprintf(stderr, "Thread %i received a finished signal from node group %i (%i of %i groups are finished)\n", thread_id, j, quit, node_final-node_base); fflush(stderr);
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
                    fprintf(stderr, "%llu reads in thread %i @ %s", local_total, thread_id, ctime_r(&now, ctime_buffer)); fflush(stderr);
                }
            }

            //Get the appropriate read from the linked list
            while(!is_ready(fastq_nodes[j], 0)) sleep(1);
            //Do we need to flush our statistics?
            if(*((char *)(fastq_nodes[j]->next->packed)) == '\2') {
                //Update!
//lock
                pthread_mutex_lock(&metrics_mutex);
                t_concordant += local_concordant;
                t_discordant += local_discordant;
                t_singletons += local_singletons;
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
                local_concordant = 0;
                local_discordant = 0;
                local_singletons = 0;
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

            //Update concordant/discordant/singleton metrics as needed
            if(config.paired) {
                if(best_node & 0xF00 && best_node & 0xFF) local_concordant++;
                else if((best_node&0x11)==0x11 || (best_node&0x22)==0x22 || (best_node&0x44)==0x44 || (best_node&0x88)==0x88) local_discordant++;
                else {
                    if(best_node & 0xF) local_singletons++; //Read#1
                    if(best_node & 0xF0) local_singletons++; //Read#2
                }
            }

            //Store the reads and free up space (N.B., the writer thread will free up the space used by the best read)
            if(!(best_node & 0xF)) { //Unmapped
                move_element(nodes[multiplier*j], to_write_sentinel_node[thread_id]);
            } else if(best_node & 0x1) { //OT#1
                local_m_reads_OT++;
                move_element(nodes[multiplier*j], to_write_sentinel_node[thread_id]);
            } else remove_element(nodes[multiplier*j]);
            if(best_node & 0x2) { //OB#1
                local_m_reads_OB++;
                move_element(nodes[multiplier*j+1], to_write_sentinel_node[thread_id]);
            } else remove_element(nodes[multiplier*j+1]);
            if(!config.directional) {
                if(best_node & 0x4) { //CTOT#1
                    local_m_reads_CTOT++;
                    move_element(nodes[multiplier*j+2], to_write_sentinel_node[thread_id]);
                } else remove_element(nodes[multiplier*j+2]);
                if(best_node & 0x8) { //CTOB#1
                    local_m_reads_CTOB++;
                    move_element(nodes[multiplier*j+3], to_write_sentinel_node[thread_id]);
                } else remove_element(nodes[multiplier*j+3]);
            }
            if(!(best_node & 0xF0) && config.paired) { //Unmapped
                move_element(nodes[multiplier*j], to_write_sentinel_node[thread_id]);
            } else if(best_node & 0x10) { //OT#2
                local_m_reads_OT++;
                move_element(nodes[multiplier*j], to_write_sentinel_node[thread_id]);
            } else if(config.paired) remove_element(nodes[multiplier*j]);
            if(best_node & 0x20) { //OB#2
                local_m_reads_OB++;
                move_element(nodes[multiplier*j+1], to_write_sentinel_node[thread_id]);
            } else if(config.paired) remove_element(nodes[multiplier*j+1]);
            if(!config.directional) {
                if(best_node & 0x40) { //CTOT#2
                    local_m_reads_CTOT++;
                    move_element(nodes[multiplier*j+2], to_write_sentinel_node[thread_id]);
                } else if(config.paired) remove_element(nodes[multiplier*j+2]);
                if(best_node & 0x80) { //CTOB#2
                    local_m_reads_CTOB++;
                    move_element(nodes[multiplier*j+3], to_write_sentinel_node[thread_id]);
                } else if(config.paired) remove_element(nodes[multiplier*j+3]);
            }
            remove_raw_element(fastq_nodes[j]);
        }
    }

    //Update the global metrics
//lock
    pthread_mutex_lock(&metrics_mutex);
    t_concordant += local_concordant;
    t_discordant += local_discordant;
    t_singletons += local_singletons;
    m_reads_OT += local_m_reads_OT;
    m_reads_OB += local_m_reads_OB;
    m_reads_CTOT += local_m_reads_CTOT;
    m_reads_CTOB += local_m_reads_CTOB;
    pthread_mutex_unlock(&metrics_mutex);
//unlock

    //Tell the writer thread that we're finished
    add_finished(to_write_sentinel_node[thread_id]);

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
    if(!config.quiet) fprintf(stderr, "Thread %i finishing!\n", thread_id); fflush(stderr);

    return NULL;
}
