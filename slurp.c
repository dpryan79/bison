#include "bison.h"

/******************************************************************************
*
*   Add an element to the end of a linked-list
*
*   struct packed_struct *last: last sentinel struct
*   void *packed: a packed read
*
*******************************************************************************/
void add_element(struct packed_struct *last, void *packed) {
    struct packed_struct *new = malloc(sizeof(struct packed_struct));
    struct packed_struct *next_to_last = last->previous;

    //Setup the new element
    new->packed = packed;
    new->next = last;
    new->previous = next_to_last;
    new->state = 0;

    //Update the sentinel struct
    last->previous = new;

    //Update the next_to_last struct
    next_to_last->next = new;
    next_to_last->state = 1;
}

/******************************************************************************
*
*   Destroy a (typically already removed) element from a linked-list
*
*   struct packed_struct *remove: element to destroy
*
*******************************************************************************/
inline void destroy_element(struct packed_struct *remove) {
    bam1_t *pbam1_t = remove->packed;
    if(pbam1_t != NULL) {
        if(pbam1_t->data != NULL) free(pbam1_t->data);
        free(pbam1_t);
    }
    free(remove);
}

/******************************************************************************
*
*   Remove an element from the start of a linked-list
*   is_ready(first, 0) must return 1!
*
*   struct packed_struct *first: first sentinel struct
*
*******************************************************************************/
void remove_element(struct packed_struct *first) {
    struct packed_struct *remove = first->next;
    struct packed_struct *new_next = remove->next;

    first->next = new_next;

    destroy_element(remove);
}


/******************************************************************************
*
*   Is the first or second element ready?
*
*   struct packed_struct *first: first sentinel struct
*   int offset: 0 (first element) or 1 (second element)
*
*    returns 1 for element ready, or 0 otherwise
*
*******************************************************************************/
inline int is_ready(struct packed_struct *first, int offset) {
    if(offset == 0) {
        if(first->next->state == 1) return 1;
    } else {
        if(first->next->next->state == 1) return 1;
    }
    return 0;
}

/******************************************************************************
*
*   Is the linked list finished?
*
*   struct packed_struct *first: first sentinel struct
*
*   returns 1 for finished, 0 otherwise
*
*******************************************************************************/
inline int is_finished(struct packed_struct *first) {
    if(first->next->packed == NULL) return 1;
    return 0;
}

/******************************************************************************
*
*   Add a finished element to a linked list
*
*   struct packed_struct *last: last sentenel struct of targeted list
*
*******************************************************************************/
void add_finished(struct packed_struct *last) {
    struct packed_struct *new = malloc(sizeof(struct packed_struct));
    struct packed_struct *next_to_last = last->previous;

    new->packed = NULL;
    new->next = last;
    new->previous = NULL;
    new->state = 1;

    //Update the sentinel struct
    last->previous = new;

    //Update the next_to_last struct
    next_to_last->next = new;
    next_to_last->state = 1;
    if(config.paired) next_to_last->previous->state = 1;
}

/******************************************************************************
*
*   Initialize a linked list, returning the last sentinel struct
*
*   struct packed_struct *first: first sentinel struct
*
*   returns first sentinel struct
*
*******************************************************************************/
struct packed_struct *initialize_list(struct packed_struct *first) {
    first = malloc(sizeof(struct packed_struct));
    struct packed_struct *last= malloc(sizeof(struct packed_struct));

    first->next = last;
    first->previous = first;
    first->packed = NULL;
    last->next = last;
    last->previous = first;
    last->packed = NULL;

    last->state = 0; //is_ready(last) should always be 0;
    first->state = 0; //is_ready(last) should always be 0;
    return first;
}

/******************************************************************************
*
*   Destroy a linked list of packed_structs
*
*   struct packed_struct *first: linked list to destroy
*
*******************************************************************************/
void destroy_list(struct packed_struct *first) {
    while(first->next->next != first->next) remove_element(first);
    free(first->next);
    free(first);
}

/******************************************************************************
*
*   The MPI receiver thread on the main node
*
*   void *a: NULL input
*
*   returns NULL
*
*******************************************************************************/
void *slurp(void *a) {
    time_t t0, t1;
#ifndef DEBUG
    void *p = NULL;
    int nnodes = (config.directional) ? 2 : 4;
    int nfinished = 0;
    int source = 0;
    int size = 0;
    struct packed_struct *target_node = NULL;
    MPI_Status status;
    int start = 1;
    int ntasks = (config.directional) ? 3: 5;
    int i;
    for(i=1; i<ntasks; i++) {
        if(!config.quiet) printf("Sending start to node %i\n", i); fflush(stdout);
        MPI_Ssend(&start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    //Get the header
    if(MPI_Recv((void *) &size, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
        printf("Received an error when trying to receive header size.\n");
        fflush(stdout);
        quit(3, -2);
    }
    p = malloc((size_t) size);
    if(MPI_Recv(p, size, MPI_BYTE, 1, 2, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
        printf("Received an error when trying to receive header.\n");
        fflush(stdout);
        quit(3, -2);
    }
    global_header = bam_header_init();
    unpack_header(global_header, p);
    free(p);
#else
    char *iname = malloc(sizeof(char) * (1+strlen(config.odir)+strlen(config.basename)+strlen("_CTOT.bam")));
    sprintf(iname, "%s%s_OT.bam", config.odir, config.basename);
    fp1 = bam_open(iname, "r");
    sprintf(iname, "%s%s_OB.bam", config.odir, config.basename);
    fp2 = bam_open(iname, "r");
    if(!config.directional) {
        sprintf(iname, "%s%s_CTOT.bam", config.odir, config.basename);
        fp3 = bam_open(iname, "r");
        sprintf(iname, "%s%s_CTOB.bam", config.odir, config.basename);
        fp4 = bam_open(iname, "r");
    }
    free(iname);
    global_header = bam_header_read(fp1);
    if(!config.quiet) printf("header written\n"); fflush(stdout);
    bam1_t *read = bam_init1();
    MPI_read *packed = calloc(1, sizeof(MPI_read));
    packed->packed = NULL;
    packed->size = 0;
    bam_header_t *tmp;
    tmp = bam_header_read(fp2);
    bam_header_destroy(tmp);
    if(!config.directional) {
        tmp = bam_header_read(fp3);
        bam_header_destroy(tmp);
        tmp = bam_header_read(fp4);
        bam_header_destroy(tmp);
    }
#endif

    //Write a header
    bam_header_write(OUTPUT_BAM, global_header);

    t0 = time(NULL);
    if(!config.quiet) printf("Started slurping @%s", ctime(&t0)); fflush(stdout);
#ifndef DEBUG
    while(nfinished < nnodes) {
        MPI_Probe(MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
        source = status.MPI_SOURCE;
        MPI_Get_count(&status, MPI_BYTE, &size);
        if(source == 1) target_node = node1_last_sentinel;
        else if(source == 2) target_node = node2_last_sentinel;
        else if(source == 3) target_node = node3_last_sentinel;
        else if(source == 4) target_node = node4_last_sentinel;

        if(size > 1) {
            p = malloc((size_t) size);
            MPI_Recv(p, size, MPI_BYTE, source, 5, MPI_COMM_WORLD, &status);
            add_element(target_node, p);
        } else {
            p = malloc((size_t) size);
            MPI_Recv(p, size, MPI_BYTE, source, 5, MPI_COMM_WORLD, &status);
            free(p);
            add_finished(target_node);
            nfinished++;
        }
    }
#else
    //OT
    while(bam_read1(fp1, read) > 1) {
        packed->size = 0;
        packed = pack_read(read, packed);
        add_element(node1_last_sentinel, packed->packed);
        if(config.paired) {
            bam_read1(fp1, read);
            packed->size = 0;
            packed = pack_read(read, packed);
            add_element(node1_last_sentinel, packed->packed);
        }
        //OB
        bam_read1(fp2, read);
        packed->size = 0;
        packed = pack_read(read, packed);
        add_element(node2_last_sentinel, packed->packed);
        if(config.paired) {
            bam_read1(fp2, read);
            packed->size = 0;
            packed = pack_read(read, packed);
            add_element(node2_last_sentinel, packed->packed);
        }
        if(!config.directional) {
            //CTOT
            bam_read1(fp3, read);
            packed->size = 0;
            packed = pack_read(read, packed);
            add_element(node3_last_sentinel, packed->packed);
            if(config.paired) {
                bam_read1(fp3, read);
                packed->size = 0;
                packed = pack_read(read, packed);
                add_element(node3_last_sentinel, packed->packed);
            }
            //CTOB
            bam_read1(fp4, read);
            packed->size = 0;
            packed = pack_read(read, packed);
            add_element(node4_last_sentinel, packed->packed);
            if(config.paired) {
                bam_read1(fp4, read);
                packed->size = 0;
                packed = pack_read(read, packed);
                add_element(node4_last_sentinel, packed->packed);
            }
        }
    }
    bam_destroy1(read);
    free(packed);
    
    add_finished(node1_last_sentinel);
    add_finished(node2_last_sentinel);
    if(!config.directional) {
        add_finished(node3_last_sentinel);
        add_finished(node4_last_sentinel);
    }
#endif
    t1 = time(NULL);
    if(!config.quiet) printf("Finished slurping @%s\t(%f seconds elapsed)\n", ctime(&t1), difftime(t1, t0)); fflush(stdout);
    return NULL;
}
