#include "bison.h"

/******************************************************************************
*
*   Take a BAM header and pack it into a single contiguous memory block. Store
*   the resulting block and its size in an MPI_Header structure.
*
*   THE RESULT MUST BE free()d
*
*   bam_hdr_t *header: The header to store
*
*******************************************************************************/
MPI_Header * pack_header(bam_hdr_t *header) {
    size_t size = sizeof(int32_t); //n_targets
    int32_t *pint32_t;
    uint32_t *puint32_t;
    char *pchar;
    int *pint;
    int i;
    void *p;
    MPI_Header *output = malloc(sizeof(MPI_Header));
    assert(output);

    //target_name
    for(i=0; i<header->n_targets; i++) {
        size += (sizeof(char) * (1+strlen(header->target_name[i])));
    }

    //target_len
    size += sizeof(uint32_t) * header->n_targets;

    //l_text
    size += sizeof(int);

    //text
    size += sizeof(char) * (1 + header->l_text);

    //Start copying, layout is n_targets,target_name[s],target_len[s],l_text,text
    output->size = (int) size;
    output->packed = malloc(size);
    assert(output->packed);
    p = output->packed;

    //n_targets
    memcpy(p, (void *) &(header->n_targets), sizeof(int32_t));
    pint32_t = (int32_t *) p;
    p = (void *) (++pint32_t);

    //target_name
    for(i=0; i<header->n_targets; i++) {
        memcpy(p, (void *) header->target_name[i], sizeof(char) * (1 + strlen(header->target_name[i])));
        pchar = (char *) p;
        p = (void *) (pchar+1+strlen(header->target_name[i]));
    }
    //target_len
    memcpy(p, (void *) header->target_len, sizeof(uint32_t)*(header->n_targets));
    puint32_t = (uint32_t *) p;
    p = (void *) (puint32_t + header->n_targets);

    //l_text
    memcpy(p, (void *) &(header->l_text), sizeof(int));
    pint = (int *) p;
    p = (void *) ++pint;

    //text
    memcpy(p, (void *) (header->text), sizeof(char) * (header->l_text));

    return output;
}

/******************************************************************************
*
*   Unpack a header packed into an initialized bam_hdr_t
*
*   bam_hdr_t *header: The header to unpack into
*   void *packed: The packed header
*
*******************************************************************************/
void unpack_header(bam_hdr_t *header, void *packed) {
    void *p = packed;
    int i;
    int *pint;
    int32_t *pint32_t;
    uint32_t *puint32_t;
    char *pchar;
    size_t strlength;

    //n_targets
    header->n_targets = *((int32_t *) packed);
    pint32_t = (int32_t *) p;
    p = (void *) (++pint32_t);

    //**target_name
    header->target_name = (char **) malloc(sizeof(char *) * (header->n_targets));
    assert(header->target_name);
    for(i=0; i<header->n_targets; i++) {
        strlength = strlen((char *) p)+1;
        header->target_name[i] = malloc(sizeof(char) * strlength);
        assert(header->target_name[i]);
        memcpy((void *) (header->target_name[i]), p, sizeof(char)*strlength);
        pchar = (char *) p;
        p = (void *) (pchar+strlength);
    }

    //target_len
    header->target_len = malloc(sizeof(uint32_t) * (header->n_targets));
    assert(header->target_len);
    for(i=0; i<header->n_targets; i++) {
        header->target_len[i] = *((uint32_t *) p);
        puint32_t = (uint32_t *) p;
        p = (void *) ++puint32_t;
    }

    //l_text
    header->l_text = *((int *) p);
    pint = (int *) p;
    p = (void *) ++pint;

    //text
    header->text = (char *) malloc(sizeof(char) * (header->l_text+1));
    assert(header->text);
    memcpy((void *) (header->text), p, sizeof(char) * (header->l_text + 1));
}

/******************************************************************************
*
*   Take a BAM read and pack it into a single contiguous memory block. Store
*   the resulting block and its size in an MPI_Read structure.
*
*   THE RESULT MUST BE free()d
*
*   bam1_t *read: The read to store
*
*******************************************************************************/
MPI_read * pack_read(bam1_t *read, MPI_read *output) {
    bam1_t *pbam1_t;
    int needed_size, m_data = read->m_data;

    needed_size = (int) (sizeof(bam1_t) + m_data);
    if(output->size == 0) {
        output->packed = malloc((size_t) needed_size);
        assert(output->packed);
        output->size = needed_size;
    } else if(needed_size > output->size) {
        output->packed = realloc(output->packed, (size_t) needed_size);
        assert(output->packed);
        output->size = needed_size;
    }
    memcpy((void *) output->packed, (void *) read, sizeof(bam1_t));
    pbam1_t = output->packed;
    pbam1_t++;
    memcpy((void *) pbam1_t, (void *) read->data, m_data);
    return output;
}

/******************************************************************************
*
*   Unpack a packed read into an initialized bam1_t read.
*
*   bam1_t *read: The read to unpack into
*   void *packed: The packed read
*
*******************************************************************************/
bam1_t *unpack_read(bam1_t *read, void *packed) {
    bam1_t *pbam1_t = packed;
    uint8_t *pdata = (uint8_t *) (pbam1_t+1);
    uint8_t *newdata;

    pbam1_t->data = pdata;
    if(read != NULL) bam_destroy1(read);
    read = bam_init1();
    read->core = pbam1_t->core;
    read->l_data= pbam1_t->l_data;
    read->m_data = pbam1_t->m_data;
    newdata = (uint8_t *) malloc(read->m_data);
    assert(newdata);
    memcpy((void *) newdata, (void *) pdata, read->m_data);
    read->data = newdata;

    return read;
}
