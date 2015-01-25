#include "../bison.h"

/******************************************************************************
*
*   Take a fastq struct and pack it for shipping
*
*   THE RESULT MUST BE free()d eventually
*
*   fastq *read: The read(s) to store
*   MPI_Fastq *output: the struct into which to pack things
*
*******************************************************************************/
MPI_Fastq * pack_fastq(fastq *read) {
    size_t size = 0;
    size_t length1, length2;
    void *p;
    char *pchar, null_char = '\0';
    MPI_Fastq *output = malloc(sizeof(MPI_Fastq));
    assert(output);

    //Calculate the size needed for read1
    length1 = sizeof(char) * (strlen(read->name1) + strlen(read->seq1) + strlen(read->qual1) + 3);
    size += length1;
    if(config.paired) {
        length2 = sizeof(char) * (strlen(read->name2) + strlen(read->seq2) + strlen(read->qual2) + 3);
        size += length2;
    }
    output->size = size;
    output->packed = malloc(size);
    assert(output->packed);

    //Set everything
    p = output->packed;

    //read1
    memcpy(p, (void *) read->name1, sizeof(char) * (strlen(read->name1)));
    pchar = (char *) p;
    p = (void *) (pchar + strlen(read->name1));
    memcpy(p, (void *) &null_char, sizeof(char));
    pchar = (char *) p;
    p = (void *) (++pchar);
    memcpy(p, (void *) read->seq1, sizeof(char) * (strlen(read->seq1)));
    pchar = (char *) p;
    p = (void *) (pchar + strlen(read->seq1));
    memcpy(p, (void *) &null_char, sizeof(char));
    pchar = (char *) p;
    p = (void *) (++pchar);
    memcpy(p, (void *) read->qual1, sizeof(char) * (strlen(read->qual1)));
    pchar = (char *) p;
    p = (void *) (pchar + strlen(read->qual1));
    memcpy(p, (void *) &null_char, sizeof(char));
    pchar = (char *) p;
    p = (void *) (++pchar);

    //read2
    if(config.paired) {
        memcpy(p, (void *) read->name2, sizeof(char) * (strlen(read->name2)));
        pchar = (char *) p;
        p = (void *) (pchar + strlen(read->name2));
        memcpy(p, (void *) &null_char, sizeof(char));
        pchar = (char *) p;
        p = (void *) (++pchar);
        memcpy(p, (void *) read->seq2, sizeof(char) * (strlen(read->seq2)));
        pchar = (char *) p;
        p = (void *) (pchar + strlen(read->seq2));
        memcpy(p, (void *) &null_char, sizeof(char));
        pchar = (char *) p;
        p = (void *) (++pchar);
        memcpy(p, (void *) read->qual2, sizeof(char) * (strlen(read->qual2)));
        pchar = (char *) p;
        p = (void *) (pchar + strlen(read->qual2));
        memcpy(p, (void *) &null_char, sizeof(char));
        pchar = (char *) p;
        p = (void *) (++pchar);
    }
    return output;
}

/******************************************************************************
*
*   Take unpack a packed fastq struct
*
*   THE RESULT MUST BE free()d
*
*   fastq *read: The fastq struct to unpack into
*   void *packed: The packed structure
*
*******************************************************************************/
fastq * unpack_fastq(fastq *read, void *packed) {
    char *pchar;
    void *p = packed;
    size_t len;

    //Read1
    len = strlen((char *) p) + 1; //name
    if(len > read->max_name1) {
        read->name1 = realloc((void *) read->name1, sizeof(char) * len);
        assert(read->name1);
        read->max_name1 = len;
    }
    strcpy(read->name1, (char *) p);
    pchar = (char *) p;
    p = (void *) (pchar + len);
    len = strlen((char *) p) + 1; //seq
    if(len > read->max_seq1) {
        read->seq1 = realloc((void *) read->seq1, sizeof(char) * len);
        assert(read->seq1);
        read->max_seq1 = len;
    }
    strcpy(read->seq1, (char *) p);
    pchar = (char *) p;
    p = (void *) (pchar + len);
    len = strlen((char *) p) + 1; //qual
    if(len > read->max_qual1) {
        read->qual1 = realloc((void *) read->qual1, sizeof(char) * len);
        assert(read->qual1);
        read->max_qual1 = len;
    }
    strcpy(read->qual1, (char *) p);
    pchar = (char *) p;
    p = (void *) (pchar + len);

    //Read2
    if(config.paired) {
        len = strlen((char *) p) + 1; //name
        if(len > read->max_name2) {
            read->name2 = realloc((void *) read->name2, sizeof(char) * len);
            assert(read->name2);
            read->max_name2 = len;
        }
        strcpy(read->name2, (char *) p);
        pchar = (char *) p;
        p = (void *) (pchar + len);
        len = strlen((char *) p) + 1; //seq
        if(len > read->max_seq2) {
            read->seq2 = realloc((void *) read->seq2, sizeof(char) * len);
            assert(read->seq2);
            read->max_seq2 = len;
        }
        strcpy(read->seq2, (char *) p);
        pchar = (char *) p;
        p = (void *) (pchar + len);
        len = strlen((char *) p) + 1; //qual
        if(len > read->max_qual2) {
            read->qual2 = realloc((void *) read->qual2, sizeof(char) * len);
            assert(read->qual2);
            read->max_qual2 = len;
        }
        strcpy(read->qual2, (char *) p);
        pchar = (char *) p;
        p = (void *) (pchar + len);
    }

    return read;
}
