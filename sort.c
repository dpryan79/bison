#include "bison.h"
#include <math.h>

typedef struct {
    uint32_t l;
    bam1_t **buf;
    int offset;
    char *opref;
} worker_t;

//How many digits are in a number?
//This is only useful if we need to create more than 10000 temp files...
int ndigits(int n) {
    if(n==0) return 1;
    return floor(log10(abs(n)))+1;
}

/*******************************************************************************
*
*  Parse a memory amount specified like 2G, or 1.5M
*
*******************************************************************************/
uint64_t str2Mem(char *s) {
    char *p;
    uint64_t o = strtoull(s, &p, 10);
    if(toupper(*p) == 'K') o = o<<10;
    else if(toupper(*p) == 'M') o = o<<20;
    else if(toupper(*p) == 'G') o = o<<30;
    return o;
}

int alignmentCmp(const void *pa, const void *pb) {
    bam1_t *a = *((bam1_t **) pa), *b = *((bam1_t **) pb);
    int rv = 0;

    if(a->core.tid != b->core.tid) {
        if(a->core.tid<b->core.tid) rv = -1;
        else rv = 1;
    } else {
        if(a->core.pos != b->core.pos) rv = a->core.pos-b->core.pos;
        else if(bam_endpos(a) != bam_endpos(b)) rv = bam_endpos(a) - bam_endpos(b);
    }
    if(!rv) {
        if(((a->core.flag | b->core.flag)&0xC0)==0xC0) { //Read #1 comes first
            if(a->core.flag & BAM_FREAD2) rv = 1;
            else rv = -1;
        }
    }
    return rv;
}

void *worker(void *data) {
    worker_t *w = *((worker_t **) data);
    uint32_t i;
    htsFile *of;
    char *oname = malloc(sizeof(char) * (strlen(w->opref)+strlen("..bam ")+ \
        ((ndigits(w->offset)>4)?ndigits(w->offset):4)));
    assert(oname);

    //Sort
    qsort((void*) w->buf, w->l, sizeof(bam1_t*), alignmentCmp);

    //Write
    sprintf(oname, "%s.%04i.bam", w->opref, w->offset);
    of = sam_open(oname, "wb1");
    sam_hdr_write(of, global_header);
    for(i=0; i<w->l; i++) {
        sam_write1(of, global_header, w->buf[i]);
        bam_destroy1(w->buf[i]);
    }

    //Close
    sam_close(of);
    free(oname);

    return NULL;
}

int sortBuffer(bam1_t **buf, uint32_t nElements, int offset, char *opref) {
    uint32_t i, l = nElements/config.n_compression_threads;
    worker_t **w = malloc(sizeof(worker_t*) * config.n_compression_threads);
    pthread_t *tid = malloc(sizeof(worker_t*) * config.n_compression_threads);
    assert(w);
    assert(tid);

    for(i=0; i<config.n_compression_threads; i++) {
        w[i] = malloc(sizeof(worker_t));
        assert(w[i]);
        w[i]->offset = offset+i;
        w[i]->l = l;
        if(i==config.n_compression_threads-1) w[i]->l = nElements-i*l;
        w[i]->buf = buf+l*i;
        w[i]->opref = opref;
        pthread_create(&tid[i], NULL, worker, w+i);
    }

    for(i=0; i<config.n_compression_threads; i++) {
        pthread_join(tid[i], 0);
        free(w[i]);
    }

    free(w);
    free(tid);

    return offset+config.n_compression_threads;
}

//alignmentBuffer should contain uint32_t l,m as well as bam_t **buf and uint64_t curMem,maxMem
//and int offset!?!?!
alignmentBuffer *pushAlignmentBuffer(alignmentBuffer *buf, bam1_t *b) {
    uint64_t mem = sizeof(bam1_t)+b->m_data + sizeof(bam1_t*);
    if(buf->curMem+mem > buf->maxMem) {
        buf->offset = sortBuffer(buf->buf, buf->l, buf->offset, buf->opref);
        buf->l = 0;
        buf->curMem = sizeof(alignmentBuffer) + buf->m*sizeof(bam1_t*);
    }
    if(buf->l+1>buf->m) {
        buf->m++;
        kroundup32(buf->m);
        buf->buf = realloc(buf->buf, sizeof(bam1_t*)*buf->m);
        assert(buf->buf);
    }
    buf->buf[buf->l++] = b;
    buf->curMem += mem;

    return buf;
}

//Merge stuff
void mergeTemp(alignmentBuffer *buf) {
    int i, found, first, nFinished = 0, nFiles = buf->offset;
    char *opref = buf->opref;
    char *oname = malloc(sizeof(char)*(strlen(opref)+7+ \
        ((ndigits(nFiles)>4)?ndigits(nFiles):4)));
    htsFile **fps;
    bam_hdr_t *tmpHeader;
    bam1_t **bs;
    assert(oname);

    //Take care of any remaining alignments
    if(buf->l) buf->offset = sortBuffer(buf->buf, buf->l, buf->offset, buf->opref);
    nFiles = buf->offset;

    fps = malloc(sizeof(htsFile*)*nFiles);
    bs = malloc(sizeof(bam1_t *)*nFiles);
    assert(fps);
    assert(bs);

    //Open the temp files
    for(i=0; i<nFiles; i++) {
        sprintf(oname, "%s.%04i.bam", opref, i);
        fps[i] = sam_open(oname, "r");
        assert((tmpHeader=sam_hdr_read(fps[i])));
        assert(fps[i]);
        bs[i] = bam_init1();
        assert(sam_read1(fps[i], global_header, bs[i])>1);
    }

    while(nFinished < nFiles) {
        nFinished = found = 0;
        first = -1;
        for(i=0; i<nFiles; i++) {
            if(bs[i] == NULL) {
                nFinished++;
                continue;
            }
            if(!found) {
                found = 1;
                first = i;
            } else {
                if(alignmentCmp((void*)(bs+i), (void*)(bs+first))<0) first = i;
            }
        }
        assert(first >= 0);
        sam_write1(buf->fp, global_header, bs[first]);
        if(sam_read1(fps[first], global_header, bs[first]) <= 1) {
            bam_destroy1(bs[first]);
            bs[first] = NULL;
            nFinished++;
        }
    }

    //Remove the temporary files
    for(i=0; i<nFiles; i++) {
        sam_close(fps[i]);
        sprintf(oname,"%s.%04i.bam", opref, i);
        unlink(oname);
    }

    //reset the buffer
    buf->l = 0;
    buf->curMem = sizeof(alignmentBuffer) + buf->m*sizeof(bam1_t*);
    buf->offset = 0;
        
    free(buf->opref);
    free(oname);
    free(fps);
    free(bs);
}
