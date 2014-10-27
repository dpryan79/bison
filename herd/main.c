#include "../bison.h"
#include <wordexp.h>

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] -g genome_dir/ {-1 fastq_A1.gz,fastq_B1.gz -2 fastq_A2.gz,fastq_B2.gz | -U fastq.gz}\n", prog);
    printf("\n \
    N.B., Any option not listed below will be passed directly to bowtie2, so you\n\
    can specify, e.g., --very-fast if you want. If you specify --local,\n\
    --score-min is changed back to the bowtie2 default of 'G,20,6', unless you\n\
    specify otherwise.\n\
\n\
    Note also that both -1/-2 and -U can accept a comma-separated list of input\n\
    files. Unlike other aligners, the alignments from each of these files will\n\
    be output to different files. This is meant to speed alignments of multiple\n\
    samples, since the bowtie2 index and the genome sequence only need to be\n\
    loaded a single time. Inputting more than one file (or pair, when using -1\n\
    -2) implies --reorder.\n\
\n\
-g          Directory containing the genome fasta files and the\n\
            Bisulfite_Sequences directory.\n\
\n\
-1          Fastq file containing read #1 (normally named something like \n\
            foo_1.fastq.gz). Reads needn't be gzipped, but that'll be more\n\
            convenient. You may also input a comma-separated list of files to be\n\
            aligned (but see note above). Doing this implies --reorder.\n\
\n\
-2          As with -1, but with read #2.\n\
\n\
-U          For convenience, this denotes a fastq file from single-ended reads.\n\
            Alternatively, -1 can be used without using -2. As with -1, you may\n\
            also specify more than one file, in which case alignments from each\n\
            will be printed to different files.\n\
\n\
-p          How many threads bowtie2 should use on each node. Default is 11.\n\
\n\
-mp         How many processing threads should run on the master node. Default\n\
            is 1. Increasing this will be required to prevent the MPI buffer\n\
            from becoming depleted and the master node then crashing. However,\n\
            too many of these will cause resource underutilization. Keep in\n\
            mind also that there are an additional 2 threads already running to\n\
            do other things.\n\
\n\
-o          Output directory. By default, everything will be written to the\n\
            directory holding the fastq files (or the file containing read #1,\n\
            as appropriate). If you would prefer for the output BAM file and\n\
            metrics txt file to be placed elsewhere, specify that here.\n\
\n\
            N.B., the directory must exist! \n\
\n\
-tmp        Temporary directory where named pipes will be created on the worker\n\
            nodes. This just need to be a directory that is bison_herd can read\n\
            and write to. The default is \"/tmp\".\n\
\n\
--directional Denotes that the library was created in a directional, rather\n\
            than non-directional manner. This will result in 3, rather than 5\n\
            nodes being used as only alignments to 2 (rather than 4) strands are\n\
            possible.\n\
\n\
-upto       The maximum number of reads to process. This is mostly useful for\n\
            debuging and more quickly determining if a library is directional or\n\
            not. 0 is the default, meaning all reads are used. N.B., the\n\
            maximum value for this parameter is whatever an unsigned long is on\n\
            your system.\n\
\n\
--reorder   Reorder output to match the same order as the input. This will make\n\
            things slower, but enable easier comparisons. This is passed to\n\
            bowtie2 regardless of whether you specify it or not. If you use\n\
            multiple input files then this option will always be used, even if\n\
            unspecified.\n\
\n\
-@          Number of BAM compression threads to use. This is equivalent to -@\n\
            in samtools. The default is 1, but this may need to be increased as\n\
            you increase the number of alotted nodes.\n\
\n\
--no-discordant Suppress discordant alignments. This is actually a bowtie2\n\
            option.\n\
\n\
--no-mixed  Suppress singleton alignments. This is actually a bowtie2 option.\n\
\n\
--unmapped  Save unaligned reads to a file or files (as appropriate). This files\n\
            will be placed in the same directory as the source fastq files,\n\
            regardless of whether \"-o\" is used.\n\
\n\
--genome-size Many of the bison tools need to read the genome into memory. By\n\
            default, they allocate 3000000000 bases worth of memory for this and\n\
            increase that as needed. However, this can sometimes be far more\n\
            than is needed (meaning wasted memory) or far too little (in which\n\
            case the process can become quite slow). If you input the\n\
            approximate size of your genome here (in bases), then you can\n\
            maximize performance and minimize wasted space. It's convenient to\n\
            round up a little.\n\
\n");
#ifndef NOTHROTTLE
    printf("\
-queue_size The maximum difference between the number of reads that have been\n\
            read and the number that have been written. The default is 1000000\n\
            and a value of 0 (or just not compiling with -DTHROTTLE) will\n\
            disable this. Since bison_herd can have a quiet large number of\n\
            worker nodes performing alignments, it can happen that they\n\
            overwhelm the master node that must then process their results. This\n\
            option can help to prevent that (though increasing -mp is a better\n\
            solution) by pausing the sending of reads out for alignment.\n\
\n");
#endif
    printf("\
--quiet     Don't print anything but errors to the console (this is also passed\n\
            to bowtie2).\n\
\n\
-h          Print this help message.\n\
\n\
-v          Print version information.\n\
\n");
#ifdef DEBUG
    printf("\
-taskid     Which node number to act as. The default is 0, the master node.\n\
            Other possibilities are 1-4, which are the worker nodes that\n\
            process OT, OB, CTOT, and CTOB alignments, respectively.\n\
\n\
            Note that if you plan to run with taskid=0 (i.e., as the master\n\
            node), files named OT.bam, OB.bam, etc. should exist in your\n\
            working directory. These will be created automatically if you run\n\
            each pseudo-worker node first, which is recommended.\n\n");
#endif
}

int main(int argc, char *argv[]) {
    int i, taskid=0, provided;
    pthread_t *threads;
    int bowtie2_options_max = MAXREAD;
    char *p = NULL, *tmp = NULL;
    wordexp_t p_wordexp;
    unsigned long upto = 0;
    int ngroups;
    int multi_file=0;
#ifndef DEBUG
    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif

    //Deal with MPI initialization, this seems like an odd way to do things.
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    if(provided != MPI_THREAD_MULTIPLE) {
        printf("You're MPI implementation doesn't support MPI_THREAD_MULTIPLE, which is required for bison_herd to work.\n");
        return -1;
    }
#ifndef DEBUG
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(processor_name, &name_len);
#endif

    config.odir = NULL;
    config.paired = 0; //Default is single-ended
    config.directional = 0; //Default is non-directional
    config.nthreads = 11; //Default is 11 threads/node
    config.bowtie2_options = calloc(MAXREAD, sizeof(char));
    config.unmapped = 0; //By default, unmapped reads are NOT written to a fastq file
    config.scoremin_type = 'L'; //--score-min 'L,-0.6,-0.6'
    config.scoremin_intercept = -0.6;
    config.scoremin_coef = -0.6;
    config.mode = 0; //--end-to-end
    config.tmpdir = NULL; //-tmpdir
    config.nmthreads = 1; //-mp
    config.reorder = 0; //--reorder
    config.outname = NULL; //Otherwise, we'll have problems when we realloc!
    config.basename = NULL; //To handle multiple inputs
    config.n_compression_threads = 0;
    config.unmapped1 = NULL;
    config.unmapped2 = NULL;
    config.genome_dir = NULL;
    global_header = NULL;
    unmapped1 = NULL;
    unmapped2 = NULL;
#ifndef NOTHROTTLE
    config.reads_in_queue = 1000000;
    nwritten = 0;
#endif
    chromosomes.max_genome = 3000000000;
    chromosomes.nchromosomes = 0; //We need to initialize the struct

    //These are only used during cleanup and will otherwise cause an error
    config.FASTQ1CT = NULL;
    config.FASTQ1GA = NULL;
    config.FASTQ2CT = NULL;
    config.FASTQ2GA = NULL;

    //Initialize the global counts
    t_reads = 0;
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

    if(argc == 1) {
        usage(argv[0]);
        quit(0, 0);
    }

    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            quit(0, 0);
        } else if(strcmp(argv[i], "-v") == 0) {
            version();
            quit(0,0);
        } else if(strcmp(argv[i], "-1") == 0) {
            i++;
            config.FASTQ1 = argv[i];
        } else if(strcmp(argv[i], "-2") == 0) {
            i++;
            config.FASTQ2 = argv[i];
            config.paired = 1;
        } else if(strcmp(argv[i], "-U") == 0) {
            i++;
            config.FASTQ1 = argv[i];
        } else if(strcmp(argv[i], "-g") == 0) {
            i++;
            config.genome_dir = strdup(argv[i]);
            if(*(config.genome_dir+strlen(config.genome_dir)-1) != '/') {
                config.genome_dir = realloc(config.genome_dir, sizeof(char) * (strlen(config.genome_dir)+2));
                sprintf(config.genome_dir, "%s/", config.genome_dir);
            }
        } else if(strcmp(argv[i], "-p") == 0) {
            i++;
            config.nthreads = atoi(argv[i]);
        } else if(strcmp(argv[i], "-mp") == 0) {
            i++;
            config.nmthreads = atoi(argv[i]);
        } else if(strcmp(argv[i], "-o") == 0) {
            i++;
            config.odir = strdup(argv[i]);
        } else if(strcmp(argv[i], "-tmpdir") == 0) {
            i++;
            config.tmpdir= argv[i];
        } else if(strcmp(argv[i], "-upto") == 0) {
            i++;
            //upto = atoi(argv[i]);
            upto = strtoul(argv[i], NULL, 10);
        } else if(strcmp(argv[i], "--directional") == 0) {
            config.directional = 1;
        } else if(strcmp(argv[i], "--unmapped") == 0) {
            config.unmapped = 1;
        } else if(strcmp(argv[i], "--reorder") == 0) {
            config.reorder = 1;
        } else if(strcmp(argv[i], "-@") == 0) {
            config.n_compression_threads = atoi(argv[++i]);
#ifndef NOTHROTTLE
        } else if(strcmp(argv[i], "-queue_size") == 0) {
            i++;
            config.reads_in_queue = atoi(argv[i]);
#endif
#ifdef DEBUG
        } else if(strcmp(argv[i], "-taskid") == 0) {
            i++;
            global_debug_taskid = atoi(argv[i]);
            taskid = global_debug_taskid;
#endif
        } else if(strcmp(argv[i], "--genome-size") == 0) {
            i++;
            chromosomes.max_genome = strtoull(argv[i], NULL, 10);
        } else if(strcmp(argv[i], "--score-min") == 0) {
            i++;
            if(!config.quiet) printf("Changing --score-min from '%c,%f,%f' to %s!\n", config.scoremin_type, config.scoremin_intercept, config.scoremin_coef, argv[i]);
            config.scoremin_type = strtok(argv[i], ",")[0];
            config.scoremin_intercept = (float) atof(strtok(NULL, ","));
            config.scoremin_coef = (float) atof(strtok(NULL, ","));
        } else {
            if(strcmp(argv[i], "--local") == 0 || strcmp(argv[i], "--very-fast-local") == 0 || strcmp(argv[i], "--fast-local") == 0 || strcmp(argv[i], "--sensitive-local") == 0 || strcmp(argv[i], "--very-sensitive-local") == 0) {
                config.mode = 1;
                if(config.scoremin_type == 'L' && config.scoremin_intercept == -0.6f && config.scoremin_coef == -0.6f) {
                    config.scoremin_type = 'G';
                    config.scoremin_intercept = 20.0;
                    config.scoremin_coef = 8.0;
                    if(!config.quiet) printf("Since --local was specified and --score-min was not already changed, changing --score-min to the bowtie2 default of 'G,20,8' (specify --score-min to change this)\n");
                }
            }
            if(strcmp(argv[i], "--quiet") == 0) config.quiet = 1; //This also needs to be passed
            //bowtie2 option
            if(strlen(config.bowtie2_options) + 1 + strlen(argv[i]) >= bowtie2_options_max) {
                bowtie2_options_max = strlen(config.bowtie2_options) + 1 + strlen(argv[i]) + 100;
                config.bowtie2_options = realloc(config.bowtie2_options, sizeof(char) * bowtie2_options_max);
            }
            strcat(config.bowtie2_options, " ");
            strcat(config.bowtie2_options, argv[i]);
        }
    }

    if(config.FASTQ1 == NULL || config.genome_dir == NULL || (config.FASTQ2 == NULL && config.paired == 1)) {
        if(taskid == MASTER) {
            printf("No FASTQ files!\n");
            usage(argv[0]);
        }
        quit(0, -1);
    }

    //If more than one input file was specified, enable reorder
    tmp = strdup(config.FASTQ1);
    p = strtok(tmp, ",");
    if(wordexp(p, &p_wordexp, WRDE_SHOWERR | WRDE_UNDEF) != 0) {
        printf("There was an error while parsing %s.\n", p);
        free(tmp);
        wordfree(&p_wordexp);
        quit(0, -1);
    }
    multi_file += p_wordexp.we_wordc;
    p = strtok(NULL, ",");
    while(p != NULL) {
        if(wordexp(p, &p_wordexp, WRDE_SHOWERR | WRDE_UNDEF | WRDE_REUSE) != 0) {
            printf("There was an error while parsing %s.\n", p);
            free(tmp);
            wordfree(&p_wordexp);
            quit(0, -1);
        }
        multi_file += p_wordexp.we_wordc;
        p = strtok(NULL, ",");
    }
    free(tmp);
    wordfree(&p_wordexp);
    if(multi_file>1) config.reorder=1; //We need to force --reorder if there are multiple input files
#ifdef DEBUG
    if(multi_file>1) {
        printf("In DEBUG mode, you can't input multiple file-sets!\n");
        quit(0,-1);
    }
#else
    if(!config.quiet) printf("%s has rank %i\n", processor_name, taskid); fflush(stdout);
    if(taskid > effective_nodes()) {
        printf("From node %i: So long and thanks for all the bits.\n", taskid); fflush(stdout);
        return(-2); //We're an extraneous node
    }
#endif
    ngroups = effective_nodes();

    if(config.tmpdir == NULL) config.tmpdir = "/tmp";

    //Allocate room for the genome, if needed
    if(taskid == MASTER) {
        if(!config.quiet) printf("Allocating space for %llu characters\n", chromosomes.max_genome);
        fflush(stdout);
        chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
        *chromosomes.genome = '\0';
        if(chromosomes.genome == NULL) {
            printf("Could not allocate enough room to hold the genome!\n");
            return -1;
        }
    } else {
        chromosomes.max_genome = 0;
    }

    //Setup the global variables (these will need to be free()d!)
#ifndef DEBUG
    if(taskid == MASTER) {
#endif
        nwritten = calloc(multi_file, sizeof(char *));
        fnames1 = calloc(multi_file, sizeof(char *));
        fnames2 = calloc(multi_file, sizeof(char *));
        flengths = calloc(multi_file, sizeof(unsigned long long));
#ifndef DEBUG
    }
#endif

    //Append score_min, and p
    if(strlen(config.bowtie2_options) + 1000 >= bowtie2_options_max) {
        bowtie2_options_max = strlen(config.bowtie2_options) + 1000; //This should suffice
        config.bowtie2_options = realloc(config.bowtie2_options, sizeof(char) * bowtie2_options_max);
    }
    if(strlen(config.bowtie2_options) > 0) {
        sprintf(config.bowtie2_options, "%s -p %i --score-min '%c,%g,%g'", config.bowtie2_options, config.nthreads, config.scoremin_type, config.scoremin_intercept, config.scoremin_coef);
    } else {
        sprintf(config.bowtie2_options, "-p %i --score-min '%c,%g,%g'", config.nthreads, config.scoremin_type, config.scoremin_intercept, config.scoremin_coef);
    }

    //There should be as many tasks according to MPI as dictated by the library type.
    if(config.directional) {
        ngroups /= 2;
    } else {
        ngroups /= 4;
    }
    if(ngroups < 1) {
        if(taskid == MASTER) printf("There are only %i groups of nodes available!! You need to allocate more nodes (at least 3 for direcional and 5 for non-directional libraries)!\n", ngroups);
        quit(0, -1);
    }
    //Yes, these silently change user input
    if(config.nmthreads < 1) config.nmthreads = 1;
    if(config.nmthreads > ngroups) config.nmthreads = ngroups;

#ifdef DEBUG
    //DEBUG can't handle multiple files
    update_odir();
    config.basename = get_basename(config.FASTQ1);
    config.outname = malloc(sizeof(char)*(strlen(config.odir)+ strlen(config.basename)+5));
    sprintf(config.outname, "%s%s.bam", config.odir, config.basename);
    if(taskid == MASTER) {
#else
    if(taskid == MASTER) {
        //Deal with the output directory
        update_odir();
#endif

        //Store the genome into memory
        read_genome();

        //Setup the mutexes
        pthread_mutex_init(&metrics_mutex, NULL);

        //Setup the linked-lists
        nodes = malloc(sizeof(struct packed_struct *)*effective_nodes());
        last_sentinel_node = malloc(sizeof(struct packed_struct *)*effective_nodes());
        fastq_nodes = malloc(sizeof(struct packed_struct *)*ngroups);
        last_fastq_sentinel_node = malloc(sizeof(struct packed_struct *)*ngroups);
        to_write_node = malloc(sizeof(struct packed_struct *)*config.nmthreads);
        to_write_sentinel_node = malloc(sizeof(struct packed_struct *)*config.nmthreads);
        for(i=0; i<effective_nodes(); i++) {
            nodes[i] = initialize_list(nodes[i]);
            last_sentinel_node[i] = nodes[i]->next;
        }
        for(i=0; i<ngroups; i++) {
            fastq_nodes[i] = initialize_list(fastq_nodes[i]);
            last_fastq_sentinel_node[i] = fastq_nodes[i]->next;
        }
        for(i=0; i<config.nmthreads; i++) {
            to_write_node[i] = initialize_list(to_write_node[i]);
            to_write_sentinel_node[i] = to_write_node[i]->next;
        }

        //Start the master node processer threads
        threads = calloc(2+config.nmthreads, sizeof(pthread_t));
        int *threadids = malloc(sizeof(int)*config.nmthreads);
        pthread_create(&(threads[0]), NULL, &send_store_fastq, (void *) &upto);
        for(i=0; i<config.nmthreads; i++) {
            *(threadids+i) = i;
            pthread_create(&(threads[i+1]), NULL, &herd_master_processer_thread, threadids+i);
        }
        pthread_create(&(threads[1+config.nmthreads]), NULL, &bam_writer, NULL);
        herd_slurp(NULL);
        pthread_join(threads[0], NULL);
        for(i=0; i<config.nmthreads; i++) pthread_join(threads[i+1], NULL);
        pthread_join(threads[1+config.nmthreads], NULL);

        //Start freeing things up
        free(threadids);
        free(threads);
        for(i=0; i<effective_nodes(); i++) destroy_list(nodes[i]);
        free(nodes);
        free(last_sentinel_node);
        for(i=0; i<ngroups; i++) destroy_raw_list(fastq_nodes[i]);
        free(fastq_nodes);
        free(last_fastq_sentinel_node);
        free(to_write_node);
        free(to_write_sentinel_node);
        pthread_mutex_destroy(&metrics_mutex);

        //Print some metrics
        bam_header_destroy(global_header);
    } else {
        //Create a temporary directory
        char *tmpdir = malloc(sizeof(char) * (strlen(config.tmpdir) + strlen("/herd_XXXXXX") + 1));
        sprintf(tmpdir, "%s/herd_XXXXXX", config.tmpdir);
        tmpdir = mkdtemp(tmpdir);

        //Name the FIFOs
        slurp_fastq_struct *silly_struct = malloc(sizeof(slurp_fastq_struct));
        silly_struct->thread_id = taskid;
        silly_struct->fastq1 = malloc(sizeof(char) * (strlen(tmpdir) + strlen("/read1") + 1));
        silly_struct->fastq2 = malloc(sizeof(char) * (strlen(tmpdir) + strlen("/read2") + 1));
        sprintf(silly_struct->fastq1, "%s/read1", tmpdir);
        sprintf(silly_struct->fastq2, "%s/read2", tmpdir);

        mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;
        int rv = mkfifo(silly_struct->fastq1, mode);
        if(rv != 0) {
            printf("mkfifo returned with status %i!\n", rv);
            fflush(stdout);
        }
        if(config.paired) {
            rv = mkfifo(silly_struct->fastq2, mode);
            if(rv != 0) {
                printf("mkfifo returned with status %i!\n", rv);
                fflush(stdout);
            }
        }

        //Start slurping in the fastq reads and converting them so they can be aligned
#ifndef DEBUG
        threads = calloc(1, sizeof(pthread_t));
        pthread_create(&(threads[0]), NULL, &slurp_fastq, (void *) silly_struct);
#else
        threads = calloc(2, sizeof(pthread_t));
        pthread_create(&(threads[1]), NULL, &send_store_fastq, (void *) &upto);
        pthread_create(&(threads[0]), NULL, &slurp_fastq, (void *) silly_struct);
#endif
        //worker node stuff
        herd_worker_node(taskid, silly_struct->fastq1, silly_struct->fastq2);
        pthread_join(threads[0], NULL);
#ifdef DEBUG
        pthread_join(threads[1], NULL);
#endif
        if(!config.quiet) printf("Returning from worker node %i\n", taskid);
        fflush(stdout);
        free(silly_struct->fastq1); //The worker node unlinks this
        free(silly_struct->fastq2); //The worker node unlinks this
        free(silly_struct);
        if(rmdir(tmpdir) != 0) {
            printf("Couldn't remove %s directory!\n", tmpdir);
            fflush(stdout);
        }
        free(tmpdir);
        free(threads);
    }

#ifndef DEBUG
    if(taskid == MASTER) {
#endif
        free(nwritten);
        for(i=0;i<multi_file;i++) free(fnames1[i]);
        free(fnames1);
        if(config.paired) for(i=0;i<multi_file;i++) free(fnames2[i]);
        free(fnames2);
        free(flengths);
#ifndef DEBUG
    }
#endif

    //Clean up
    if(config.odir != NULL) free(config.odir);
    if(config.genome_dir != NULL) free(config.genome_dir);
    quit(3, 0);
    return 0;
}
