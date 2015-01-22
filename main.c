#include "bison.h"

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] -g genome_dir {-1 fastq.gz -2 fastq.gz | -U fastq.gz}\n", prog);
    printf("\n\
    N.B., Any option not listed below will be passed directly to bowtie2, so you\n\
    can specify, e.g., --very-fast if you want. If you specify --local,\n\
    --score-min is changed back to the bowtie2 default of 'G,20,6', unless you\n\
    specify otherwise.\n\
\n\
-g          Directory containing the genome fasta files and the\n\
            Bisulfite_Sequences directory.\n\
\n\
-1          Fastq file containing read #1 (normally named something like \n\
            foo_1.fastq.gz). Reads needn't be gzipped, but that'll be more\n\
            convenient.\n\
\n\
-2          As with -1, but with read #2.\n\
\n\
-U          For convenience, this denotes a fastq file from single-ended reads.\n\
            Alternatively, -1 can be used without using -2.\n\
\n\
-p          How many threads bowtie2 should use on each node. Default is 12.\n\
\n\
-C          Output CRAM rather than BAM.\n\
\n\
--sort      Sort the output by coordinate. Note that bison_markduplicates and\n\
            bison_methylation_extractor are not able to process coordinate\n\
            sorted files.\n\
\n\
-m          The amount of memory to use per-thread for sorting (only used if\n\
            --sort is specified). The formatting is the same as for\n\
            'samtools sort', so 2G specifies 2 gigabytes and 1K one kilobyte.\n\
            Accepted suffixes are K, M, and G. Note that without one of these\n\
            suffixes, the value is assumed to specify the number of bytes. The\n\
            default is 512M.\n\
\n\
-o          Output directory. By default, everything will be written to the\n\
            directory holding the fastq files (or the file containing read #1,\n\
            as appropriate). If you would prefer for the output BAM file and\n\
            metrics txt file to be placed elsewhere, specify that here.\n\
\n\
            N.B., the directory must exist! \n\
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
--no-discordant Suppress discordant alignments. This is actually a bowtie2\n\
            option.\n\
\n\
--no-mixed  Suppress singleton alignments. This is actually a bowtie2 option.\n\
\n\
--unmapped  Save unaligned reads to a file or files (as appropriate).\n\
\n\
--genome-size Many of the bison tools need to read the genome into memory. By\n\
            default, they allocate 3000000000 bases worth of memory for this and\n\
            increase that as needed. However, this can sometimes be far more\n\
            than is needed (meaning wasted memory) or far too little (in which\n\
            case the process can become quite slow). If you input the\n\
            approximate size of your genome here (in bases), then you can\n\
            maximize performance and minimize wasted space. It's convenient to\n\
            round up a little.\n\
\n\
--quiet     Suppress printing of anything other than errors to the console.\n\
\n\
-h          Print this help message.\n\
\n\
-v          Print version information.\n\
\n");
#ifdef DEBUG
    printf("\n\
-taskid     Which node number to act as. The default is 0, the master node.\n\
            Other possibilities are 1-4, which are the worker nodes that\n\
            process OT, OB, CTOT, and CTOB alignments, respectively. A value of\n\
            -1 will only convert the reads.\n\
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
    //These are bowtie2 options, but they'll just be pushed into the config.bowtie2_options character array
    char *cmd, *p;
    unsigned long upto = 0;
#ifndef DEBUG
    int mpi_ntasks, ntasks;
    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif

    //Deal with MPI initialization, this seems like an odd way to do things.
#ifndef DEBUG
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(processor_name, &name_len);
#else
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
#endif
    if(provided != MPI_THREAD_FUNNELED) {
        fprintf(stderr, "Your implementation does not support MPI_THREAD_FUNNELED, which is required for bison to run. This is actually quite unusual!\n");
        return -1;
    }

    config.odir = NULL;
    config.paired = 0; //Default is single-ended
    config.directional = 0; //Default is non-directional
    config.nthreads = 12; //Default is 12 threads/node
    config.bowtie2_options = calloc(MAXREAD, sizeof(char));
    config.unmapped = 0; //By default, unmapped reads are NOT written to a fastq file
    config.scoremin_type = 'L'; //--score-min 'L,-0.6,-0.6'
    config.scoremin_intercept = -0.6;
    config.scoremin_coef = -0.6;
    config.mode = 0; //--end-to-end
    config.quiet = 0;
    config.FASTQ1 = NULL;
    config.FASTQ2 = NULL;
    config.genome_dir = NULL;
    config.isCRAM = 0;
    config.fai = NULL;
    config.argc = argc;
    config.argv = argv;
    config.maxMem = 512<<20;
    config.sort = 0;
    config.n_compression_threads = 1;
    chromosomes.max_genome = 3000000000;
    chromosomes.nchromosomes = 0; //We need to initialize the struct

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
        } else if(strcmp(argv[i], "-C") == 0) {
            config.isCRAM = 1;
        } else if(strcmp(argv[i], "-m") == 0) {
            config.maxMem = str2Mem(argv[++i]);
        } else if(strcmp(argv[i], "--sort") == 0) {
            config.sort = 1;
        } else if(strcmp(argv[i], "-o") == 0) {
            i++;
            config.odir = strdup(argv[i]);
        } else if(strcmp(argv[i], "-upto") == 0) {
            i++;
            upto = strtoul(argv[i], NULL, 10);
        } else if(strcmp(argv[i], "--directional") == 0) {
            config.directional = 1;
        } else if(strcmp(argv[i], "--unmapped") == 0) {
            config.unmapped = 1;
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
            if(!config.quiet) fprintf(stderr, "Changing --score-min from '%c,%f,%f' to %s!\n", config.scoremin_type, config.scoremin_intercept, config.scoremin_coef, argv[i]);
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
                    if(!config.quiet) fprintf(stderr, "Since --local was specified and --score-min was not already changed, changing --score-min to the bowtie2 default of 'G,20,8' (specify --score-min to change this)\n");
                }
            }
            if(strcmp(argv[i], "--quiet") == 0) config.quiet = 1;
            //bowtie2 option
            if(strlen(config.bowtie2_options) + 1 + strlen(argv[i]) >= bowtie2_options_max) {
                bowtie2_options_max = strlen(config.bowtie2_options) + 1 + strlen(argv[i]) + 100;
                config.bowtie2_options = realloc(config.bowtie2_options, sizeof(char) * bowtie2_options_max);
            }
            strcat(config.bowtie2_options, " ");
            strcat(config.bowtie2_options, argv[i]);
        }
    }

#ifndef DEBUG
    if(!config.quiet) {printf("%s has rank %i\n", processor_name, taskid); fflush(stderr);}
#endif

    if(config.FASTQ1 == NULL || config.genome_dir == NULL || (config.FASTQ2 == NULL && config.paired == 1)) {
        if(taskid == MASTER) {
            fprintf(stderr, "No FASTQ files!\n");
            usage(argv[0]);
        }
        quit(0, -1);
    }

    //Allocate room for the genome, if needed
    if(taskid == MASTER) {
        if(!config.quiet) fprintf(stderr, "Allocating space for %llu characters\n", chromosomes.max_genome);
        fflush(stderr);
        chromosomes.genome = malloc(sizeof(char)*chromosomes.max_genome);
        *chromosomes.genome = '\0';
        if(chromosomes.genome == NULL) {
            fprintf(stderr, "Could not allocate enough room to hold the genome!\n");
            return -1;
        }
    } else {
        chromosomes.max_genome = 0;
    }

    //Append score_min, and p
    if(strlen(config.bowtie2_options) + 1000 >= bowtie2_options_max) {
        bowtie2_options_max = strlen(config.bowtie2_options) + 1000; //This should suffice
        config.bowtie2_options = realloc(config.bowtie2_options, sizeof(char) * bowtie2_options_max);
    }
    sprintf(config.bowtie2_options, "%s -p %i --score-min '%c,%g,%g'", config.bowtie2_options, config.nthreads, config.scoremin_type, config.scoremin_intercept, config.scoremin_coef);

    //There should be as many tasks according to MPI as dictated by the library type.
#ifndef DEBUG
    ntasks = 5;
    if(config.directional) ntasks = 3;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_ntasks);
    if(mpi_ntasks < ntasks) {
        if(taskid == MASTER) fprintf(stderr, "There are only %i nodes available but we need %i to work. You need to allocate more nodes!\n", mpi_ntasks, ntasks);
        quit(0, -1);
    }
#endif
    /***********************************************************************************************
    *
    * Convert the input reads C->T and G->A as needed. There are 4 use cases:
    *   Directional:
    *       Paired-end: FASTQ1 will be C->T converted and FASTQ2 will be G->A converted.
    *       Single-end: FASTQ1 will be C->T converted
    *   Non-directional:
    *       Paired-end: Both FASTQ1 and FASTQ2 will be C->T and G->A converted.
    *       Single-end: FASTQ1 will be both C->T and G->A converted.
    *
    *   convert_fastq() takes a single integer parameter:
    *       8 = convert FASTQ1 C->T
    *       4 = convert FASTQ1 G->A
    *       2 = convert FASTQ2 C->T
    *       1 = convert FASTQ2 G->A
    *
    ***********************************************************************************************/
    update_odir();
    create_fastq_names(config.FASTQ1, config.FASTQ2);

#ifndef DEBUG
    if(taskid == MASTER) {
#endif
        //MASTER specific procedures
        config.basename = get_basename(config.FASTQ1);
        config.outname = malloc(sizeof(char)*(strlen(config.odir)+ strlen(config.basename)+6));
        if(config.isCRAM) sprintf(config.outname, "%s%s.cram", config.odir, config.basename);
        else sprintf(config.outname, "%s%s.bam", config.odir, config.basename);
#ifdef DEBUG
        //When debugging, don't convert the files if it's already been done
        if(access(config.FASTQ1CT, F_OK) == -1) {
#endif
        if(config.directional) {
            if(config.paired) {
                convert_fastq(9, upto);
            } else {
                convert_fastq(8, upto);
            }
        } else {
            if(config.paired) {
                convert_fastq(15, upto);
            } else {
                convert_fastq(12, upto);
            }
        }
#ifdef DEBUG
    }
    //Just convert the reads
    if(taskid == -1)  {
        quit(3, 0);
        return 0;
    }
#endif

#ifdef DEBUG
    if(taskid == MASTER) {
#endif
        //Open the input reads
        cmd = malloc(sizeof(char) * (strlen(config.FASTQ1) + 7));
        p = strrchr(config.FASTQ1, '.');
        if(strcmp(p,".gz") == 0 || strcmp(p,".GZ") == 0) {
            sprintf(cmd, "zcat %s", config.FASTQ1);
        } else if(strcmp(p,".bz") == 0 || strcmp(p,".bz2") == 0) {
            sprintf(cmd, "bzcat %s", config.FASTQ1);
        } else {
            sprintf(cmd, "cat %s", config.FASTQ1);
        }
        zip1 = popen(cmd, "r");
        if(config.paired) {
            cmd = realloc(cmd, sizeof(char) * (strlen(config.FASTQ2) + 7));
            p = strrchr(config.FASTQ2, '.');
            if(strcmp(p,".gz") == 0 || strcmp(p,".GZ") == 0) {
                sprintf(cmd, "zcat %s", config.FASTQ2);
            } else if(strcmp(p,".bz") == 0 || strcmp(p,".bz2") == 0) {
                sprintf(cmd, "bzcat %s", config.FASTQ2);
            } else {
                sprintf(cmd, "cat %s", config.FASTQ2);
            }
            zip2 = popen(cmd, "r");
        }

        //Open the output file handles
        if(config.unmapped) {
            cmd = realloc(cmd, sizeof(char) * (strlen(config.unmapped1) + 8));
            if(!config.quiet) fprintf(stderr, "Writing unmapped reads to %s\n", config.unmapped1);
            sprintf(cmd, "gzip > %s", config.unmapped1);
            unmapped1 = popen(cmd, "w");
            if(config.paired) {
                cmd = realloc(cmd, sizeof(char) * (strlen(config.unmapped2) + 8));
                if(!config.quiet) fprintf(stderr, "Writing unmapped reads to %s\n", config.unmapped2);
                sprintf(cmd, "gzip > %s", config.unmapped2);
                unmapped2 = popen(cmd, "w");
            }
        }
        free(cmd);

        //Store the genome into memory
        read_genome();

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
        if(!config.quiet) fprintf(stderr, "Alignment metrics will be printed to %s%s.txt\n",config.odir,config.basename);
        fflush(stderr);

        //Setup the linked-lists
        node1 = initialize_list(node1);
        node1_last_sentinel = node1->next;
        node2 = initialize_list(node2);
        node2_last_sentinel = node2->next;
        node3 = initialize_list(node3);
        node3_last_sentinel = node3->next;
        node4 = initialize_list(node4);
        node4_last_sentinel = node4->next;

        //Start the master node processer threads
        threads = calloc(1, sizeof(pthread_t));
        pthread_create(&(threads[0]), NULL, &master_processer_thread, NULL);
        slurp(NULL);
        pthread_join(threads[0], NULL);

        //Start freeing things up
        if(!config.quiet) fprintf(stderr, "Closing input files\n");
        free(threads);
        pclose(zip1);
        if(config.paired) pclose(zip2);

        //Print some metrics
        print_metrics();
    } else {
        //worker node stuff, wait for the master
        worker_node(taskid);
        if(!config.quiet) fprintf(stderr, "Returning from worker node %i\n", taskid);
        fflush(stderr);
    }

    //Clean up
    if(config.odir != NULL) free(config.odir);
    if(config.genome_dir != NULL) free(config.genome_dir);
    quit(3, 0);
    return 0;
}
