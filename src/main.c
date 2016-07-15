#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>

#include "defs.h"
#include "sort.h"
#include "kmers.h"
#include "partition.h"

void createSortSave(char *buffer, ThreadPart *threadParts, size_t memory_for_kmers, int kmer_length, char out_path[1024], int myid, int core_part_num, int nthreads, int tid, int nmers_int);

int main(int argc, char *argv[]) 
{    
    struct timeval t_start, t_end;
    struct timeval t_start_createSortSave, t_end_createSortSave;
    gettimeofday(&t_start, NULL);
    
    int i, j, c, t;
    int nps, myid;

    // MPI Init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nps);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    int ARG_CORE_THREADS = 2;
    int ARG_CORE_MEMORY = 2048; 
    int kmer_length = 31;
    int nmers_bits = 16;
    int nmers_int = (int)pow((double)2, (double)nmers_bits);

    char in_path[1024] = "/home/3xkuban/work/freqCountData/data/SRX040485/raw_part_aa";
    char out_path[1024] = "/home/3xkuban/work/freqCountData/outs/parker/";
    char out_filename_base[1024] = "kmers_counted";

    int min_threshold = 2; // min amount o k-mers to be saved

    int out_path_set = 0;
    int kmer_length_set = 0;
    opterr = 0;
    while ((c = getopt (argc, argv, "t:m:o:k:l:h")) != -1) {
        switch (c) {
        case 't':
            ARG_CORE_THREADS = atoi(optarg);
            break;
        case 'm':
            ARG_CORE_MEMORY = atoi(optarg);
            break;
        case 'o':
            sprintf(out_path, "%s", optarg);
            out_path_set = 1;
            break;
        case 'k':
            kmer_length = atoi(optarg);
            kmer_length_set = 1;
            break;
        case 'l':
            min_threshold = atoi(optarg);
            break;
        case 'h':
            if (myid == 0) {
                printf("Usage: mpirun -n <num_of_processors> parker [options] <input_file_path>\n");
                printf("Options (default value in (), * - required):\n\n");
                printf(" -k=uint32\t *Length of k-mer (31)\n");
                printf(" -t=uint32\t Number of threads (2)\n");
                printf(" -m=uint32\t Available memory for core/processor in MB (2048)\n");
                printf(" -o=string\t *Output file name (full path)\n");
                printf(" -l=uint32\t Don't output k-mers with count lower than this value (2)\n");
            }
            MPI_Finalize();
            return 0;
            break;
        case '?':
        //     if (optopt == 'c')
        //     fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        //     else if (isprint (optopt))
        //     fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        //     else
        //     fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            ;
        }
    }
    // k-mer length vs nmers_int
    // TODO: checkDirectoryExists
    // TODO: checkFileExists
    // TODO: fixPath ? => add slash /

    if (kmer_length_set == 0) {
        fprintf(stderr, "Error: Kmer length is required (-k <value>)\n");
        exit(1);
    }
    if (out_path_set == 0) {
        fprintf(stderr, "Error: Output file name (full path) is required (-o <path>)\n");
        exit(1);
    }
    // nastavenie in_path    
    sprintf(in_path, "%s", argv[optind]);
    FILE *input_file = fopen(in_path, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Error: input file missing\n");
        exit(1);
    }
    size_t input_file_size = getFileSize(input_file);

    size_t MEMORY_MB = ARG_CORE_MEMORY;
    MEMORY_MB -= 16; // nejaka rezerva // debug uncomment
    int NUM_THREADS_PER_CORE = ARG_CORE_THREADS;

    // core required memories
    size_t max_memory = (size_t)MEMORY_MB*1024*1024;
    size_t memory_for_kmers = max_memory / 2;
    size_t memory_for_raw_part = memory_for_kmers / 8;
    
    FilePart *fileParts = (FilePart*)malloc(sizeof(FilePart) * nps);
    getFilePartsStartsAndEnds(input_file, input_file_size, END_READ_RESERVE, fileParts, nps, memory_for_kmers, memory_for_raw_part);
    
    fclose(input_file);
    
    int nthreads = NUM_THREADS_PER_CORE;
    ThreadPart *threadParts;
    char *buffer = (char*)malloc(sizeof(char) * memory_for_raw_part);

    #pragma omp parallel num_threads(nthreads) shared(nthreads, threadParts, buffer, in_path) private(i)
    {
        nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        
        FILE *input_file = fopen(in_path, "r");
        
        for (i = 0; i < fileParts[myid].num_of_core_parts; i++) {
            CorePart corePart = fileParts[myid].coreParts[i];
            threadParts = (ThreadPart*)malloc(sizeof(ThreadPart) * nthreads);  
            #pragma omp master
            {
                fseek(input_file, corePart.start_index, SEEK_SET);
                size_t buf_size = fread(buffer, sizeof(char), corePart.size, input_file);
                getThreadPartsStartsAndEnds(buffer, corePart, END_READ_RESERVE, nthreads, threadParts);
            }
            #pragma omp barrier
            createSortSave(buffer, threadParts, memory_for_kmers, kmer_length, out_path, myid, i, nthreads, tid, nmers_int);
            #pragma omp barrier //
        }
    }
    
    #pragma omp barrier // nanic?
    MPI_Barrier(MPI_COMM_WORLD);
    free(buffer);
    free(threadParts);
    
    int **nmers_parts = (int**)malloc(sizeof(int*) * nps * nthreads * fileParts[myid].num_of_core_parts);     // zapamatat si pre jednotlive casti (Part <=> subory)
    // printf("|nmers_parts| = %d\n", nps * nthreads * coreParts[myid].num_of_parts);
    int *nmers_total = (int*)malloc(sizeof(int) * nmers_int);  // zapamatat si kolko je vo vsetkych partoch (suboroch)
    for (i = 0; i < nmers_int; i++) {nmers_total[i] = 0;}
    
    for (i = 0; i < fileParts[myid].num_of_core_parts * nps * nthreads; i++) {   // TODO: skontrolovat
        nmers_parts[i] = (int*)malloc(sizeof(int) * nmers_int);
        for (j = 0; j < nmers_int; j++) {
            nmers_parts[i][j] = 0;
        }
    }
    
    // prejst subory s nmermi
    size_t total_kmers_count = 0;
    int *nmers = (int*)malloc(sizeof(int) * nmers_int);
    for (c = 0; c < nps; c++) {
        for (t = 0; t < nthreads * fileParts[c].num_of_core_parts; t++) {
            char filename_nmers[1024];
            sprintf(filename_nmers, "%snmers_%d_%d", out_path, c, t);

            FILE *fr_nmers = fopen(filename_nmers, "r");
            
            if (fr_nmers == NULL) {
                fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_nmers);
                exit(1);
            }
            fread(nmers, sizeof(int), nmers_int, fr_nmers);
            fclose(fr_nmers);
            
            for (i = 0; i < nmers_int; i++) {
                nmers_parts[c * nthreads * fileParts[c].num_of_core_parts + t][i] = nmers[i];
                nmers_total[i] += nmers[i];
                total_kmers_count += nmers[i];
            }
        }
    }

    // COUNT PARTS
    size_t max_kmers_per_count_part = (memory_for_kmers * nthreads) / sizeof(uint64_t);
    size_t count_parts_cca = (total_kmers_count / max_kmers_per_count_part) + 1;
    
    // count parts count dividable by nps ! 
    size_t count_parts_tmp = 0;
    while (count_parts_tmp < count_parts_cca) {
        count_parts_tmp += nps;
    }
    count_parts_cca = count_parts_tmp;
    size_t balanced_kmers_per_count_part = total_kmers_count / count_parts_cca;
    // postupne pridavat 1% k balanced kmers poctu
    size_t one_percent_kmers_count = (max_kmers_per_count_part - balanced_kmers_per_count_part) / 100;
    
    while (balanced_kmers_per_count_part * count_parts_cca < total_kmers_count) {
        balanced_kmers_per_count_part += one_percent_kmers_count;
    }
    CountPart *countParts = (CountPart*)malloc(sizeof(CountPart) * count_parts_cca * 2); // aj nejaka rezerva
    
    int count_part_id = 0;
    size_t count_kmers = 0;
    size_t total_count_kmers = 0;

    countParts[count_part_id].id = count_part_id;
    countParts[count_part_id].start_nmer = 0;
    countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
    countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    
    for (i = 0; i < nmers_int; i++) {

        if (i == nmers_int - 1) { // ked posledny n-mer je moc velky 
            count_kmers += nmers_total[i];
            total_count_kmers += count_kmers;
            countParts[count_part_id].id = count_part_id;
            countParts[count_part_id].kmers_count = count_kmers;
            countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
            countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_nmer = i;
        } else if (count_kmers + nmers_total[i] > balanced_kmers_per_count_part/* || i == nmers_int - 1 */) {
            total_count_kmers += count_kmers;
            countParts[count_part_id].id = count_part_id;
            countParts[count_part_id].kmers_count = count_kmers;
            countParts[count_part_id].end_nmer = i - 1;
            countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
            countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            
            count_part_id++;
            count_kmers = nmers_total[i];
            countParts[count_part_id].start_nmer = i;
        } else {
            count_kmers += nmers_total[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // 1.6
    // [CountedParams init]
    CountedParams countedParams;
    countedParams.nmers_bits = nmers_bits;
    countedParams.counted_parts_count = count_part_id + 1;
    countedParams.start_nmers = (int*)malloc(sizeof(int) * (count_part_id + 1));
    countedParams.end_nmers = (int*)malloc(sizeof(int) * (count_part_id + 1));
    // [END CountedParams init]
    
    for (i = 0; i <= count_part_id; i++) {
        CountPart cp = countParts[i];
        size_t offset = 0;

        countedParams.start_nmers[i] = cp.start_nmer;
        countedParams.end_nmers[i] = cp.end_nmer;
        
        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[myid].num_of_core_parts;
            for (t = 0; t < num_of_core_parts * nthreads; t++) {
                char filename_kmers[1024];
                sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
                FILE *fr_kmers = fopen(filename_kmers, "r");
                
                if (fr_kmers == NULL) {
                    fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_kmers);
                    exit(1);
                }
                // printf("fileSize of '%s' is %lu\n", filename_kmers, getFileSize(fr_kmers));
                // fprintf(fdebug, "fileSize of '%s' is %lu\n", filename_kmers, getFileSize(fr_kmers));
                // countParts[i].files_size[c * num_of_core_parts * nthreads + t] = getFileSize(fr_kmers);
                fclose(fr_kmers);
            }
        }
    }
    
    if (myid == 0) {
        saveCountedParams(countedParams, out_path);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // finding cummulative counts of nmers in each file (where those nmers are, start_offsets and end_offsets)
    
    int count_parts_count = count_part_id + 1;
    for (i = 0; i < count_parts_count; i++) {
        CountPart cp = countParts[i];
        
        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[c].num_of_core_parts;
            for (t = 0; t < nthreads * num_of_core_parts; t++) {
                size_t start_nmer = cp.start_nmer;
                size_t end_nmer = cp.end_nmer;

                cp.start_offsets[c * nthreads * num_of_core_parts + t] = 0;
                cp.end_offsets[c * nthreads * num_of_core_parts + t] = 0;
                int n;

                for (n = 0; n <= end_nmer; n++) {
                    
                    if (n < start_nmer) {
                        cp.start_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                        cp.end_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                    } else if (n <= end_nmer) {
                        cp.end_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                    }
                }
            }
        }   
    }
    
    uint64_t *kmers = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));
    uint64_t *aux = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));
    for (i = 0; i < count_parts_count; i++) {
//     	printf("[%d]: count_part_id = %d\n", myid, i);
        int cp_id = i;
        if (cp_id % nps != myid) { // kazdy svoju cast
            continue;            
        }
        aux = (uint64_t*)aux;
        CountPart cp = countParts[cp_id];
        size_t offset = 0;

        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[myid].num_of_core_parts;
            for (t = 0; t < num_of_core_parts * nthreads; t++) {
                char filename_kmers[1024];
                sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
                FILE *fr_kmers = fopen(filename_kmers, "r");
                
                if (fr_kmers == NULL) {
                    fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_kmers);
                    exit(1);
                }
                int offset_id = c * num_of_core_parts * nthreads + t; // 11.4
                fseek(fr_kmers, cp.start_offsets[offset_id] * sizeof(uint64_t), SEEK_SET);
                size_t count_to_read =  cp.end_offsets[offset_id] - cp.start_offsets[offset_id];
                size_t really_read_kmers = fread(&kmers[offset], sizeof(uint64_t), count_to_read, fr_kmers);
				offset += really_read_kmers;
                fclose(fr_kmers);
            }
        }

        if (offset > max_kmers_per_count_part) {
            printf("CHYBA offset > max_kmers_per_count_part!\n >>> TODO <<<"); // TODO !!! => s balanced by sa to nemalo stat
            offset = max_kmers_per_count_part;  // TODO: inak (pridat dalsiu cast?)
        }

        bucketSortSerial_uint64(aux, kmers, offset, 0, 12, 1);

//////////
        char filename_kmers_counted[1024];
        sprintf(filename_kmers_counted, "%s%s_%d_%d", out_path, out_filename_base, cp.start_nmer, cp.end_nmer);
        FILE *fw_kmers_counted = fopen(filename_kmers_counted, "wb");

        uint64_t kmer = kmers[0];
        size_t diff_kmers_count = 0;
        size_t total_kmers_counted = 0;
        size_t threshold_kmers_count = 0;
        
        size_t kmer_count = 1;
        size_t max_pairs = (max_kmers_per_count_part * sizeof(uint64_t)) / sizeof(KmerCountPair); // sizeof(aux) / sizeof(KmerCountPair)
//         max_pairs /= 2; // test 4.4
// 				size_t max_pairs = 10000000;

        KmerCountPair *kcp = (KmerCountPair*)aux;
// 				KmerCountPair *kcp = (KmerCountPair*)malloc(sizeof(KmerCountPair) * (max_pairs + 1000));
        size_t pair_count = 0;
        size_t flushed = 0;
                
        size_t ii;
       
        for (ii = 1; ii <= offset; ii++) {

            if (kmers[ii] != kmer) {
                diff_kmers_count++;
                
                if (kmer_count >= min_threshold) {
                    kcp[pair_count].kmer = kmer;
                    kcp[pair_count].count = kmer_count;
                    total_kmers_counted += kmer_count;
                    kmer_count = 0;
                    kmer = kmers[ii];
                    pair_count++;
                    flushed = 0;
                } else {
                    threshold_kmers_count++;
                }
                        
                if (pair_count == max_pairs) {
                    fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);                            
                    pair_count = 0;
                    flushed = 1;
                }             
                // break;
            }
            kmer_count++;
        } 
        // printf("diff_kmers_count = %lu\n", diff_kmers_count);
        // printf("threshold_kmers_count = %lu\n", threshold_kmers_count);

        if (flushed == 0 && pair_count > 0) {
            fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);
        }
        fclose(fw_kmers_counted);
    } // [end] count_parts_count
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // free each corePart in fileParts?
    
    free(fileParts);
    free(countParts);
    free(kmers);
    free(aux);

    // remove tmp files
    for (c = 0; c < nps; c++) {
        for (t = 0; t < nthreads * fileParts[c].num_of_core_parts; t++) {
            char filename_kmers[1024];
            char filename_nmers[1024];
            sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
            sprintf(filename_nmers, "%snmers_%d_%d", out_path, c, t);
            remove(filename_kmers);
            remove(filename_nmers);
        }
    }

    // final time
    gettimeofday(&t_end, NULL);
    double time_total = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    printf("[%d]: total time = %lf\n", myid, time_total);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

void createSortSave(char *buffer, ThreadPart *threadParts, size_t memory_for_kmers, int kmer_length, char out_path[1024], int myid, int core_part_num, int nthreads, int tid, int nmers_int) {
    struct timeval t_start, t_end;
    
    size_t kmers_limit = memory_for_kmers / sizeof(uint64_t);  // max k-mers for core
    size_t kmers_limit_thread = kmers_limit / nthreads;        // max k-mers for thread // TODO: deleno 2 -> pre aux?
    size_t kmers_count = 0;
    uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);
    uint64_t *aux = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);

    int k;
    int *nmers_counts = (int*)malloc(sizeof(int) * nmers_int);
    int *nmers_counted = (int*)malloc(sizeof(int) * nmers_int);
    for (k = 0; k < nmers_int; k++) {nmers_counts[k] = 0; nmers_counted[k] = 0;}

    createKmers(&buffer[threadParts[tid].start_index], threadParts[tid].size, kmer_length, kmers_limit_thread, kmers, &kmers_count, nmers_counts, nmers_int);
    bucketSortSerial_uint64(aux, kmers, kmers_count, 0, 12, 1);
    
    saveKmers(kmers, kmers_count, out_path, myid, core_part_num * nthreads + tid);
    saveNmersCounts(nmers_counts, nmers_int, out_path, myid, core_part_num * nthreads + tid);

    free(nmers_counts);
    free(nmers_counted);

    free(aux);
    free(kmers);
}