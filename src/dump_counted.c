#include <stdio.h>
#include <unistd.h> /* getopt */
#include <ctype.h>  /* getopt */

#include "kmers.h"

int BASES[256];

void getKmerChar(uint64_t kmer_uint, int kmer_length, char *kmer_char) {
    int i;
    for (i = 0; i < kmer_length; i++) {
        kmer_char[i] = BASES[ (kmer_uint >> ((kmer_length * 2) - i * 2)) & ((uint64_t)3) ];
    }
}

int main(int argc, char *argv[])
{
    BASES[0] = 'A';    BASES[1] = 'C';    BASES[2] = 'G';    BASES[3] = 'T';
    int kmer_length = 31;
    int pairs_in_once = 1000000;

    char in_path[1024] = "";
    char out_path[1024] = "./";
    char out_filename_base[1024] = "kmers_counted";


    int out_path_set = 0;
    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "t:m:o:k:l:h")) != -1) {
        switch (c) {
        case 'o':
            sprintf(out_path, "%s", optarg);
            out_path_set = 1;
            break;
        case 'h':
            printf("Usage: dump_kmers [options] <input_file_path>\n");
            printf("* - required\n\n");
            printf(" -o=string\t *Output file name (full path)\n");
            return 0;
            break;
        case '?':
            return 1;
        default:
            ;
        }
    }

    if (out_path_set == 0) {
        fprintf(stderr, "Error: Output file name (full path) is required (-o <path>)\n");
        exit(1);
    }

    sprintf(in_path, "%s", argv[optind]);
    FILE *input_file = fopen(in_path, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Error: input file missing\n");
        exit(1);
    }

    FILE *fr_kmers_counted = fopen(in_path, "r");
    FILE *fw = fopen(out_path, "w");

    KmerCountPair *kcp = (KmerCountPair*)malloc(sizeof(KmerCountPair) * pairs_in_once);
    char *kmer_char = (char *) malloc(sizeof(char) * kmer_length);
    size_t r, i;
    while ((r = fread(kcp, sizeof(KmerCountPair), pairs_in_once, fr_kmers_counted)) > 0) {

        for (i = 0; i < r; i++) {
            getKmerChar(kcp[i].kmer, kmer_length, kmer_char);
            // printf("%llu => %llu | %s => %llu\n", kcp[i].kmer, kcp[i].count, kmer_char, kcp[i].count);
            fprintf(fw, "%s\t%lu\n", kmer_char, kcp[i].count);
        }
    }
    free(kmer_char);
    free(kcp);
    
    fclose(fw);
    fclose(fr_kmers_counted);
}
    