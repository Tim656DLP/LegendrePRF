/*
MIT License

Copyright (c) 2019 Novak Kaluđerović

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "gmp.h"
#include <time.h> 
#include <inttypes.h>
#include <string.h>
#include "common.h"
#include <unistd.h>


/* Static parameters (should be read from a project file or input) */
static char* folder;
static int num_threads = 1;



static void usage()
{
    fprintf(stderr, "precomp -b bname [ options ]\n");
    fprintf(stderr, "  -b bname: project file name\n");
    fprintf(stderr, "  -t num_threads: number of threads. default: 1\n");

    fprintf(stderr, "  -h help\n");
    exit(1);
}





static void get_options(int argc, char **argv)
{
    char c;

    folder = NULL;

    while ((c = getopt(argc, argv, "b:t:a:h")) != (char)(-1))
    {
        switch (c)
        {
            case 'b':
                folder = optarg;
                break;
            case 't':
                if (sscanf(optarg, "%d", &num_threads) != 1)
                    complain("Bad argument to -t!\n");
                break;
            case 'h':
                usage();
                break;
            default:
                fprintf(stderr, "Bad option %c\n", (char)c);
                usage();
        }
    }

    if (folder == NULL)
    {
        fprintf(stderr, "argument -b bname is necessary\n");
        usage();
    }
}
















int main(int argc, char* argv[])
{
    char ppath[200], stringspath[200], occurrencespath[200], positionspath[200], bitmappath[200], allseqspath[200];
    

    get_options(argc, argv);


    //Get options and set initial parameters
    omp_set_num_threads(num_threads);       //number of threads

    sprintf(ppath, "%s/p", folder);
    sprintf(stringspath, "%s/%s.bin", folder, folder);
    sprintf(occurrencespath, "%s/occurrencesL%dA%d", folder, CONST_L,  ADDRESS_NUM_BITS);
    sprintf(positionspath, "%s/positionsL%dA%d", folder, CONST_L, ADDRESS_NUM_BITS);
    sprintf(bitmappath, "%s/bitmapL%dB%d", folder, CONST_L, BITMAP_NUM_BITS);
    sprintf(allseqspath, "%s/allseqsL%dA%d", folder, CONST_L, ADDRESS_NUM_BITS);


    //DEFINE VARIABLES
    mpz_t p_global;

    unsigned char *k_symbols = (unsigned char*) malloc(CONST_M_NUM_BYTES * sizeof(unsigned char));    
    position_t *positions = (position_t*) malloc((OCC_LEN + 1) * sizeof(position_t));
    sequence_t *precomp_seqs = (sequence_t*) malloc(CONST_A * sizeof(sequence_t));
    unsigned char* bitmap = (unsigned char*) calloc( BITMAP_NUM_BYTES, sizeof(unsigned char));



    // READ INPUT DATA (common.c)
    readprime(ppath, p_global);
    readsymbols(stringspath, k_symbols);



    // MAKE AND SORT SEQUENCES (common.c)
    makepositions(k_symbols, positions, p_global);
    make_allseqs_and_bitmap(precomp_seqs, k_symbols, bitmap, positions, p_global);
    sortallseqs(precomp_seqs, positions, num_threads);



    // WRITE DATA (common.c)
    writepositions(positionspath, positions);
    writebitmap(bitmappath, bitmap);
    writeallseqs(allseqspath, precomp_seqs);


    free(k_symbols);
    free(positions);
    free(precomp_seqs);
    free(bitmap);
    
    return 0;
}
