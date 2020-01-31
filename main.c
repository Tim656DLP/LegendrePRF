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

static uint64_t const_N;
static uint64_t const_k_N;
static int bits_N = 24;
static int seed = 0;
static int num_threads = 6;
static int read_data = 1;



static void usage()
{
    fprintf(stderr, "precomp -b bname [ options ]\n");
    fprintf(stderr, "  -b bname: project file name\n");
    fprintf(stderr, "  -n bits: number of bits of N, default: 24\n");
    fprintf(stderr, "  -s seed: randomness seed, default: time(NULL)\n");
    fprintf(stderr, "  -t num_threads: number of threads. default: 6\n");
    fprintf(stderr, "  -r read_data: read precomputed data (1) or create it (0). default: 1\n");

    fprintf(stderr, "  -h help\n");
    exit(1);
}




static void get_options(int argc, char **argv)
{
    char c;

    folder = NULL;

    while ((c = getopt(argc, argv, "b:n:s:t:r:h")) != (char)(-1))
    {
        switch (c)
        {
            case 'b':
                folder = optarg;
                break;
            case 'n':
                if (sscanf(optarg, "%d", &bits_N) != 1)
                    complain("Bad argument to -n!\n");
                break;
            case 's':
                if (sscanf(optarg, "%d", &seed) != 1)
                    complain("Bad argument to -s!\n");
                break;
            case 't':
                if (sscanf(optarg, "%d", &num_threads) != 1)
                    complain("Bad argument to -t!\n");
                break;
            case 'r':
                if (sscanf(optarg, "%d", &read_data) != 1)
                    complain("Bad argument to -r!\n");
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





static inline uint64_t generate_sequence(mpz_t R, mpz_t p, uint16_t seq_len)
{
    assert(seq_len <= 64);

    mpz_t tmp;
    mpz_init_set(tmp, R);

    uint64_t seq = 0;
    for (uint16_t j = 0; j < seq_len; j++)
    {
        unsigned char bit = (1 - mpz_legendre(tmp, p)) >> 1;
        seq |= ((uint64_t)bit << j);
        mpz_add_ui(tmp,tmp,1);
    }

    mpz_clear(tmp);
    return seq;
}





static inline uint64_t static_get_position(uint64_t i, position_t position)
{
    int64_t appx_pos = floor( (((double)CONST_A)/((double)(1ULL << ADDRESS_NUM_BITS))) * ((int64_t)i) );
    int64_t pos = (int64_t)position;

    while (1)
    {
        if(labs(appx_pos - pos) < (1ULL << (ADDRESS_NUM_BITS - 1)))
            return ((uint64_t)pos);
        pos += (1ULL << ADDRESS_NUM_BITS);
    }
}




static inline int binary_search_t(sequence_t *seqs, sequence_t val, uint16_t seqs_len)
{
    int32_t left, right, idx;

    idx = seqs_len >> 1;
    left = 0;
    right = seqs_len - 1;

    while (left <= right)
    {
        if (seqs[idx] < val)
            left = idx + 1;
        else if (seqs[idx] == val)
            return 0;

        else 
            right = idx - 1;

        idx = (left + right) >> 1;
    }

    return -1;
}



static inline int binary_search64(sequence_J *seq_J, uint64_t val, uint64_t seqs_len, uint32_t *dd, uint32_t *ii)
{
    int32_t left, right, idx;

    idx = seqs_len >> 1;
    left = 0;
    right = seqs_len - 1;

    while (left <= right)
    {
        if (seq_J[idx].seq < val)
            left = idx + 1;
        else if (seq_J[idx].seq == val)
        {
            *ii = seq_J[idx].i;
            *dd = seq_J[idx].d;
            return 0;
        }

        else 
            right = idx - 1;

        idx = (left + right) >> 1;
    }

    return -1;
}



static inline uint64_t static_get_sequence_from_id(unsigned char *k_symbols, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result)
{
    assert(seq_len <= 64);

    uint64_t seq = 0;
    for (int j = 0; j < seq_len; j++)
    {
        uint32_t idx = i + j*d;
        uint32_t byte_idx = idx >> 3;
        uint32_t bit_idx  = 7 - (idx & 0x7);

        unsigned char bit = (k_symbols[byte_idx] >> bit_idx) & 1;
        seq |= ((uint64_t)bit << j);
    }

    if (negate_result)
        seq ^= MASK(seq_len, 0);

    return seq;
}





static inline void static_get_two_sequences_from_id(unsigned char *J_list, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result, uint64_t *seq, uint64_t *seq_minus, unsigned char symbol_minus1)
{
    assert(seq_len <= 64);

    for (int j = 0; j < seq_len; j++)
    {
        uint64_t idx = i + j*d;
        uint64_t byte_idx = idx >> 3;
        uint64_t bit_idx  = 7 - (idx & 0x7);

        unsigned char bit = (J_list[byte_idx] >> bit_idx) & 1;
        *seq |= ((uint64_t)bit << j);
        *seq_minus |= ((uint64_t)bit << (seq_len - 1 - j));
    }

    if (negate_result)
        *seq ^= MASK(seq_len, 0);

    if (symbol_minus1^negate_result)
        *seq_minus ^= MASK(seq_len, 0);
}






static inline unsigned char bitmap_hit(unsigned char *bitmap, uint64_t seq)
{
    uint64_t bmp_pos = seq & MASK(BITMAP_NUM_BITS,0);
    uint64_t pos_byt = bmp_pos >> 3;
    unsigned char pos_bit = bmp_pos & 0x7;

    unsigned char hit = (bitmap[pos_byt] >> pos_bit) & 1;
    return hit;
}








static inline int is_the_key_correct(mpz_t key, mpz_t p, uint64_t *tests)
{
    mpz_t tmp;
    mpz_init_set(tmp, key);
    mpz_add_ui(tmp, tmp, (NUM_TESTS - 1) * TEST_BIT_CHUNKS);

    for (int i = NUM_TESTS - 1; i >= 0; i--)
    {
        uint64_t seq = generate_sequence(tmp, p, TEST_BIT_CHUNKS);
        if (seq != tests[i])
            return 0;
        mpz_sub_ui(tmp, tmp, TEST_BIT_CHUNKS);
    }

    mpz_clear(tmp);
    return 1;
}







int main(int argc, char* argv[])
{
    char ppath[200], stringspath[200], occurrencespath[200], positionspath[200], bitmappath[200], allseqspath[200], jspath[200];
    struct timespec start, end, parallel_end, parallel_start;

    get_options(argc, argv);


    //Get options and set initial parameters
    omp_set_num_threads(num_threads);       //number of threads
    const_N = (1ULL << bits_N);             //length of J_list
    const_k_N = (uint64_t)floor((double)(const_N - 1)/((double)(CONST_L - 1)));
    if(seed==0) seed=time(NULL);


    sprintf(ppath, "%s/p", folder);
    sprintf(stringspath, "%s/%s.bin", folder, folder);
    sprintf(positionspath, "%s/positionsL%dA%d", folder, CONST_L, ADDRESS_NUM_BITS);
    sprintf(bitmappath, "%s/bitmapL%dB%d", folder, CONST_L, BITMAP_NUM_BITS);
    sprintf(allseqspath, "%s/allseqsL%dA%d", folder, CONST_L, ADDRESS_NUM_BITS);


    //DEFINE VARIABLES
    mpz_t p_global;

    unsigned char *k_symbols = (unsigned char*) malloc(CONST_M_NUM_BYTES * sizeof(unsigned char));    
    position_t *positions = (position_t*) malloc((OCC_LEN + 1) * sizeof(position_t));
    unsigned char* bitmap = (unsigned char*) malloc(BITMAP_NUM_BYTES);
    sequence_t *precomp_seqs = (sequence_t*) malloc(CONST_A * sizeof(sequence_t));
    sequence_J *J_hits = (sequence_J*) calloc(MAX_J_HITS, sizeof(sequence_J)); 


    // READ INPUT DATA AND MAKE AUXILIARY DATA (common.c)
    readprime(ppath, p_global);
    readsymbols(stringspath, k_symbols);



    if(read_data)
    {
        readpositions(positionspath, positions);
        readbitmap(bitmappath, bitmap);
        readallseqs(allseqspath, precomp_seqs);
    }
    else
    {
        makepositions(k_symbols, positions, p_global);
        make_allseqs_and_bitmap(precomp_seqs, k_symbols, bitmap, positions, p_global);
        sortallseqs(precomp_seqs, positions, num_threads);       
    }








    // Start randomness
    srand(seed);
    gmp_randstate_t main_rstate;
    gmp_randinit_default(main_rstate);
    gmp_randseed_ui(main_rstate, rand());



    // WORK
    int stop = 0;
    int num_js = 0;

    uint64_t tot_trials = 0;
    uint64_t tot_hits = 0;


    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

while (!stop)
{
    printf("\nStarting the J sequence...\n");

    mpz_t rand_j; 
    mpz_init(rand_j);
    mpz_urandomm(rand_j, main_rstate, p_global);
    mpz_set_str(rand_j, "4533551310507856472906002", 10);
    gmp_printf("J = %Zd\n", rand_j);
    
    ++num_js;


    uint64_t num_trials = 0;
    uint64_t num_hits = 0;

    uint32_t curr_d = 0;





    // COMPUTE LIST OF SYMBOLS [J, J+N-1]
    unsigned char *J_list_g = (unsigned char*) malloc((const_N  + const_k_N ) >> 3);        
    generate_long_list(J_list_g, rand_j, p_global, const_N + const_k_N); 

    clock_gettime(CLOCK_MONOTONIC_RAW, &parallel_start);






#pragma omp parallel shared(curr_d, stop, num_hits) reduction(+: num_trials)
{
    mpz_t p, J, tmp;
    mpz_init_set(p, p_global);
    mpz_init_set(J, rand_j);
    mpz_init(tmp);
    uint32_t num_ds;
    uint32_t d, d_index;
    uint64_t seq_J, hit_index;


    // TIMING AND COUNTING
    int my_id = omp_get_thread_num();
    struct timespec thread_start, thread_end;
    uint64_t thread_trials = 0;


    // VARIABLES
    unsigned char symb_minus1;
    uint32_t *list_of_ds;


    // FUNCTIONS (common.h)
    symbolofmin1(p, &symb_minus1);
    denominators(const_N, &list_of_ds, &num_ds);


    // COMPUTE LIST OF SYMBOLS [J, J+N-1]
    unsigned char *J_list = (unsigned char*) malloc((const_N  + const_k_N) >> 3);        
    memcpy(J_list, J_list_g, (const_N + const_k_N)>>3); 


    int stoppriv = 0;
    clock_gettime(CLOCK_MONOTONIC_RAW, &thread_start);

    while (!stoppriv)
    {
        #pragma omp atomic capture
        d_index = curr_d++;

        if (d_index >= num_ds)
            break;
        else
            d = list_of_ds[d_index];


        // Compute the Legendre symbol of d
        mpz_set_ui(tmp, d);
        unsigned char d_symbol = (1 - mpz_legendre(tmp, p)) >> 1;


        // Generate the sequences and check for collisions
        for (uint32_t i = 0; i < d; i++)
        {
            #pragma omp atomic read
            stoppriv = stop;

            if (stoppriv)
                break;

            seq_J = static_get_sequence_from_id(J_list, i, d, CONST_L, d_symbol);

            if(bitmap_hit(bitmap, seq_J))
            {
                uint64_t addr = seq_J & MASK(ADDRESS_NUM_BITS, 0);
                sequence_t val = (sequence_t)((seq_J >> ADDRESS_NUM_BITS) & MASK(CONST_L - ADDRESS_NUM_BITS, 0));

                uint64_t pos_addrplus1 = static_get_position(addr + 1, positions[addr + 1]);
                uint64_t pos_addr = static_get_position(addr, positions[addr]);

                if(binary_search_t(&precomp_seqs[pos_addr], val, pos_addrplus1 - pos_addr) == 0)
                {
                    #pragma omp atomic capture
                    hit_index = num_hits++;
                    
                    if(hit_index >= MAX_J_HITS)
                        stop = 1;
                    else
                    {
                        J_hits[hit_index].seq = seq_J;
                        J_hits[hit_index].d = d;
                        J_hits[hit_index].i = i;
                    }
                }
            }

            thread_trials += 1;

        }

    }

    #pragma omp atomic update
    num_trials += thread_trials;


      
    clock_gettime(CLOCK_MONOTONIC_RAW, &thread_end);
    //printf("Thread %d finished %ld trials in %lf seconds\n", my_id, thread_trials, delta(&thread_end, &thread_start));

    mpz_clear(tmp);
    mpz_clear(J);
    mpz_clear(p);
    free(J_list);
    free(list_of_ds);
}

    clock_gettime(CLOCK_MONOTONIC_RAW, &parallel_end);

    if (num_hits == MAX_J_HITS) num_hits--;

    printf("J done!\n");
    printf("Number of trials = %" PRIu64 "\n", num_trials);
    printf("Number of hits is %" PRIu64 "\n", num_hits);
    printf("Finished the trials in %lf seconds\n", delta(&parallel_end, &parallel_start));
    printf("Now testing if hits are good...\n");

    tot_trials += num_trials;
    tot_hits += num_hits;



    gmp_sprintf(jspath, "%s/Jlist%Zd", folder, rand_j);

    // Sort J_hits
    qsort(J_hits, num_hits, sizeof(sequence_J), compareJ_seq);
    
    // Write rand_j and J_hits to a file
    FILE *fp = fopen(jspath, "wb");
    assert((fp != NULL) && "Error opening the J_hits output file");

    fwrite((sequence_J *)J_hits, num_hits * sizeof(sequence_J), 1, fp);
    fclose(fp);


    stop = 0;
    curr_d = 0;

    clock_gettime(CLOCK_MONOTONIC_RAW, &parallel_start);

if(num_hits)
#pragma omp parallel shared(curr_d, stop, k_symbols, num_hits)
{
    uint32_t ii, d, dd;
    mpz_t k, J, dinv, p, tmp;
    mpz_inits(k, dinv, tmp, NULL);
    mpz_init_set(J, rand_j);
    mpz_init_set(p, p_global);



    // TIMING AND COUNTING
    int my_id = omp_get_thread_num();
    struct timespec thread_start, thread_end;
    uint64_t thread_trials = 0;



    // VARIABLES
    unsigned char symbol_minus1;
    uint64_t tests[NUM_TESTS];


    // FUNCTIONS (common.h)
    symbolofmin1(p, &symbol_minus1);
    maketests(k_symbols, tests);


    int stoppriv = 0;
    clock_gettime(CLOCK_MONOTONIC_RAW, &thread_start);

    while(!stoppriv)
    {
        #pragma omp atomic capture
        d = curr_d++;

        if (d > CONST_k)
            break;

        mpz_set_ui(tmp, d);
        unsigned char d_symbol = (1 - mpz_legendre(tmp, p)) >> 1;

        for(uint32_t j = 0; j < d; j++)
        {
            #pragma omp atomic read
            stoppriv = stop;

            if (stoppriv)
                break;
     
            uint64_t seq = 0;
            uint64_t seq_minus = 0;

            for (uint32_t i = j; i < CONST_M-(CONST_L-1)*d; i+=d)
            {
                if (i == j)
                    static_get_two_sequences_from_id(k_symbols, i, d, CONST_L, d_symbol, &seq, &seq_minus, symbol_minus1);

                else
                {
                    uint32_t idx = i + (CONST_L - 1)*d;
                    uint32_t byte_idx = idx >> 3;
                    uint32_t bit_idx  = 7 - (idx & 0x7);

                    unsigned char bit = (k_symbols[byte_idx] >> bit_idx) & 1;

                    bit ^= d_symbol;
                    seq = (seq >> 1) | ((uint64_t)bit << (CONST_L - 1));
                    seq = seq & MASK(CONST_L, 0);

                    bit ^= symbol_minus1;
                    seq_minus = (seq_minus << 1) | ((uint64_t)bit);
                    seq_minus = seq_minus & MASK(CONST_L, 0);
                }



                if(binary_search64(J_hits, seq, num_hits, &dd, &ii) == 0)
                {
                    mpz_set_ui(dinv, dd);
                    mpz_invert(dinv, dinv, p);

                    mpz_set(k,J);
                    mpz_add_ui(k,k,ii);
                    mpz_mul(k,k,dinv);
                    mpz_mul_ui(k,k,d);
                    mpz_sub_ui(k,k,i);
                    mpz_mod(k,k,p);

                    if (is_the_key_correct(k, p, tests))
                    {
                        stop = 1;
                        gmp_printf("k = %Zd\n", k);
                        break;
                    }

                }

                if(binary_search64(J_hits, seq_minus, num_hits, &dd, &ii) == 0)
                {
                    mpz_set_ui(dinv, dd);
                    mpz_invert(dinv, dinv, p);

                    mpz_set(k,J);
                    mpz_add_ui(k,k,ii);
                    mpz_mul(k,k,dinv);
                    mpz_mul_ui(k,k,d);
                    mpz_mul_si(k,k,-1);
                    mpz_sub_ui(k,k,(CONST_L-1)*(uint64_t)d);
                    mpz_sub_ui(k,k,i);
                    mpz_mod(k,k,p);

                    if (is_the_key_correct(k, p, tests))
                    {
                        stop = 1;
                        gmp_printf("k = %Zd\n", k);
                        break;
                    }

                }


            }

        }

    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &thread_end);
    //printf("Thread %d finished testing hits in %lf seconds\n", my_id, delta(&thread_end, &thread_start));

    mpz_clears(tmp, J, p, k, dinv, NULL);
}

    clock_gettime(CLOCK_MONOTONIC_RAW, &parallel_end);
    printf("Finished testing hits in %lf seconds\n", delta(&parallel_end, &parallel_start));





} // END WHILE J LOOP

    
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    printf("\n\nTime taken to solve is %lf seconds\n", delta(&end, &start));
    printf("Total number of trials is %" PRIu64 "\n", tot_trials);
    printf("Total number of hits is %" PRIu64 "\n", tot_hits);
    printf("Total number of Js is %d\n", num_js);




    






    free(k_symbols);
    free(positions);
    free(precomp_seqs);
    free(bitmap);

    return 0;
}
