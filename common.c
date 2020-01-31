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
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "gmp.h"
#include <time.h> 
#include <inttypes.h>
#include <string.h>
#include "common.h"




void complain(char *fmt, ...)
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    exit(1);
}




int compare8_seq(const void *a, const void *b)
{
    unsigned char a1 = *((sequence_t *)a);
    unsigned char b1 = *((sequence_t *)b);

    if (a1 < b1)
        return -1;
    else if (a1 > b1)
        return 1;
    else
        return 0;
}



int compare16_seq(const void *a, const void *b)
{
    uint16_t a1 = *((sequence_t *)a);
    uint16_t b1 = *((sequence_t *)b);

    if (a1 < b1)
        return -1;
    else if (a1 > b1)
        return 1;
    else
        return 0;
}



int compare32_seq(const void *a, const void *b)
{
    uint32_t a1 = *((sequence_t *)a);
    uint32_t b1 = *((sequence_t *)b);

    if (a1 < b1)
        return -1;
    else if (a1 > b1)
        return 1;
    else
        return 0;
}



int compareJ_seq(const void *a, const void *b)
{
    sequence_J a0 = *((sequence_J *)a);
    sequence_J b0 = *((sequence_J *)b);

    uint64_t a1 = a0.seq;
    uint64_t b1 = b0.seq;

    if (a1 < b1)
        return -1;
    else if (a1 > b1)
        return 1;
    else
        return 0;
}









uint64_t get_sequence_from_id(unsigned char *k_symbols, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result)
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










void get_two_sequences_from_id(unsigned char *J_list, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result, uint64_t *seq, uint64_t *seq_minus, unsigned char symbol_minus1)
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
        *seq = ~(*seq);

    if (symbol_minus1^negate_result)
        *seq_minus = ~(*seq_minus);
}













void generate_long_list(unsigned char *out_seq, mpz_t start, mpz_t p, uint64_t seq_len)
{
    //struct timespec end, time_start;
    //clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);

    mpz_t tmp;
    mpz_init_set(tmp, start);

    for (uint64_t i = 0; i < (seq_len >> 3); i++)
    {
        unsigned char res = 0;
        for (int j = 7; j >= 0; j--)
        {
            int sym = (1 - mpz_legendre(tmp, p)) >> 1;
            res |= ((unsigned char)sym << j);
            mpz_add_ui(tmp, tmp, 1);
        }
        out_seq[i] = res;
    }
    mpz_clear(tmp);


    //clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    //printf("Computed the J list in %lf seconds\n\n", delta(&end, &time_start));
}

























// AUXILIARY FUNCTIONS





void maketests(unsigned char* k_symbols, uint64_t *tests)
{
    for (uint32_t i = 0; i < NUM_TESTS; i++)
        tests[i] = get_sequence_from_id(k_symbols, i * TEST_BIT_CHUNKS, 1, TEST_BIT_CHUNKS, 0);  
}




void symbolofmin1(mpz_t p, unsigned char* symb)
{
    if (mpz_fdiv_ui(p, 4) == 1)
        *symb = 0;
    else
        *symb = 1;
}




void denominators(uint64_t N, uint32_t **list_of_ds, uint32_t *num_ds)
{
    // Make a list of prime d's in the range [(M-1)/(L-1) + 1, (N-1)/(L-1)]

    uint64_t k_N = (uint64_t)floor((double)(N - 1)/((double)(CONST_L - 1)));
    int Lbound = CONST_k + 1;
    int Ubound = k_N;
    int sieverange = Ubound - Lbound;
    uint32_t *all_ds_to_be_sieved = malloc( sieverange * sizeof(uint32_t));


    for (int j = 0; j < sieverange; j++)
        all_ds_to_be_sieved[j] = j + Lbound;

    

    for(int sieveprime = 2; sieveprime < sqrt(CONST_M); sieveprime++)
    {    
        int start = ((((-Lbound) % sieveprime) + sieveprime) % sieveprime);
     
       while (start < sieverange)
        {
            all_ds_to_be_sieved[start] = 0;
            start += sieveprime;
        }
    }


    uint32_t num_ds_global = 0;
    for (int j = 0; j < sieverange; j++)
        if (all_ds_to_be_sieved[j])
            num_ds_global++;


    int i = 0;
    *list_of_ds = malloc( num_ds_global * sizeof(uint32_t));
    
    for (int j = 0; j < Ubound - Lbound; j++)
        if (all_ds_to_be_sieved[j])
            (*list_of_ds)[i++] = all_ds_to_be_sieved[j];

    *num_ds = num_ds_global;
    free(all_ds_to_be_sieved);    
}




uint64_t get_position(uint64_t i, position_t position)
{
    int64_t appx_pos = floor( (((double)CONST_A)/((double)(1ULL << 32))) * ((int64_t)i) );
    int64_t pos = (int64_t)position;

    while (1)
    {
        if(labs(appx_pos - pos) < (1ULL << 31))
            return ((uint64_t)pos);
        pos += (1ULL << 32);
    }
}





void get_bitmap_index(uint64_t seq, uint64_t *byte, unsigned char *bit)
{
    uint64_t bmp_pos = seq & MASK(BITMAP_NUM_BITS,0);
    *byte = bmp_pos >> 3;
    *bit = bmp_pos & 0x7;
}

















// READ DATA

void readprime(char* primepath, mpz_t p)
{
    char pp[200];

    FILE *fp = fopen(primepath, "r");
    assert((fp != NULL) && "Error opening the prime number file");

    assert( fscanf(fp, "%s", pp) != 0) ;
    fclose(fp);

    mpz_init_set_str (p, pp, 16);
}





void readsymbols(char* stringspath, unsigned char* k_symbols)
{

    FILE *fp = fopen(stringspath, "rb");
    assert((fp != NULL) && "Error opening the symbols file");


    assert( fread(k_symbols, CONST_M_NUM_BYTES, 1, fp) != 0 );
    fclose(fp);

    // Flip the bits because they are computing Legendre as
    // (1+Leg(x,p))/2 instead of (1-Leg(x,p))/2
    // NOTE: remember to remove this when using data generated correctly
    for (uint64_t i = 0; i < CONST_M_NUM_BYTES; i++)
        k_symbols[i] = ~k_symbols[i];
}








void readpositions(char *positionspath, position_t *positions)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(positionspath, "rb");
    assert((fp != NULL) && "Error opening the positions file");

    assert( fread(positions, (OCC_LEN + 1) * sizeof(position_t), 1, fp) != 0);
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to read positions is %lf seconds\n", delta(&end,&start));
}





void readbitmap(char *bitmappath, unsigned char *bitmap)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(bitmappath, "rb");
    assert((fp != NULL) && "Error opening the input file with bitmap");
 
    assert(fread(bitmap, BITMAP_NUM_BYTES, 1, fp));
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to read bitmap is %lf seconds\n", delta(&end,&start));
}





void readallseqs(char *allseqspath, sequence_t *precomp_seqs)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(allseqspath, "rb");
    assert((fp != NULL) && "Error opening the input file with precomputed sequences");
 
    assert(fread(precomp_seqs, CONST_A * sizeof(sequence_t), 1, fp));
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to read precomputed sequences is %lf seconds\n", delta(&end,&start));
}


















// MAKE DATA

void makepositions(unsigned char *k_symbols, position_t *positions, mpz_t p_global)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    if (!positions)
        printf("positions malloc failed! Probably out of ram\n");

    uint32_t curr_d = 1;
    #pragma omp parallel shared(curr_d)
    {
        mpz_t tmp, p;
        mpz_init(tmp);
        mpz_init_set(p,p_global);
        uint32_t d;

        unsigned char symbol_minus1;
        symbolofmin1(p, &symbol_minus1);

        while(1)
        {
            #pragma omp atomic capture
            d = curr_d++;

            if (d > CONST_k)
                break;

            mpz_set_ui(tmp, d);
            unsigned char d_symbol = (1 - mpz_legendre(tmp, p)) >> 1;


            for(uint32_t j = 0; j < d; j++)
            {
                uint64_t seq = 0;
                uint64_t seq_minus = 0;

                for (uint32_t i = j; i < CONST_M-(CONST_L-1)*d; i+=d)
                {

                    if (i == j)
                        get_two_sequences_from_id(k_symbols, i, d, CONST_L, d_symbol, &seq, &seq_minus, symbol_minus1);
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

                    uint64_t addr1 = seq & MASK(ADDRESS_NUM_BITS, 0);
                    #pragma omp atomic
                    positions[addr1]++;

                    uint64_t addr2 = seq_minus & MASK(ADDRESS_NUM_BITS, 0);
                    #pragma omp atomic
                    positions[addr2]++;
                }

            }
        }
        mpz_clear(p);
        mpz_clear(tmp);
    }


    for (uint64_t i = OCC_LEN; i > 0 ; i--)
        positions[i] = positions[i-1];
    positions[0] = 0;

    for (uint64_t i = 1; i <= OCC_LEN ; i++)
        positions[i] += positions[i-1];


    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to make list of positions is %lf seconds\n", delta(&end, &start));
}











void make_allseqs_and_bitmap(sequence_t *precomp_seqs, unsigned char *k_symbols, unsigned char *bitmap, position_t *positions, mpz_t p_global)
{
    struct timespec start, end;

    if (!precomp_seqs)
        printf("Precomp_seqs malloc failed! Probably out of ram\n");

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    uint32_t curr_d = 1;   
    #pragma omp parallel shared(curr_d)
    {
        mpz_t tmp, p;
        mpz_init(tmp);
        mpz_init_set(p,p_global);
        uint32_t d;
        position_t pos0, pos1;
        uint64_t bmp_byte0, addr0, val0, index0;
        uint64_t bmp_byte1, addr1, val1, index1;
        unsigned char bmp_bit0;
        unsigned char bmp_bit1;
        unsigned char symbol_minus1;
        symbolofmin1(p, &symbol_minus1);

        while(1)
        {
            #pragma omp atomic capture
            d = curr_d++;

            if (d > CONST_k)
                break;

            mpz_set_ui(tmp, d);
            unsigned char d_symbol = (1 - mpz_legendre(tmp, p)) >> 1;


            for(uint32_t j = 0; j < d; j++)
            {
                uint64_t seq = 0;
                uint64_t seq_minus = 0;

                for (uint32_t i = j; i < CONST_M-(CONST_L-1)*d; i+=d)
                {
                    if (i == j)
                        get_two_sequences_from_id(k_symbols, i, d, CONST_L, d_symbol, &seq, &seq_minus, symbol_minus1);

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


                    //MAKE BITMAP
                    get_bitmap_index(seq, &bmp_byte0, &bmp_bit0);
                    #pragma omp atomic
                    bitmap[bmp_byte0] |= (unsigned char)(1 << bmp_bit0);

                    get_bitmap_index(seq_minus, &bmp_byte1, &bmp_bit1);
                    #pragma omp atomic
                    bitmap[bmp_byte1] |= (unsigned char)(1 << bmp_bit1);



                    // MAKE ALLSEQ
                    addr0 = seq & MASK(ADDRESS_NUM_BITS,0);
                    val0  = seq >> ADDRESS_NUM_BITS;
                    #pragma omp atomic capture
                    pos0 = positions[addr0]++;
                    index0 = get_position(addr0, pos0);

                    #pragma omp atomic write
                    precomp_seqs[index0] = ((sequence_t)val0);



                    addr1 = seq_minus & MASK(ADDRESS_NUM_BITS,0);
                    val1  = seq_minus >> ADDRESS_NUM_BITS;
                    #pragma omp atomic capture
                    pos1 = positions[addr1]++;
                    index1 = get_position(addr1, pos1);

                    #pragma omp atomic write
                    precomp_seqs[index1] = ((sequence_t)val1);


                }

            }

        }
        mpz_clear(tmp);
        mpz_clear(p);
    }



    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to make list of sequences and bitmaps is %lf seconds\n", delta(&end, &start));

    for(uint64_t i = OCC_LEN; i>0; i--)
        positions[i] = positions[i-1];
    positions[0] = 0;
}





void sortallseqs(sequence_t *precomp_seqs, position_t *positions, int num_threads)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    uint64_t number_of_threads = (uint64_t) num_threads;
    uint64_t my_id = omp_get_thread_num();

    #pragma omp parallel shared(num_threads)
    for (uint64_t i = my_id; i < OCC_LEN; i += number_of_threads)
    {
        uint64_t pos_iplus1 = get_position(i+1, positions[i+1]);
        uint64_t pos_i = get_position(i, positions[i]);

        if ((CONST_L-ADDRESS_NUM_BITS)==32)
            qsort(&(precomp_seqs[pos_i]), pos_iplus1 - pos_i, sizeof(sequence_t), compare32_seq);
        else if ((CONST_L-ADDRESS_NUM_BITS)==24)
            qsort(&(precomp_seqs[pos_i]), pos_iplus1 - pos_i, sizeof(sequence_t), compare32_seq);
        else if ((CONST_L-ADDRESS_NUM_BITS)==16)
            qsort(&(precomp_seqs[pos_i]), pos_iplus1 - pos_i, sizeof(sequence_t), compare16_seq);
        else if ((CONST_L-ADDRESS_NUM_BITS)==8)
            qsort(&(precomp_seqs[pos_i]), pos_iplus1 - pos_i, sizeof(sequence_t), compare8_seq);
        
    }   
    

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    printf("Time taken to sort list of sequences is %lf seconds\n\n", delta(&end, &start));
    printf("If all is correct %" PRIu32" = %" PRIu32" \n\n", (uint32_t)CONST_A, positions[OCC_LEN]);
}

















// WRITE DATA


void writepositions(char *positionspath, position_t *positions)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(positionspath, "wb");
    assert((fp != NULL) && "Error opening the positions output file");

    fwrite((position_t *)positions, (OCC_LEN + 1) * sizeof(position_t), 1, fp);
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to write list of positions is %lf seconds\n", delta(&end, &start));
}



void writebitmap(char *bitmappath, unsigned char *bitmap)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(bitmappath, "wb");
    assert((fp != NULL) && "Error opening the bitmap output file");

    fwrite((unsigned char *)bitmap, BITMAP_NUM_BYTES, 1, fp);
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to write bitmap is %lf seconds\n", delta(&end, &start));
}



void writeallseqs(char* allseqspath, sequence_t *precomp_seqs)
{
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    FILE *fp = fopen(allseqspath, "wb");
    assert((fp != NULL) && "Error opening the allseqs output file");

    fwrite((sequence_t *)precomp_seqs, CONST_A * sizeof(sequence_t), 1, fp);
    fclose(fp);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Time taken to write list of sequences is %lf seconds\n", delta(&end, &start));
}
















// TIMING

double delta(struct timespec *end, struct timespec *start)
{
   double delta_time = (double)(end->tv_sec - start->tv_sec) + (double)(end->tv_nsec - start->tv_nsec) / 1E9;
   return delta_time;
}

















/*







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










int is_the_key_correct(mpz_t key, mpz_t p, uint64_t *tests)
{
    mpz_t tmp;
    mpz_init_set(tmp, key);

    for (int i = 0; i < NUM_TESTS; i++)
    {
        uint64_t seq = generate_sequence(tmp, p, TEST_BIT_CHUNKS);
        if (seq != tests[i])
            return 0;
        mpz_add_ui(tmp, tmp, TEST_BIT_CHUNKS);
    }

    mpz_clear(tmp);
    return 1;
}






static inline int search_for_hit(sequence_t *precomp_seqs, uint32_t value, uint64_t bin_search_range, uint64_t *tests, mpz_t J, uint32_t i, uint32_t d, unsigned char is_minus, mpz_t p)
{
    uint64_t addr_id = 0;
    uint64_t index = 0;
    
    if (binary_search(precomp_seqs, value, bin_search_range, &addr_id, &index) == 0)
    {

        
        int move = 0;
        uint32_t ik, dk;
        mpz_t k, J_id, d_inv;
        mpz_init(k);
        mpz_init_set(J_id, J);
        mpz_init_set_ui(d_inv, d);
        mpz_invert(d_inv, d_inv, p);

        mpz_add_ui(J_id, J, i);
        mpz_mul(J_id, J_id, d_inv);
        mpz_mod(J_id, J_id, p);

        if(is_minus)
        {
            mpz_sub_ui(J_id, J_id, CONST_L - 1);
            mpz_mul_si(J_id, J_id, -1);
        }

    
        get_id_from_address(addr_id, &ik, &dk);

        mpz_set(k, J_id);
        mpz_mul_ui(k, k, dk);
        mpz_sub_ui(k, k, ik);
        mpz_mod(k, k, p);

        if (is_the_key_correct(k, p, tests)){
            gmp_printf("\nk = %Zd\n", k);
            return 2;
        }
    

        mpz_clear(J_id);
        mpz_clear(d_inv);
        mpz_clear(k);
        return 1;
    }
    
    return 0;

}





uint64_t get_address_from_id(uint32_t i, uint32_t d)
{   
    uint64_t address = i*CONST_k + (d-1);

    if (address >= CONST_A)
        address = 2*CONST_A - 1 - address;

    return address;
}







void get_id_from_address(uint64_t address, uint32_t *i, uint32_t *d)
{
    uint32_t ii;
    uint32_t dd;

    ii = address / CONST_k;
    dd = address % CONST_k;

    if (ii > CONST_M-1-(CONST_L-1)*(dd+1))
    {
        ii = 2*CONST_M - (CONST_L-1)*(CONST_k+1) - 1 - ii;
        dd = CONST_k - 1 - dd;
    }
    dd = dd + 1;

    *i = ii;
    *d = dd;
}

*/














