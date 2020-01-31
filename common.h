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


//#define CONST_L (64)          //GIVEN IN MAKE
//#define ADDRESS_NUM_BITS (32) //GIVEN IN MAKE
//#define BITMAP_NUM_BITS (38) //GIVEN IN MAKE


#define CONST_M  (1 << 20) // Number of oracle calls
#define CONST_k  ((uint64_t)(floor((double)(CONST_M-1)/(double)(CONST_L-1))))
#define CONST_A  (2*(CONST_M*CONST_k - (((CONST_L-1)*CONST_k*(CONST_k+1))>>1)) )
#define CONST_M_NUM_BYTES  (CONST_M >> 3) 

#define OCC_LEN (1ULL << ADDRESS_NUM_BITS)
#define BITMAP_NUM_BYTES (1ULL << (BITMAP_NUM_BITS - 3))

#define MASK(end, start) (((-(1ULL)) >> (64 - (end - start))) << start) // Compute mask from start bit to end-1 bit

#define MAX_J_HITS (1ULL << 15)

#define NUM_TESTS  (3)
#define TEST_BIT_CHUNKS (64)




#if (ADDRESS_NUM_BITS==32)
typedef uint32_t position_t;
#elif (ADDRESS_NUM_BITS==24)
typedef uint32_t position_t;
#elif (ADDRESS_NUM_BITS==16)
typedef uint16_t position_t;
#elif (ADDRESS_NUM_BITS==8)
typedef unsigned char position_t;
#endif


#if((CONST_L-ADDRESS_NUM_BITS)==32)
typedef uint32_t sequence_t;
#elif((CONST_L-ADDRESS_NUM_BITS)==24)
typedef uint32_t sequence_t;
#elif((CONST_L-ADDRESS_NUM_BITS)==16)
typedef uint16_t sequence_t;
#elif((CONST_L-ADDRESS_NUM_BITS)==8)
typedef unsigned char sequence_t;
#endif

typedef struct{
    uint64_t seq;
    uint32_t d;
    uint32_t i;
} sequence_J;




void complain(char *fmt, ...);

int compare8_seq(const void *a, const void *b);
int compare16_seq(const void *a, const void *b);
int compare32_seq(const void *a, const void *b);
int compareJ_seq(const void *a, const void *b);


//  AUXILIARY
void maketests(unsigned char* k_symbols, uint64_t *tests);
void symbolofmin1(mpz_t p, unsigned char* symb);
void denominators(uint64_t N, uint32_t **list_of_ds_global, uint32_t *num_ds);
void get_bitmap_index(uint64_t seq, uint64_t *byte, unsigned char *bit);
uint64_t get_position(uint64_t i, position_t position);



uint64_t get_sequence_from_id(unsigned char *k_symbols, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result);
void get_two_sequences_from_id(unsigned char *J_list, uint32_t i, uint32_t d, uint16_t seq_len, unsigned char negate_result, uint64_t *seq, uint64_t *seq_minus, unsigned char symbol_minus1);
void generate_long_list(unsigned char *out_seq, mpz_t start, mpz_t p, uint64_t seq_len);




// INPUT DATA
void readprime(char* primepath, mpz_t p);
void readsymbols(char* stringspath, unsigned char* k_symbols);




// READING
void readpositions(char *positionspath, position_t *positions);
void readbitmap(char *bitmappath, unsigned char *bitmap);
void readallseqs(char *allseqspath, sequence_t *precomp_seqs);


// MAKING
void makepositions(unsigned char *k_symbols, position_t *positions, mpz_t p_global);
void make_allseqs_and_bitmap(sequence_t *precomp_seqs, unsigned char *k_symbols, unsigned char *bitmap, position_t *positions, mpz_t p_global);
void sortallseqs(sequence_t *precomp_seqs, position_t *positions, int num_threads);


// WRITING
void writepositions(char *positionspath, position_t *positions);
void writebitmap(char *bitmappath, unsigned char *bitmap);
void writeallseqs(char* allseqspath, sequence_t *precomp_seqs);


// TIMING
double delta(struct timespec *end, struct timespec *start);




// OLD
//static inline uint64_t generate_sequence(mpz_t R, mpz_t p, uint16_t seq_len);
//static inline int is_the_key_correct(mpz_t key, mpz_t p, uint64_t *tests);
//static inline int search_for_hit(sequence_t *precomp_seqs, uint32_t value, uint64_t bin_search_range, uint64_t *tests, mpz_t J, uint32_t i, uint32_t d, unsigned char is_minus, mpz_t p);
//static inline uint64_t get_address_from_id(uint32_t i, uint32_t d);
//static inline void get_id_from_address(uint64_t address, uint32_t *i, uint32_t *d);
//int binary_search(sequence_t *seqs, uint32_t val, uint16_t seqs_len, uint64_t *address, uint64_t *index);
