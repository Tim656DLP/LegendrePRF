# LegendrePRF

LegendrePRF is an algorithm for key recovery of the Legendre pseudorandom function written in C.
You can find more information on the algorithm in https://eprint.iacr.org/2020/098 .

## Prerequisites

You will need [gmp](https://gmplib.org) and [openmp](https://www.openmp.org).

## Installation


Call the makefile by using

```bash
make L=64 A=32 B=37
```
Two executables will be created, precomp.out and leg.out.
You will need at least 97GB or RAM.
## Usage

To run the precomputation phase call "precomp.out". It will generate the hash table together with the auxiliary data and write it in the project folder.

      -b bname: project file name
      -t num_threads: number of threads. default: 1
      -h help

```bash
./precomp.out -b 0 -t 12
```

To run the search phase call "leg.out". It will read the precomputed data and start looking for collisions.
Found collisions will be written in the project folder.
The program halts after finding the key k and printing it out.

    -b bname: project file name
    -n bits: number of bits of N, default: 24
    -s seed: randomness seed, default: time(NULL)
    -t num_threads: number of threads. default: 6
    -r read_data: read precomputed data (1) or create it (0). default: 1
    -h help

```bash
./leg.out -b 0 -n 23 -t 12
```


If do not wish to write the precomputed data to a disk run only the "leg.out" script with -r 0.
```bash
./leg.out -b 0 -n 23 -t 12 -r 0
```

## Challenge 2
For Challenge 2 you may run the code with -n 26 and -s 122. The 204th value of j=4533551310507856472906002 will output a correct solution.

```bash
./leg.out -b 2 -n 26 -s 122 -t 12 -r 0
```
Challenge 2 was solved on the evening between 27th and 28th of November after less than two days of computation - https://tinyurl.com/84bitprf . The resulting key is 187320452088744099523844.
## License
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
