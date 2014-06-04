Kontaminant
===========
Kontaminant is a kmer-based screening and filtering tool for next-gen sequencing reads.
It shares its main data structures with [cortex_con](http://cortexassembler.sourceforge.net)

Compiling 
---------
Kontaminant is developed and tested on gcc. Other C compilers are unlikely to work due the use of nested-functions in the use of data structures. To complile just do:

```sh
$ make
```
The default maximum k-mer size is 31. You can specify the following maximum kmer sizes: 31, 63, 95, 127, 160, 192, 223 and 255. with the following flag:

```sh
$ make MAXK=63
```
