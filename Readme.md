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
The drawback of having a bigger kmer size is that Kontaminant will take more memory to run. 

If the compilation was sucesful, you wil find in your bin folder 3 programs: 
 * *kmer_hash_build_MAXK*  Builds the kmer reference
 * *kmer_contamination_MAXK* Detects if the contaminant is present or not (and in what level) in the sample
 * *kmer_filter_MAXK* Filters the reads with the contaminant out of the reference. 

Building the hash
-----------------
To build the database of kmers from a fasta file:

```
 ./bin/kmer_hash_build_31 --input contaminant_reference.fasta --output contamiant_K21.kmers --format FASTA 
```

*Notes*:

If you get the error:
```
too much rehashing!! Rehash=26
```
try increasing the --mem_height

The output file contains the kmers in the contaminant and it is the database to be used in the rest of the programs. 

Flags:

	 --mem_width = Size of hash table buckets (default 100).
	 --help = This help message
	 --input = File of reads in fastq or fasta format. If not privided, STDIN is used
	 --file_format [FASTQ | FASTA] = Format of the reads. fasta or fastq. Default fastq 
	 --kmer_size = Size of the kmer in the reference file
	 --mem_height = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).
	 --output = A kmers file in binary format. To be used with the kmer_stats and kmer_contamination programs.
	 --quality_score_offset = Offset for the fastq file. Default 33.
	 --quality_score_threshold = Minimum quality across a kmer to be used



