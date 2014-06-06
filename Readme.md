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

**Notes**

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

Filtering reads
===============

Finally, to filter the reads, use the kmer_filter program. 

To produce a range of fastq files binned by number of reads that have the contaminant use: 
```
./bin/kmer_filter_31 --reference reference.kmers --read_1 read_1.fq --output_prefix in_bins
```

Print a file with the "clean" reads. The flag --threashold is used as cutoff. Anything with more than --threashold kmers is considered a contaminant. If you also use --print_contaminated  you can get a fastq file with the reads that are contaminated (For example, you only want the reads from the same organism than your reference). 

```
./bin/kmer_filter_31 --reference referece.kmers --read_1 read_1.fq --read_2 read_2.fq  --output_prefix filtered --threashold 10
```


Flags:

	 --mem_width = Size of hash table buckets (default 100).
	 --help = This help message
	 --read_1 = File of reads in fastq format.
	 --read_2 = File of reads in fastq format. Optional. Should be in the same reads in the same order than the reads in read_1
	 --file_format = Format of the read. fasta or fastq. Default fastq
	 --kmer_size = Size of the kmer in the reference file
	 --mem_heigh = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).
	 --output_prefix = The prefix for the files with bins.
	 --quality_score_offset = Offset for the fastq file. Default 33.
	 --quality_score_threshold = Minimum quality across a kmer to be used
	 --reference = File containing the kmers to compare. This are built with kmer_hash_build
	 --print_contaminated = Print contaminated reads in a new fastq file
	 --threashold = Number of kmers to be considered as a above or below the cutoff
	 --output_folder = Folder to write the filtered reads.

**Notes**
You get a file with output_prefix.hist with the histogram of how many kmers have the reads to the contaminant. It is usally two big peaks with 0 and all the kmers. It can be useful to define the threashold. 

Screening for contamination
===========================
To know how many kmers you have in your sample, you need to

```
./bin/kmer_contamination_31 --reference_kmers contamiant_K21.kmers --input somple_k21 --output kmers.stat 
```
Flags:

	 --mem_width = Size of hash table buckets (default 100).
	 --help = This help message	 --input = File of reads in fastq format.
	 --file_format [KMERS] = Format of the read. fasta or fastq. Default KMERS
	 --kmer_size = Size of the kmer in the reference file
	 --mem_heigh = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).
	 --output = A tab separated file with the results of the analysis. The header is printed if the file doesn't exist. If it exists, the new stats are appeneded. We don't validate that previous entries are valid.
	 --quality_score_offset = Offset for the fastq file. Default 33. 
	 --quality_score_threshold = Minimum quality across a kmer to be used
	 --reference_kmers = File containing the kmers to compare. This are built with kmer_hash_build
	 --histogram_output = Output file where the histogram of encounter kmers per read is to be written. 

**Notes**

The columns in the output file are:
 *  **REFERENCE** Filename of the reference
 *	**SAMPLE** Filename of the sample
 *	**REFERENCE_READS** Number of reads in the reference
 *  **SAMPLE_READS** Number of reads in the sample
 *	**KMER_SIZE** Kmer size. MUST be consistent between sample and reference
 *	**REFERENCE_KMERS** Total count of kmers in the reference
 *	**REFERENCE_COVERAGE** Average coverage of the reference
 *  **SAMPLE_KMERS** Total count of kmers in the sample
 *	**SAMPLE_COVERAGE** Average coverage of the sample
 *	**PERCENTAGE** Percemtage of kmers sample in reference. Useful to know how much of the contaminant is precent. If it is low, (less than 5%, depending on how similar are the genomes), the overlapping kmers may represent common regions. 
 *	**COMMON** count of common kmers
 *	**SAMPLE_READS** Reads in the sample
 *	**CONTAMINTED_READS** Number of contaminated reads (see BUGS)
 *	**PERCENTAGE_CONTAMINATED_READS**  Number of contaminated reads (see BUGS)

**BUG:** The FASTA and FASTQ file reader is not loading the kmers in the hash, hence only the KMERS file format is supported
**BUG:** The histogram_output depends on being able to read from fastq and fasta. 
