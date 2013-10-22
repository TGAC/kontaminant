#$: << File.expand_path(File.dirname(__FILE__) + '../../../../lib')
lib = File.expand_path(File.dirname(__FILE__) + '../../../../lib')
puts "Library: " + lib
$: << lib
require 'bio/cortex' 

kmer_size = 21
bucket_size = 100
max_rehash_tries = 10
buckets=16
fasta="/Users/ramirezr/Documents/TGAC/Workshops/ecoli_sample/ecoli/ecoli.fasta"
bad_reads= FFI::MemoryPointer.new :long_long
bad_reads.put_long_long(0,5)
puts "Initial bad reads: " + bad_reads.get_long_long(0).to_s + "(" + bad_reads.to_s + ")\n"

Bio::Cortex::analysis_init_cortex_environment(kmer_size)
#HashTable * hash_table_new(int number_bits, int bucket_size,int max_rehash_tries, short kmer_size);
graph = Bio::Cortex::hash_table_new(buckets, bucket_size, max_rehash_tries, kmer_size)
Bio::Cortex::db_graph_print_status(graph)
#    load_fasta_from_filename_into_graph(char *filename, short colour,long long *bad_reads, int max_chunk_length, dBGraph * db_graph)
#    load_fasta_from_filename_into_graph(char *filename, short colour,
#    				    long long *bad_reads,
#    				    int max_chunk_length, dBGraph * db_graph)
reads=Bio::Cortex::load_fasta_from_filename_into_graph(fasta, 0, bad_reads, 2000, graph)

puts "\nLoaded " + reads.to_s + " with " + bad_reads.get_long_long(0).to_s  + " bad reads\n"