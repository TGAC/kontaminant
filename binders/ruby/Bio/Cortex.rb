require 'rubygems'
require'ffi'


module Bio
  module Cortex
    extend FFI::Library
    #ffi_lib 'libctx'
    puts "Library: "+  File.dirname(__FILE__) + '/../../../lib/libctx_31'
    ffi_lib File.dirname(__FILE__) + '/../../../lib/libctx_31.dylib'
    attach_function :analysis_init_cortex_environment, [ :short ], :void
#    HashTable * hash_table_new(int number_bits, int bucket_size,int max_rehash_tries, short kmer_size);
    attach_function :hash_table_new, [:int, :int, :int, :short], :pointer
    attach_function :db_graph_print_status, [:pointer], :void
#       long long
#    load_fasta_from_filename_into_graph(char *filename, short colour,
#    				    long long *bad_reads,
#    				    int max_chunk_length, dBGraph * db_graph)
    attach_function :load_fasta_from_filename_into_graph, [:string, :short, :pointer, :int, :pointer], :long_long
  end
end