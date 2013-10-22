#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <peptide.h>
#include <binary_tree.h>   
#include <alignment.h>     


typedef enum {HOMOPOLYMER, JOINER, PRIMER, CHIMERAS, REFERENCE_FILTER} Program;
typedef enum {SANGER, ILLUMINA} Format;


/** The options we understand. */


/* Used by main to communicate with parse_opt. */
struct arguments
{
	char * alignments_file;
	char * alignments_fasta_file;
	char * input_file;
	char * output_file;
	char * peptide_file;
	char * reference_file;
	char * stats_file;
	
	Program program;
	Format alignments_input_format;
	
	int max_mismatches;
	int max_name_length;
	int max_read_length;
	int min_alignment_length;
	int min_count;
	int number_of_references;
	int quality_offset;
	int quality_threashold;
	int score_threshold;
	int silent, verbose;
	int primer_size;
	int max_length;
	
	boolean extreme;
	boolean peptide;
	char qual_window;
};

/* Parse a single option. */
void 
parse_opt (int key, char *arg,struct arguments *arguments )
{
	
	switch (key)
	{
		case 'a':
			arguments->alignments_file = arg;
			break;
		case 'c':
			arguments->min_count = atoi(arg);
			break;
			
		case 'f':
			if(strcmp(arg, "sanger") == 0 ){ 
				arguments->alignments_input_format = SANGER;
				arguments->quality_offset  = 33;
			}else if(strcmp(arg, "ilumina") == 0 ){ 
				arguments->alignments_input_format = ILLUMINA;
				arguments->quality_offset  = 64;
			}else {
				fprintf(stderr, "[-f sanger | ilumina] invalid input format (%s). \n", arg);
				exit(-1);
			}
			
			break;
			
		case 'i':
			arguments->input_file = arg;
			break;
		case 'l':	
			arguments->min_alignment_length = atoi(arg);
			break;
		case 'm':	
			arguments->max_mismatches = atoi(arg);
			break;
		case 'q':
			arguments->silent = 1;
			break;
			
		case 'o':
			arguments->output_file = arg;
			break;
			
		case 'p':
			arguments->peptide = true;
			arguments->peptide_file = arg;
			break;
		case 'r':
			//printf("%s\n",arg);
			arguments->reference_file = arg;
			break;	
		case 's':
			arguments->primer_size = atoi(arg);
			break;		
		case 'w':
			arguments->qual_window = atoi(arg);
		case 'v':
			arguments->verbose = 1;
			break;
		case 'A':
			arguments->alignments_fasta_file = arg;
			break;
		case 'C':
			arguments->program = CHIMERAS;
			break;
		case 'F':
			arguments->program = REFERENCE_FILTER;
			break;
		case 'H':
			arguments->program = HOMOPOLYMER;
			break;		
		case 'J':
			arguments->program = JOINER;
			break;
		case 'O':
			arguments->quality_offset  = atoi(arg);
			break;
		case 'P':
			arguments->program = PRIMER;
			break;
		case 'R':
			arguments->number_of_references = atoi(arg);
			break;
		case 'Q':
			arguments->quality_threashold  = atoi(arg);
			break;
		case 'S':
			arguments->stats_file = arg;
			break;
		case 'T':
			arguments->score_threshold = atoi(arg);
			break;
		case 'X':
			arguments->extreme = true;
			break;
			
		default:
			fprintf(stderr, "Invalid argument %c", key);
			abort();
			break;
	}
	
}

void align_parse(int argc, char **argv, struct arguments * arg){
	opterr = 0;
	int c;
	while ((c = getopt (argc, argv, "a:c:d:e:f:i:l:m:o:p:qr:s:vw:A:CFHJO:PQ:RXS:T:")) != -1)
	{
		
		parse_opt (c, optarg, arg);
	}
	
	
}

/*void joiner_get_next_sequences(Sequence * seq, struct arguments * args, FILE * first, FILE * second){
 
 switch (args->alignments_input_format) {
 case ILLUMINA:
 
 break;
 default:
 break;
 }
 }*/

void do_join(AlignmentToReference *al,struct arguments *args,FILE *align,Sequence *a,Sequence *b,FILE *peptide,FILE *stats,FILE *out) {
	if(args->extreme){
		alignment_to_reference_align_all(al, &extreme_nucleotide_score);
	}else{
		alignment_to_reference_align_all(al, &default_nucleotide_score);
	}
	
	
	if(align != NULL){
		alignment_to_reference_print(align, al);
	}
	
	/*fprintf(stdout, " %s\n",al->alignments[0]->b_aligned->seq);
	 fprintf(stdout, " %s\n",al->alignments[0]->b_aligned->qual);
	 fprintf(stdout, " %s\n",sequence_array_get_sequence(0, al->sequences)->seq);
	 fprintf(stdout, " %s\n",al->alignments[1]->b_aligned->seq);
	 fprintf(stdout, " %s\n",al->alignments[1]->b_aligned->qual);
	 fprintf(stdout, " %s\n",sequence_array_get_sequence(1, al->sequences)->seq);
	 */
	alignment_to_reference_consense_by_quality(al);
	
	sequence_set_name(a->name, al->concensus);
	sequence_append_name(" - ", al->concensus);
	sequence_append_name(b->name, al->concensus);
	
	alignment_to_reference_quality_filter_window(al);
	
	peptide_sequence_clean(al->concensus_protein);
	peptide_sequence_translate(al->concensus, 0, al->concensus_protein);	
	strcpy(al->concensus_protein->name, al->concensus->name);
	
	
	if(stats != NULL)
		alignment_to_reference_print_stats(stats, al);
	
	int correct ;
	correct = alignment_to_reference_count_correct_bases(al);
	
	if(correct > al->min_correct_bases){
		sequence_print_fasta(out, al->concensus);
		if(peptide != NULL){	
			peptide_sequence_print_fasta(peptide, al->concensus_protein);
		}
		
	}else if (args->verbose == 1){
		fprintf(stderr, "%s:Not enough bases(%d of %d) with good quality\n", al->concensus->name, correct, al->min_correct_bases );
	}
	alignment_to_reference_clean(al);
	
}
/*
 struct argp argp = { options, parse_opt, args_doc, doc };
 */




void  joiner( struct arguments  * args ){
	FILE * in;
	FILE * out;
	FILE * reference = NULL;
	FILE * peptide = NULL;
	FILE * stats = NULL;
	FILE * align = NULL;
	FILE * current_file = NULL;
	FILE * current_file2 = NULL;
	char filename[1000];
	
	int i = 0;
	int j = 0;
	
	if(strcmp(args->reference_file, "-") == 0){
		fprintf(stderr, "You have to specify a reference file with the option -r\n");
	}else{
		reference = fopen(args->reference_file, "r"); //open file of file names
		if(reference == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->reference_file);
			exit(-1);
		}
	}
	if(strcmp(args->input_file, "-") == 0){
		in = stdin;
	}else{
		in = fopen(args->input_file, "r"); //open file of file names
		if(in == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->input_file);
			exit(-1);
		}
	}
	
	//TODO Use it with more files, at the moment it only with one reference for two values...
	AlignmentToReference * al = alignment_to_reference_new(2, args->max_read_length, args->max_name_length);
	alignment_to_reference_set_threadshodlds(args->quality_threashold, args->score_threshold, al);
	al->concensus->qual_offset = args->quality_offset;//TODO, make a safer API that encapsulates and validates this offset through all the analysis
	al->qual_window = args->qual_window;
	//reference = fopen(args->reference_file, "r");
	boolean full_entry;
	if(read_sequence_from_fasta(reference, alignment_to_reference_get_reference(al), args->max_read_length, true,&full_entry ,0) <= 0){
		fprintf(stderr, "Unable to read reference in fasta format from the file %s\n", args->reference_file);
		exit(-1);
	}
	fclose(reference);
	
	al->min_correct_bases = (al->reference->length * args->min_alignment_length) / 100;
	
	if(args->peptide){//Do the translation to peptide if we want to do the analysis with peptides. 
		peptide_sequence_translate(alignment_to_reference_get_reference(al), 0, al->reference_protein);
		
		if(strcmp(args->output_file, "-") == 0){
			peptide = stdout;
		}else{
			peptide = fopen(args->peptide_file, "w"); //open file for peptides
		}
	}
	
	if(strcmp(args->output_file, "-") == 0){
		out = stdout;
	}else{
		out = fopen(args->output_file, "w"); //open file of alignments
	}
	
	Sequence * a =  alignment_to_reference_get_sequence(0, al);
	Sequence * b =   alignment_to_reference_get_sequence(1, al);
	a->qual_offset = args->quality_offset;
	b->qual_offset = args->quality_offset;
	//sequence_print_fasta(out, alignment_to_reference_get_reference(al));
	
	if(args->alignments_file != NULL){
		align = fopen(args->alignments_file, "w");
    }
    
    if(args->stats_file != NULL){
		if(strcmp(args->stats_file, "-") == 0){
			stats = stdout;
		}
		else{
			stats= fopen(args->stats_file, "w");
		}
		alignment_to_reference_print_stats_header(stats);
    }
    if (args->alignments_input_format == SANGER) {
		while (!feof(in)) {
			boolean both = false;
			i++;
			fscanf(in, "%s\n", filename);
			
			
			current_file = fopen(filename, "r");
			if(current_file==NULL){
				fprintf(stderr, "unable to open %s\n", filename);
				exit(-1);
			}
			if(read_sequence_from_fastq(current_file, a, args->max_read_length) > 0){
				//	sequence_remove_low_quality(a, args->quality_threashold);
				both = true;
			}else{
				both = false;
			}
			fscanf(in, "%s\n", filename);
			current_file = fopen(filename, "r");
			if(current_file==NULL){
				fprintf(stderr, "unable to open %s\n", filename);
				exit(-1);
			}
			if(read_sequence_from_fastq(current_file, b, args->max_read_length) > 0){
				//	sequence_remove_low_quality(b, args->quality_threashold);
				both &= true;
			}else{
				both &= false;
			}
			
			
			do_join(al,args,align,a,b,peptide,stats,out);
			
			//fprintf(out, "%s", filename);
			fclose(current_file);
		}
	}else if (args->alignments_input_format == ILLUMINA) {
		i = 0;
		j = 0;
		while (!feof(in)) {
			
			
			fscanf(in, "%s\n", filename);
			
			if (j%2 == 0) {
				current_file = fopen(filename, "r");
				
				if(current_file == NULL){
					fprintf(stderr, "unable to open %s\n", filename);
					exit(-1);
				}
			}else {
				current_file2 = fopen(filename, "r");
				
				if(current_file2 == NULL){
					fprintf(stderr, "unable to open %s\n", filename);
					exit(-1);
				}
			}
			
			
			
			do_join(al,args,align,a,b,peptide,stats,out);
			
			j++;
		}
	}
	
	if(align != NULL)
		fclose(align);
	if(peptide != NULL)
		fclose(peptide);
	fprintf(stderr, "%d Pairs aligned\n", i);
	fclose(in);
	fclose(out);
	
	alingment_to_reference_destroy(&al);
}

void  primer( struct arguments  * args ){
	FILE * in;
	FILE * out;
	FILE * reference = NULL;
	FILE * peptide = NULL;
	FILE * stats = NULL;
	FILE * align = NULL;
	
	int i = 0;
	if(strcmp(args->reference_file, "-") == 0){
		fprintf(stderr, "You have to specify a reference file with the option -r\n");
	}else{
		reference = fopen(args->reference_file, "r"); //open file of file names
		if(reference == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->reference_file);
			exit(-1);
		}
	}
	if(strcmp(args->input_file, "-") == 0){
		in = stdin;
	}else{
		in = fopen(args->input_file, "r"); //open file of file names
		if(in == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->input_file);
			exit(-1);
		}
	}
	
	//TODO Use it with more files, at the moment it only with one reference for two values...
	AlignmentToReference * al = alignment_to_reference_new(2, args->max_read_length, args->max_name_length);
	alignment_to_reference_set_threadshodlds(args->quality_threashold, args->score_threshold, al);
	al->concensus->qual_offset = args->quality_offset;//TODO, make a safer API that encapsulates and validates this offset through all the analysis
	al->qual_window = args->qual_window;
	//reference = fopen(args->reference_file, "r");
	boolean full_entry;
	
	Sequence * a =  alignment_to_reference_get_sequence(0, al);
	Sequence * b =   alignment_to_reference_get_sequence(1, al);
	if(read_sequence_from_fasta(reference, a, args->max_read_length, true,&full_entry ,0) <= 0){
		fprintf(stderr, "Unable to read first primer in fasta format from the file %s\n", args->reference_file);
		exit(-1);
	}
	if(read_sequence_from_fasta(reference, b, args->max_read_length, true,&full_entry ,0) <= 0){
		fprintf(stderr, "Unable to read second primer in fasta format from the file %s\n", args->reference_file);
		exit(-1);
	}
	fclose(reference);
	
	al->min_correct_bases = (al->reference->length * args->min_alignment_length) / 100;
	
	if(args->peptide){//Do the translation to peptide if we want to do the analysis with peptides. 
		peptide_sequence_translate(alignment_to_reference_get_reference(al), 0, al->reference_protein);
		peptide = fopen(args->peptide_file, "w");
	}
	
	if(strcmp(args->output_file, "-") == 0){
		out = stdout;
	}else{
		out = fopen(args->output_file, "w"); //open file of alignments
	}
	
	
	a->qual_offset = args->quality_offset;
	b->qual_offset = args->quality_offset;
	sequence_print_fasta(out, alignment_to_reference_get_reference(al));
	
	if(args->alignments_file != NULL){
		align = fopen(args->alignments_file, "w");
    }
    
    if(args->stats_file != NULL){
		stats = fopen(args->stats_file, "w");
		alignment_to_reference_print_stats_header(stats);
    }
	Sequence * ref =   alignment_to_reference_get_reference(al);
	while (read_sequence_from_fasta(in, ref, args->max_read_length,  true,&full_entry ,0) > 0) {
		
		i++;
		
		
		if(args->extreme){
			alignment_to_reference_align_all(al, &extreme_nucleotide_score);
		}else{
			alignment_to_reference_align_all(al, &default_nucleotide_score);
		}
		//
		//	alignment_to_reference_print(stdout, al);
		if(align != NULL)
			alignment_to_reference_print(align, al);
		
		alignment_to_reference_mask_primer(al);
		
		sequence_print_fasta(out, al->concensus);
		
		
		alignment_to_reference_clean_reference(al); 
		//fprintf(out, "%s", filename);
	}
	if(align != NULL)
		fclose(align);
	if(peptide != NULL)
		fclose(peptide);
	fprintf(stderr, "%d Pairs aligned\n", i);
	fclose(in);
	fclose(out);
	
	alingment_to_reference_destroy(&al);
}


void  homopolymer( struct arguments  * args ){
	FILE * in;
	FILE * out;
	FILE * align = NULL;
	FILE * reference = NULL;
	FILE * stats = NULL;
	
	Sequence * ref = NULL;
	MaskAlignment * ma = NULL;
	Sequence * tmp = sequence_new(args->max_read_length, args->max_name_length, 0);
	//BinaryTree * bt = binary_tree_new(200, &mask_alignment_compare);
    //BinaryTree * bt = binary_tree_new(200, &mask_alignment_simple_compare);
	boolean full_entry;
	if(strcmp(args->reference_file, "-") == 0){
		fprintf(stderr, "You have to specify a reference file with the option -r\n");
	}else{
		ref = sequence_new(args->max_read_length, args->max_name_length,0 );
		
		reference = fopen(args->reference_file, "r"); //open file of file names
		
		if(0 > read_sequence_from_fasta(reference, ref, args->max_read_length, true,&full_entry ,0) ){
			fprintf(stderr, "Unable to read reference in fasta format from the file %s\n", args->reference_file);
			exit(-1);
		}
		fclose(reference);
		
		
	}
	if(strcmp(args->input_file, "-") == 0){
		in = stdin;
	}else{
		in = fopen(args->input_file, "r"); //open file of file names
		if(in == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->input_file);
			exit(-1);
		}
	}
	
	if(strcmp(args->output_file, "-") == 0){
		out = stdout;
	}else{
		out = fopen(args->output_file, "w"); //open file of alignments
		if(out == NULL){
			fprintf(stderr, "Unable to open file for writing %s\n", args->output_file);
			exit(-1);
		}
	}
	
	
	int longest = 0, tmp_l = 0;
	int count = 0;
	fprintf(stderr, "min length %d\n", args->min_alignment_length);
	sequence_print_fasta(out, ref);
	while (read_sequence_from_fasta(in, tmp, args->max_read_length, true,&full_entry ,0) ) {
		count ++;
		if(sequence_get_length(tmp) > args->min_alignment_length){ 
			ma = mask_alignment_new(ref);
			ma->max_mismatches = args->max_mismatches;
			ma->search_window = args->qual_window;
			//		for (int i=0; i<args->primer_size; i++) {
			sequence_remove_base_up_to_limit(0, args->primer_size, tmp);
			//		}
			mask_alignment_set_sequence(tmp, ma);
			//mask_alignment_align(ma); This alignment is not required. The rest of the steps don't use it. It may be nice to use it to infer the starting point, but it is a unnecessary overhead. 
			//mask_alignment_fix_homopolymers(ma);
			//		
			
			mask_alignment_fix_insertions(ma);
			mask_alignment_fix_local_del_in(ma);
			mask_alignment_fix_local_in_del(ma);
			mask_alignment_fix_homopolymers(ma);
			mask_alignment_fix_deletions(ma);
			//mask_alignment_fix_homopolymers(ma);
			tmp_l = mask_alignment_longest_match(ma);
			
			if(tmp_l > longest){
				longest = tmp_l;
				fprintf(stderr, "Longest length: %d\n", longest);
			}
			sequence_set_name(tmp->name, ma->unmatched);
			if(sequence_get_length(ma->sequence)){
				sequence_print_fasta(out, ma->sequence);
				
			}
		}else{
			fprintf(stderr, "%d ignored\r", count);
			//	sequence_print_fasta(out, tmp);	
		}
		if(count %100 == 0){
			fprintf(stderr, "Processed %d sequences\r", count);
			fflush(stderr);
		}
		//	binary_tree_add_element(ma, bt);
		//	tmp = sequence_new(args->max_read_length, args->max_name_length);
	}
	fprintf(stderr, "Processed %d sequences\n", count);
	fclose(in);
	if(args->alignments_file != NULL){
		align = fopen(args->alignments_file, "w");
	}
	
	if(args->stats_file != NULL){
		stats = fopen(args->stats_file, "w");
		mask_alignment_print_stats_header(stats);
    }
	//   sequence_print_fasta(out, ref);
	void print_al(void * al){
		ma = (MaskAlignment * )al;
		
		if(align != NULL)
			mask_alignment_print_alignment_fasta(align, ma);		
		if(stats != NULL)
			mask_alignment_print_print_stats(stats, ma);
		sequence_print_fasta(out, ma->sequence);
		
		
	}
	
	//	binary_tree_sorted_walk(&print_al, bt);
	if(stats != NULL)
		fclose(stats);
	if(align !=  NULL)
		fclose(align);
	if(out != stdout)
		fclose(out);
	
	
}


int compare_count(void * aa, void * bb){
	Sequence * a = (Sequence * )aa;
	Sequence * b = (Sequence * )bb;
	if (a->count > b->count) {
		return -1;
	}else if (a->count < b->count) {
		return 1;
	}else {
		return strcmp(a->name, b->name);
	}
	
	
}

void  chimeras( struct arguments  * args ){
	FILE * in;
	FILE * out;
	FILE * align = NULL;
	FILE * reference = NULL;
	FILE * stats = NULL;
	
	Sequence * ref = NULL;
	
	Sequence * tmp = sequence_new(args->max_read_length, args->max_name_length,0);
	BinaryTree * seqs = binary_tree_new(200, &compare_count);
    //BinaryTree * bt = binary_tree_new(200, &mask_alignment_simple_compare);
	boolean full_entry;
	if(strcmp(args->reference_file, "-") == 0){
		fprintf(stderr, "You have to specify a reference file with the option -r\n");
	}else{
		ref = sequence_new(args->max_read_length, args->max_name_length, 0);
		
		reference = fopen(args->reference_file, "r"); //open file of file names
		
		if(0 > read_sequence_from_fasta(reference, ref, args->max_read_length, true,&full_entry ,0) ){
			fprintf(stderr, "Unable to read reference in fasta format from the file %s\n", args->reference_file);
			exit(-1);
		}
		fclose(reference);
		
		
	}
	if(strcmp(args->input_file, "-") == 0){
		in = stdin;
	}else{
		in = fopen(args->input_file, "r"); //open file of file names
		if(in == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->input_file);
			exit(-1);
		}
	}
	
	if(strcmp(args->output_file, "-") == 0){
		out = stdout;
	}else{
		out = fopen(args->output_file, "w"); //open file of alignments
		if(out == NULL){
			fprintf(stderr, "Unable to open file for writing %s\n", args->output_file);
			exit(-1);
		}
	}
	
	
	int seq_count;
	int count = 0;
	while (read_sequence_from_fasta(in, tmp, args->max_read_length, true,&full_entry ,0) ) {
		count++;
		seq_count = atoi(strrchr(tmp->name, '_')+1);
		//printf("%s: %d\n", tmp->name, seq_count);
		tmp->count = seq_count;
		Sequence * new_seq = sequence_new(args->max_read_length, args->max_name_length, 0);
		sequence_copy(new_seq, tmp);
		binary_tree_add_element(new_seq, seqs);
	}
	
	//RelatedSequences * rels[] = calloc(count, sizeof (RelatedSequences*));
	RelatedSequences* rels[count];
	
	count = 0;
	int i = 0;
	void do_alignments(void * v){
		Sequence * seq = (Sequence *) v;
		rels[count] = related_sequences_new(seq);
		
		void align_sequence(void * v2){
			Sequence * seq2 = (Sequence *) v2;
			if (v == v2) {
				return;
			}
			Alignment * al = alignment_new(seq->length+seq2->length, strlen(seq->name));
			//alignment_set_sequences(seq, seq2, al);
			//alignment_align_nucleotides(al);
			
			
			//	printf("aligning... %d\n", i++);
			
			alignment_set_sequences(seq, seq2, al);
			alignment_align_nucleotides(al);
			if (nump_score(al, 0.9, .1) != 0) {
				printf("Added_nump \n");
				sequence_print_fasta(stdout, seq);
				sequence_print_fasta(stdout, seq2);
				binary_tree_add_element(al, rels[count]->numps);
			}
			if (chimera_score(al, 0.95, 0.30, 0.80)!= 0) {
				printf("Added_chim \n");
				sequence_print_fasta(stdout, seq);
				sequence_print_fasta(stdout, seq2);
				binary_tree_add_element(al, rels[count]->chimeras);
			}
		}
		
		binary_tree_sorted_walk(&align_sequence, seqs);
		
		count++;
		
	}
	
	binary_tree_sorted_walk(&do_alignments, seqs);
	
	fclose(in);
	if(args->alignments_file != NULL){
		align = fopen(args->alignments_file, "w");
	}
	
	if(args->stats_file != NULL){
		stats = fopen(args->stats_file, "w");
		mask_alignment_print_stats_header(stats);
    }
    sequence_print_fasta(out, ref);
	
	for (i =0; i<count; i++) {
		related_sequences_to_fasta(out, rels[i]);
	}
	/*void print_al(void * al){
	 ma = (MaskAlignment * )al;
	 
	 if(align != NULL)
	 mask_alignment_print_alignment_fasta(align, ma);		
	 if(stats != NULL)
	 mask_alignment_print_print_stats(stats, ma);
	 sequence_print_fasta(out, ma->sequence);
	 
	 
	 }
	 
	 binary_tree_sorted_walk(&print_al, seqs);*/
	if(stats != NULL)
		fclose(stats);
	if(align!= NULL)
		fclose(align);
	if(out != stdout)
		fclose(out);
	
	
}

int compare_sequences(void * aa, void * bb){
	Sequence * a = (Sequence * )aa;
	Sequence * b = (Sequence * )bb;
	
	int compared = sequence_compare_with_ambiguity(a, b);
	
	if(compared == 0){
		sequence_merge_removing_ambiguity(a,b);
	}
	
	return compared;
}

void  filter_reference( struct arguments  * args ){
	FILE * in;
	FILE * out;
	FILE * reference = NULL;
	Sequence * ref = NULL;
	
	Sequence * tmp = sequence_new(args->max_read_length, args->max_name_length, 0);
	BinaryTree * seqs = binary_tree_new(100, &compare_sequences);
    //BinaryTree * bt = binary_tree_new(200, &mask_alignment_simple_compare);
	boolean full_entry;
	
	if(strcmp(args->reference_file, "-") == 0){
		fprintf(stderr, "You have to specify a reference file with the option -r\n");
	}else{
		ref = sequence_new(args->max_read_length, args->max_name_length, 0);
		
		reference = fopen(args->reference_file, "r"); //open file of file names
		
		if(0 > read_sequence_from_fasta(reference, ref, args->max_read_length, true,&full_entry ,0) ){
			fprintf(stderr, "Unable to read reference in fasta format from the file %s\n", args->reference_file);
			exit(-1);
		}
		fclose(reference);
		
		
	}
	
	if(strcmp(args->input_file, "-") == 0){
		in = stdin;
	}else{
		in = fopen(args->input_file, "r"); //open file of file names
		if(in == NULL){
			fprintf(stderr, "Unable to open file %s\n", args->input_file);
			exit(-1);
		}
	}
	
	if(strcmp(args->output_file, "-") == 0){
		out = stdout;
	}else{
		out = fopen(args->output_file, "w"); //open file of alignments
		if(out == NULL){
			fprintf(stderr, "Unable to open file for writing %s\n", args->output_file);
			exit(-1);
		}
	}
	
	int seq_count;
	int count = 0;
	
	int max_l = args->primer_size;
	int max_mismatches = args->max_mismatches;
	
	while (read_sequence_from_fasta(in, tmp, args->max_read_length, true,&full_entry ,0) ) {
		count++;
		seq_count = atoi(strrchr(tmp->name, '_')+1);
		if(seq_count == 0) seq_count = 1;
		if(DEBUG)
		fprintf(stderr,"%s: %d\n", tmp->name, seq_count);
			tmp->count = seq_count;
		
		sequence_trim(max_l, tmp);
		if(sequence_get_length(tmp)==max_l){
			
			Sequence * found_seq ;
			TreeElement * te=  binary_tree_find(tmp, seqs);
			
			if(te != NULL){
				found_seq = (Sequence *) te->element;
				found_seq->count += seq_count;
				//	fprintf(stderr, "FOUND!!!\n");
				//	fprintf(stderr, "ADDR %p\n",found_seq);
			}else {
				int mismatches = sequence_differences_with_mask(tmp, ref); 
				if(mismatches < max_mismatches){
					found_seq = sequence_new(args->max_read_length, args->max_name_length, 0);
					snprintf(tmp->name,tmp->max_name_length, "SEQ%d",  binary_tree_get_size(seqs) +1);
					sequence_copy(found_seq, tmp); 
					binary_tree_add_element(found_seq,seqs);
				}else{
					if(DEBUG)			fprintf(stderr, "mismatches: %d\n", mismatches);
				}
				
				//		fprintf(stderr, "ADDR %p\n",found_seq);
			}
			
			if(DEBUG){
				//
				//te=  binary_tree_find(tmp, seqs);
				if(found_seq!=NULL){
					//	found_seq =  (Sequence *) te->element;
				//	fprintf(stderr,"%s: %d\n", found_seq->name, found_seq->count);
				}
				
			}
			found_seq = NULL;
		}
	}
	int id = 0;
	void print(void * s){
		Sequence * seq = (Sequence *) s;
		if(seq->count >= args->min_count){
			//char ext[10];
			//			sprintf(&ext,"_%d", seq->count);
			//snprintf(name_ext,seq->max_name_length-seq->length, "_%d",  seq->count);
			//strcat(seq->name,ext);
			sprintf(seq->name, "Seq%d_%d",++id,seq->count);
			sequence_print_fasta(out, seq);
		}
	}
	
	binary_tree_sorted_walk(&print, seqs);
	
	if(out != stdout)
		fclose(out);
	
}

int main (int argc, char **argv)
{
	struct arguments arg;
	arg.silent = 0;
	arg.verbose = 0;
	arg.verbose = 0;
	arg.output_file = "-";
	arg.input_file = "-";
	arg.reference_file = "-";
	arg.quality_offset = 33;//Sanger fastq
	arg.max_read_length = 2000;
	arg.max_name_length = 500;
	
	arg.program = JOINER;
	
	//This are the default parameters for sanger. if selected Illumina,
	//the whole set will change, however 
	arg.alignments_input_format = SANGER;
	arg.quality_threashold = 25;
	arg.score_threshold = 200;
	arg.qual_window = 5;
	arg.min_alignment_length = 90;//If sanger, expressed in % of the reference sequence, if illumina in no. of bases. 
	//For homopolymers, is the minimum lenght of the sequence to be processed
	
	arg.stats_file = NULL;
	arg.alignments_file = NULL;
	arg.alignments_fasta_file = NULL;
	arg.max_mismatches = 1;
	arg.min_count = 2;
	
	arg.primer_size = 31;
	
	
	align_parse(argc, argv, &arg);
	
	//argp_parse (&argp, argc, argv, 0, 0, &arg);
	// print_args(&arg);
	//fprintf(stderr, "INPUT:%s\nOUPUT:%s\n", arg.input_file, arg.output_file);
	switch(arg.program){
		case JOINER:
			joiner(&arg);
			break;
		case HOMOPOLYMER:
			homopolymer(&arg);
			break;
		case PRIMER:
			primer(&arg);
			break;
		case CHIMERAS:
			chimeras(&arg);
			break;
		case REFERENCE_FILTER:
			filter_reference(&arg);
			break;
		default:
			fprintf(stderr, "Unknown program");
			exit(-1);
			break;
	}
	
	return 0;
}

