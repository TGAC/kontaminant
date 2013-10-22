/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team: 
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <err.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include <global.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <seq.h>
#include <limits.h>
#include <logger.h>

static char BASES[10] = { 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n' };

/* Setup expected quality values */
void sequence_set_quality_parameters(Sequence *seq, int quality_offset)
{
    seq->qual_offset = quality_offset;
    if (quality_offset == 33) {
        seq->check_quality_values = true;
        seq->lowest_expected_value = 0;
        seq->highest_expected_value = 93;
        seq->highest_raw_read_value = 60;
    } else if (quality_offset == 64) {
        seq->check_quality_values = true;
        seq->lowest_expected_value = 0;
        seq->highest_expected_value = 62;
        seq->highest_raw_read_value = 40;
    } else {
        seq->check_quality_values = false;
        log_and_screen_printf("Warning: unexpected quality offset passed to sequence_set_quality_parameters.\n");
    }
}

/*
 * Read sequence from file "fp" in FASTA format. it returns the length of
 * the sequence, 0 if no sequence is left in file 
 */
boolean good_base(char c);
boolean nucleotide_good_base(char c);

char current_entry_name[LINE_MAX + 1] = "";
int last_end_coord;

// new_entry tells the parser to expect a new fasta entry (ie starts with
// >)
// full_entry tells the caller if the end of the entry has been reached
// This function can be use with big fasta entries. Successive calls on
// the same
// big entry remember the name of the entry in seq->name.
// seq->start, seq->end are the coordinate with respect to the full fasta
// entry (eg chr1 start: 1000 end: 3000)
// offset defines the length of seq->seq that is preserved, this allows to 
// append new sequence from the entry to a prefix in seq->seq
// until max_chunk_length is reached. 
// ir returns 0 when reaches EOF (in this case full_entry is true)

int read_sequence_from_fasta(FILE * fp, Sequence * seq, int max_chunk_length,
                             boolean new_entry, boolean * full_entry, int offset)
{
    
	char line[LINE_MAX];	// LINE_MAX is defined in limits.h
	int i;
	int j = 0;		// length of sequence
	boolean good_read;
    
	if (fp == NULL) {
		fputs("File not defined\n", stderr);
		exit(1);
	}
    
	if (seq == NULL) {
		fputs("Cannot pass a NULL pointer for seq\n", stderr);
		exit(1);
	}
    
	if (seq->seq == NULL) {
		fputs
        ("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->qual == NULL) {
		fputs
        ("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->name == NULL) {
		fputs
        ("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
	// get name
    boolean name_set = false;
	if (new_entry == true) {	// need a '>' followed by a name
		while (name_set == false && fgets(line, LINE_MAX, fp) != NULL) {
            
            if(line[0] == '#'){
                
                //A comment, we ignore it. 
            }else if (line[0] == '>') {
				i = 1;
				while (i < LINE_MAX
				       && !(line[i] == '\n' || line[i] == ' '
                            || line[i] == '\t'
                            || line[i] == '\r')) {
                           seq->name[i - 1] = line[i];
                           i++;
                       }
				seq->name[i - 1] = '\0';
                name_set = true; //RHRG: As we know the name, we can continue. 
			} else {
				fprintf(stderr,
                        "syntax error in fasta entry %s\n",
                        line);
				assert(0);
                exit(1);
			}
		}
        
		strcpy(current_entry_name, seq->name);
		seq->start = 1;	// new entry
	} else {		// old entry
		strcpy(seq->name, current_entry_name);
		seq->start = last_end_coord + 1;
	}
    
	// get sequence
    
	boolean ready_to_return = false;
	long file_pointer = 0;
	long prev_file_pointer = 0;
    
	j = offset;		// global length
	*full_entry = false;
    
	file_pointer = ftell(fp);
    
	while (ready_to_return == false) {
        
		file_pointer = ftell(fp);
        
		if (fgets(line, LINE_MAX, fp) != NULL) {
            
			size_t length = strlen(line);
            
			// sanity check
            if (line[0] == '>') {
				fseek(fp, file_pointer, SEEK_SET);
				ready_to_return = true;
				*full_entry = true;
			} else {
				if (j == max_chunk_length) {
					fseek(fp, prev_file_pointer, SEEK_SET);
					ready_to_return = true;
					*full_entry = false;
				} else {
					// check line is fine
					i = 0;	// counter within line
					while (i < length
					       && !(line[i] == '\n'
                                || line[i] == ' '
                                || line[i] == '\t'
                                || line[i] == '\r')
					       && ready_to_return == false) {
						good_read = good_base(line[i]);
						if (!good_read) {
							fprintf(stderr,
                                    "Invalid symbol [%c] pos:%i in entry %s\n",
                                    line[i], i,
                                    seq->name);
                            
						}
                        
						seq->seq[j] = line[i];
						seq->qual[j] = '\0';
						i++;
						file_pointer++;
						j++;
                        
						if (j == max_chunk_length) {
							// check if line is not complete
							if (i < length - 1) {
								ready_to_return = true;
								*full_entry =false;
								fseek(fp, file_pointer, SEEK_SET);
							}
							// else 3 cases might happen with line
							// complete
							// 1. next line starts ">" -> *full_entry =
							// true
							// 2. next line is EOF -> *full_entry = true
							// 3. next line is more sequence ->
							// *full_entry = false
						}
					}
				}
			}
		} else {
			*full_entry = true;
			ready_to_return = true;
			// printf("complete line - complete entry - complete file\n");
		}
        
		prev_file_pointer = file_pointer;
	}
    
	seq->end = seq->start + j - 1 - offset;
	last_end_coord = seq->end;
    
	seq->seq[j] = '\0';
	seq->qual[j] = '\0';
	seq->length = j;
    
	return j;
}

/*
 * Read sequence from file "fp" in FASTQ format. it returns the length of
 * the sequence, 0 if no sequence is left in file read with bad chars are
 * skipped -- this is a different behaviour to read_sequence_from_fasta  
 *  
 */
static int read_sequence_from_fastq_from_stream(FILE * fp, Sequence * seq, int max_read_length)
{
    // printf("READING FROM STDIN NOT IMPLEMENTED\n");
    
	char line[LINE_MAX];	// LINE_MAX is defined in limits.h
	int i;
	int j = 0;		// length of sequence
	int q = 0;		// length of qualities
    
    int offset = seq->qual_offset;    
	if (offset == 0) {
		offset = 64;	// Ilumina fastq format, default to keep
		// backwards compatibility. 
	}  
    
    assert(fp != NULL);
    assert(seq != NULL);  
	assert(seq->seq != NULL);
    assert(seq->name != NULL);
    assert(seq->qual != NULL); 
    sequence_clean(seq);
    char * readed;
    //Read the name
    readed = fgets(line, LINE_MAX, fp);
    if (readed == NULL) {
        return 0;
    }
    
    if (line[0] != '@') {
        fputs("syntax error in fastq file -- it misses @\n", stderr);
        assert(0);
    }
    for (i = 1; i < LINE_MAX; i++) {
        if (line[i] == '\n'// || line[i] == ' ' We alow the space to be able to process Illumina 1.8+
            || line[i] == '\t'
            || line[i] == '\r') {
            break;
        }
        if (i > LINE_MAX) {
            fputs("Name too long\n",
                  stderr);
            exit(1);
        }
        seq->name[i - 1] = line[i];
    }
    seq->name[i - 1] = '\0';
    
    //Read the sequence
    readed = fgets(line, LINE_MAX, fp);
    if (readed == NULL) {
        sequence_clean(seq);
        return 0;
    }
    
    j=0;
    for (i = 0; i < LINE_MAX; i++) {
        if (line[i] == '\n'
            || line[i] == ' '
            || line[i] == '\t'
            || line[i] == '\r') {
            break;	// fine but nothing to add
        }
        
        if(line[i] == '.'){//In case you have dots instad of Ns, that come from converting a qseq to fastq in a "quick and dirty" way
            line[i] = 'N';
        }
        
        if (!good_base(line[i])) {
            //            good_read = false;
            fprintf(stderr,
                    "Invalid symbol [%c] pos:%i in entry %s - skip read\n",
                    line[i], i,
                    seq->name);
        }
        
        seq->seq[j] = line[i];
        j++;
        
        if (j == max_read_length) {
            fprintf(stdout,
                    "read [%s] too long [%i]. Exiting...\n",
                    seq->name, j);
            assert(NULL);
            exit(1);
        }
        
    }
    
    //Line with the name for quality
    readed = fgets(line, LINE_MAX, fp);
    if (readed == NULL) {
        sequence_clean(seq);
        return 0;
    }    
    // if (line[0] != '+' || line[0] != '-') {
    //     fputs("syntax error in fastq file -- missing + or - \n", stderr);
    //    assert(0);
    //}
    
    
    //Line with the qualities. 
    readed = fgets(line, LINE_MAX, fp);
    if (readed == NULL) {
        sequence_clean(seq);
        return 0;
    }
    q=0;
    for (i = 0; i < LINE_MAX; i++) {
        if (line[i] == '\n'
            || line[i] == ' '
            || line[i] == '\t'
            || line[i] == '\r') {
            break;	// fine but nothing to add
        }
        
        seq->qual[q] = line[i] - offset;
        
        if (seq->check_quality_values) {
            if (seq->qual[q] < seq->lowest_expected_value) {
                log_and_screen_printf("Warning: Quality [%d] for [%s] lower than expected for specified quality offset [%d-%d].\n", seq->qual[q], seq->name, seq->lowest_expected_value, seq->highest_expected_value);
                seq->check_quality_values = false;
            } else if (seq->qual[q] > seq->highest_expected_value) {
                log_and_screen_printf("Warning: Quality [%d] for [%s] higher than expected for specified quality offset [%d-%d].\n", seq->qual[q], seq->name, seq->lowest_expected_value, seq->highest_expected_value);
                seq->check_quality_values = false;
            } else if (seq->qual[q] > seq->highest_raw_read_value) {
                log_and_screen_printf("Warning: Quality [%d] for [%s] higher than expected for raw reads [%d].\n", seq->qual[q], seq->name, seq->highest_raw_read_value);
                seq->check_quality_values = false;
            }         
        }
        
        q++;
        
        if (q == max_read_length) {
            fprintf(stdout,
                    "qualities for [%s] longer than the max read length  [%i]. Exiting...\n",
                    seq->name, q);
            exit(1);
        }
        
    }
    
    if (j != q) {
        fprintf(stdout,
                "Lengths of quality [%i] and sequence [%i]  strings don't coincide for this read: [%s]. Skip it\n",
                q, j, seq->name);
        sequence_clean(seq);
        //       good_read = false;
    }else{
        seq->length = j;
    }
    return seq->length;
}

/*
 * Read sequence from file "fp" in FASTQ format. it returns the length of
 * the sequence, 0 if no sequence is left in file read with bad chars are
 * skipped -- this is a different behaviour to read_sequence_from_fasta 
 * This only works with files, as it relies on moving the pointer. 
 */

static int read_sequence_from_fastq_from_file(FILE * fp, Sequence * seq, int max_read_length)
{
    
	char line[LINE_MAX];	// LINE_MAX is defined in limits.h
	int i;
	int j = 0;		// length of sequence
	int q = 0;		// length of qualities
	long file_pointer;
	boolean good_read = true;
	//int offset = fastq_ascii_offset; // this is usually 33, the offset to convert from the ascii in a fastq code, to quality value in standard Sanger format    
    int offset = seq->qual_offset;    

	if (offset == 0) {
		offset = 64;	// Ilumina fastq format, default to keep backwards compatibility. 
	}    
    
	if (fp == NULL) {
		fputs("File not defined\n", stderr);
		exit(1);
	}
    
	if (seq == NULL) {
		fputs("Cannot pass a NULL pointer for seq\n", stderr);
		exit(1);
	}
    
	if (seq->seq == NULL) {
		fputs
        ("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->name == NULL) {
		fputs
        ("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->qual == NULL) {
		fputs
        ("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	boolean end_of_file = false;
    
	do {
        
		good_read = true;
		q = 0;
		j = 0;
        
		// read name of fastq entry
		if (fgets(line, LINE_MAX, fp) != NULL) {
            //if(fscanf(fp, "%s", line)){
			if (line[0] == '@') {
				for (i = 1; i < LINE_MAX; i++) {
					if (line[i] == '\n' /*|| line[i] == ' '*/
					    || line[i] == '\t'
					    || line[i] == '\r') {
						break;
					}
					if (i > LINE_MAX) {
						fputs("Name too long\n",
						      stderr);
						exit(1);
					}
					seq->name[i - 1] = line[i];
				}
				seq->name[i - 1] = '\0';
                
				// read sequence 
                
				while (fgets(line, LINE_MAX, fp) != NULL) {
                    
					if ((line[0] == '+') || (line[0] == '-')) {	// go to
						// get
						// qualities
						break;
					}
					// check line is fine
					for (i = 0; i < LINE_MAX; i++) {
						if (line[i] == '\n'
						    || line[i] == ' '
						    || line[i] == '\t'
						    || line[i] == '\r') {
							break;	// fine but nothing to add
						}
                        
						if(line[i] == '.'){//In case you have dots instad of Ns, that come from converting a qseq to fastq in a "quick and dirty" way
                            line[i] = 'N';
						}
						
						if (!nucleotide_good_base(line[i])) {
							good_read = false;
							fprintf(stderr,
                                    "Invalid symbol [%c] pos:%i in entry %s - skip read\n",
                                    line[i], i,
                                    seq->name);
						}
                        
						seq->seq[j] = line[i];
						j++;
                        
						if (j == max_read_length) {
							fprintf(stdout,
                                    "read [%s] too long [%i]. Exiting...\n",
                                    seq->name, j);
                            assert(NULL);
							exit(1);
						}
                        
					}
				}
                
				// read qualities -- verify position first
				file_pointer = ftell(fp);
                
				while (fgets(line, LINE_MAX, fp) != NULL) {
                    
					if (line[0] == '@' && j > 0 && (j <= q))	// then we have
						// gone on to the
						// next read 
						// allowing q>j in case where qualities longer
						// than j
					{
						fseek(fp, file_pointer,
						      SEEK_SET);
						break;	// goto next read
					}
                    
					for (i = 0; i < LINE_MAX; i++) {
						if (line[i] == '\n'
						    || line[i] == ' '
						    || line[i] == '\t'
						    || line[i] == '\r') {
							break;	// fine but nothing to add
						}
                        
						seq->qual[q] = line[i] - offset;
                        
                        if (seq->check_quality_values) {
                            if (seq->qual[q] < seq->lowest_expected_value) {
                                log_and_screen_printf("Warning: Quality [%d] for [%s] lower than expected for specified quality offset [%d-%d].\n", seq->qual[q], seq->name, seq->lowest_expected_value, seq->highest_expected_value);
                                seq->check_quality_values = false;
                            } else if (seq->qual[q] > seq->highest_expected_value) {
                                log_and_screen_printf("Warning: Quality [%d] for [%s] higher than expected for specified quality offset [%d-%d].\n", seq->qual[q], seq->name, seq->lowest_expected_value, seq->highest_expected_value);
                                seq->check_quality_values = false;
                            } else if (seq->qual[q] > seq->highest_raw_read_value) {
                                log_and_screen_printf("Warning: Quality [%d] for [%s] higher than expected for raw reads [%d].\n", seq->qual[q], seq->name, seq->highest_raw_read_value);
                                seq->check_quality_values = false;
                            }                                     
                        }
                        
						q++;
                        
						if (q == max_read_length) {
							fprintf(stdout,
                                    "qualities for [%s] longer than the max read length  [%i]. Exiting...\n",
                                    seq->name, q);
							exit(1);
						}
                        
					}
                    
					file_pointer = ftell(fp);
				}
                
				if (j != q) {
					fprintf(stdout,
                            "Lengths of quality [%i] and sequence [%i]  strings don't coincide for this read: [%s]. Skip it\n",
                            q, j, seq->name);
					good_read = false;
				}
                
			}	// if line starts with @
			else {
				fputs
                ("syntax error in fastq file -- it misses @\n",
                 stderr);
                assert(0);
				exit(1);
			}
		} else {
			end_of_file = true;
		}
        
	}
	while (!good_read && !end_of_file);
    
	seq->seq[j] = '\0';
	seq->qual[q] = '\0';	// this is not technically necessary but
	// simplifies many checks downstream
	seq->length = j;
	return j;
}
/*
 * Wrapper for the fastq parser, if fp is stdin, it will use an alternate parser
 */
int read_sequence_from_fastq(FILE * fp, Sequence * seq, int max_read_length){
    int ret;
    header_function * f;
    struct stat file_stat;
    fstat(fileno(fp), &file_stat);
    mode_t file_mode = file_stat.st_mode;
    
    if (S_ISFIFO(file_mode) || S_ISSOCK(file_mode)) {
        ret =  read_sequence_from_fastq_from_stream(fp, seq, max_read_length);
    }else{
        ret = read_sequence_from_fastq_from_file(fp, seq, max_read_length);
    }
    
    if (ret && seq->header != NULL) {
        
        f = seq->header;
        if (f->header_parser != NULL) {
            f->header_parser(seq);
        }
        
        
        
    }
    return ret;
}

/*
 * Read sequence from file "fp" in FASTA format. Read the qualities file
 * "fq" in Qual format (454) it returns the length of the sequence, 0 if
 * no sequence is left in file read with bad chars are skipped -- this is
 * a different behaviour to read_sequence_from_fasta 
 */

int read_sequence_from_fasta_and_qual(FILE * fp, FILE * fq, Sequence * seq, int max_read_length)
{
    
	char line[LINE_MAX];	// LINE_MAX is defined in limits.h
	int i;
	int j = 0;		// length of sequence
	int q = 0;		// length of qualities
	long file_pointer;
	boolean good_read = true;
	int offset = seq->qual_offset;
	if (offset == 0) {
		offset = 64;	// Ilumina fastq format, default to keep
        seq->qual_offset = offset;
		// backwards compatibility. 
	}
	// int offset = 64; // this is the offset to convert from the ascii in 
	// a fastq code, to quality value in standard Ilumina 1.3+ format
	// TODO parametize
	// int offset = 33; // this is the offset to convert from the ascii in 
	// a fastq code, to quality value in standard Sanger format
    
	if (fp == NULL) {
		fputs("File not defined\n", stderr);
		exit(1);
	}
    
	if (fq == NULL) {
		fputs("Quality file not defined\n", stderr);
		exit(1);
	}
    
	if (seq == NULL) {
		fputs("Cannot pass a NULL pointer for seq\n", stderr);
		exit(1);
	}
    
	if (seq->seq == NULL) {
		fputs
        ("Dont give me a null pointer for seq->seq - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->name == NULL) {
		fputs
        ("Dont give me a null pointer for seq->name - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	if (seq->qual == NULL) {
		fputs
        ("Dont give me a null pointer for seq->qual - alloc memory yourself and give me that\n",
         stderr);
		exit(1);
	}
    
	boolean end_of_file = false;
    
	do {
        
		good_read = true;
		q = 0;
		j = 0;
		file_pointer = ftell(fp);
		// read name of fastq entry
		if (fgets(line, LINE_MAX, fp) != NULL) {
			if (line[0] == '>') {
				for (i = 1; i < LINE_MAX; i++) {
					if (line[i] == '\n' || line[i] == ' '
					    || line[i] == '\t'
					    || line[i] == '\r') {
						break;
					}
					if (i > LINE_MAX) {
						fputs("Name too long\n",
						      stderr);
						exit(1);
					}
					seq->name[i - 1] = line[i];
				}
				seq->name[i - 1] = '\0';
				
				// read sequence 
				while (fgets(line, LINE_MAX, fp) != NULL) {
					if (line[0] == '>') {	// go to get qualities
						break;
					}
					// check line is fine
					for (i = 0; i < LINE_MAX; i++) {
						if (line[i] == '\n'
						    || line[i] == ' '
						    || line[i] == '\t'
						    || line[i] == '\r') {
							break;	// fine but nothing to add
						}
						if (!good_base(line[i])) {
							good_read = false;
							fprintf(stderr,
                                    "Invalid symbol [%c] pos:%i in entry %s - skip read\n",
                                    line[i], i,
                                    seq->name);
						}
						seq->seq[j] = line[i];
						j++;
						if (j == max_read_length) {
							fprintf(stdout,
                                    "read [%s] too long [%i]. Exiting...\n",
                                    seq->name, j);
							exit(1);
						}
                        
					}
				}
				// read qualities -- verify position first
				file_pointer = ftell(fq);
				while (fgets(line, LINE_MAX, fq) != NULL) {
                    
					if (line[0] == '>' )	// then we have
						// gone on to the
						// next read 
						// allowing q>j in case where qualities longer
						// than j
					{
						fgets(line, LINE_MAX, fq);
						if(j <= q)
                            break;
					}
					
					char * tok;//TODO use strtock_r
					char * exp = " ";
					tok = strtok(line, exp);
					while(tok != NULL) {
						
						
						seq->qual[q] = atoi(tok);
						printf("%s (%d)\n", tok, q );
						q++;
                        
						if (q == max_read_length) {
							fprintf(stdout,
                                    "qualities for [%s] longer than the max read length  [%i]. Exiting...\n",
                                    seq->name, q);
							exit(1);
						}
						tok = strtok(NULL, exp);
					}
                    
					file_pointer = ftell(fq);
				}
                
				if (j != q) {
					fprintf(stdout,
                            "Lengths of quality [%i] and sequence [%i]  strings don't coincide for this read: [%s]. Skip it\n",
                            q, j, seq->name);
					good_read = false;
				}
                
			}	// if line starts with @
			else {
				fputs
                ("syntax error in fastq file -- it misses >\n",
                 stderr);
				exit(1);
			}
		} else {
			end_of_file = true;
		}
        
	}
	while (!good_read && !end_of_file);
    
	seq->seq[j] = '\0';
	seq->qual[q] = '\0';	// this is not technically necessary but
	// simplifies many checks downstream
	seq->length = j;
	return j;
}



/**
 * Wrapper to alloc_sequence to make the allocations shorter.
 * The usual care of pointers can't be forgoten
 */
Sequence *sequence_new(int max_read_length, int max_name_length, char offset)
{
	Sequence *tmp = calloc(1, sizeof(Sequence));
    if (tmp == 0) {
        printf("Error: failed to allocated memory for sequence.\n");
        exit(-1);
    }    
	alloc_sequence(tmp, max_read_length, max_name_length, offset);
	return tmp;
}

/**
 * WARNING The pointer has to be preallocated where you are calling it. However, 
 * when calling the destroy method, you send the pointer to the pointer itself. 
 * 
 * It would be nice to have the "constructor" allocating the pointer.  
 */
void alloc_sequence(Sequence * seq, int max_read_length, int max_name_length, char offset)
{
    
	if (seq == NULL) {
		fputs("Cannot pass a null seq to alloc_sequence\n", stderr);
		exit(1);
	}
	seq->name = malloc(sizeof(char) * max_name_length);
	if (seq->name == NULL) {
		fputs("Out of memory trying to allocate string\n", stderr);
		exit(1);
	}
	seq->seq = malloc(sizeof(char) * (max_read_length + 1));
	if (seq->seq == NULL) {
		fputs("Out of memory trying to allocate string\n", stderr);
		exit(1);
	}
	// printf("Allocated seq! %p\n", seq->seq);
	seq->qual = malloc(sizeof(char) * (max_read_length + 1));
	if (seq->qual == NULL) {
		fputs("Out of memory trying to allocate string\n", stderr);
		exit(1);
	}
	seq->seq[max_read_length] = '\0';
	seq->qual[max_read_length] = '\0';
	seq->max_length = max_read_length;
	seq->max_name_length = max_name_length;
    seq->check_quality_values = false;
    seq->header = NULL;
    sequence_set_quality_parameters(seq, offset);
    
} 

boolean good_base(char c)
{   
    boolean ret = true; //We can have sequences of different stuff. If we want to validate this, we shall generate different "good base", according to the type of sequences we are reading. 
	/*boolean ret;
     if (c != 'A' && c != 'a' &&
     c != 'C' && c != 'c' &&
     c != 'G' && c != 'g' && c != 'T' && c != 't' && c != 'N'
     && c != 'n') {
     ret = false;
     } else {
     ret = true;
     }*/
    
	return ret;
}

boolean nucleotide_good_base(char c)
{   
    boolean ret = true; //We can have sequences of different stuff. If we want to validate this, we shall generate different "good base", according to the type of sequences we are reading. 
	
     if (c != 'A' && c != 'a' &&
         c != 'C' && c != 'c' &&
         c != 'G' && c != 'g' && c != 'T' && c != 't' && c != 'N'
         && c != 'n') {
            ret = false;
     } else {
         ret = true;
     }
    
	return ret;
}

void free_sequence(Sequence ** sequence)
{
	// fprintf(stderr, "sequence free: %p\n", &((*sequence)));
	free((*sequence)->name);	// TODO! make it safe, there are strings
	// lost somewhere!
	free((*sequence)->seq);
	free((*sequence)->qual);
	free(*sequence);
	*sequence = NULL;
}

void shift_last_kmer_to_start_of_sequence(Sequence * sequence, int length,
                                     short kmer_size)
{
    
	int i;
    
	if (length - kmer_size < kmer_size) {
		puts("kmer_size too long\n");
        assert(false);
		exit(1);
	}
    
	for (i = 0; i < kmer_size; i++) {
		sequence->seq[i] = sequence->seq[length - kmer_size + i];
		sequence->qual[i] = '\0';
	}
    
}

static int din_index(int n1, int n2)
{
	return (n1 * 5) + n2;
}

void clean_stats(SequenceStats * total)
{
	int i, j, k;
	for (i = 0; i < total->max_length; i++) {
		for (j = 0; j < 5; j++) {
			total->nucleotides[j][i] = 0;
			total->qual_nucleotides[j][i] = 0;
			for (k = 0; k < 5; k++) {
				total->dinculeotides[din_index(j, k)][i] = 0;
				total->qual_dinculeotides[din_index(j, k)][i] =
                0;
			}
		}
		total->total_bases = 0;
	}
	total->max_length = 0;
	total->reads = 0;
}

void sequence_remove_low_quality(Sequence * seq, char threshold)
{
	int i;
	for (i = 0; i < seq->length; i++) {
		if (seq->qual[i] < threshold)
			seq->seq[i] = 'N';
        
	}
    
}

void sequence_stats(SequenceStats * total, Sequence * sequence)
{
	int i, j, k;
	total->reads++;
	if (sequence->length > total->max_length) {
		total->max_length = sequence->length;
	}
	for (i = 0; i < sequence->length; i++) {
        
		boolean added = false;
        
		for (j = 0; j < 5; j++) {
			if (sequence->seq[i] == BASES[j]
			    || sequence->seq[i] == BASES[j + 5]
			    || (j == 4 && !added)) {
                
				total->nucleotides[j][i]++;
				total->qual_nucleotides[j][i] +=
                (long long)sequence->qual[i];
				added = true;
				boolean din_added = false;
				for (k = 0; k < 5; k++) {
                    
					if (i < sequence->length - 1) {	// looking only
						// where the
						// dinucleotide
						// exists. 
                        
						if (sequence->seq[i + 1] ==
						    BASES[k]
						    || sequence->seq[i + 1] ==
						    BASES[k + 5]
						    || (k == 4 && !din_added)) {
							// printf("Dinucleotide Found:
							// %c%c(%d)\n",BASES[j],BASES[k],i);
							total->dinculeotides
                            [din_index(j, k)]
                            [i]++;
							total->qual_dinculeotides
                            [din_index(j, k)][i]
                            +=
                            (sequence->qual[i] +
                             sequence->qual[i +
                                            1]);
							din_added = true;
						}
                        
					}
				}
			}
		}
        
		/*
		 * if(!added){ total->nucleotides[j][i]++;
		 * total->qual_nucleotides[j][i] += (long long) sequence->qual[i];
		 * } 
		 */
		total->total_bases++;
	}
    
}

void sequence_clean(Sequence * seq)
{
	seq->name[0] = '\0';
	seq->start = 0;
	seq->end = 0;
	seq->length = 0;
	seq->seq[0] = '\0';
	seq->qual[0] = '\0';
	seq->upper_case = false;
}

void sequence_append(char *s, char *q, Sequence * seq)
{
	long len = strlen(s);
	if (seq->max_length < len + seq->length) {
		fprintf(stderr,
                "[sequence_append] The sequence where appending doesn't have enought space\n");
		exit(-1);
	}
	// TODO validate that all the values are consistent... soft
	// validation?
    
	strcpy(&seq->seq[seq->end], s);
	strcpy(&seq->qual[seq->end], q);
	seq->length += (int)len;
	if (seq->start == 0) {
		seq->start = 1;
		seq->end = 1;
	}
	seq->end += (int)len;
    
}

void sequence_set_name(char *name, Sequence * seq)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (seq == NULL) {
		fprintf(stderr, "[sequence_set_name] Null sequence\n");
		exit(-1);
	}
	if (name == NULL) {
		fprintf(stderr, "[sequence_set_name]Null name\n");
		exit(-1);
	}
#endif
	long len = strlen(name);
	if (len > seq->max_name_length) {
		fprintf(stderr,
                "[sequence_set_name]The name %s  is to big (Max lenght name:%d\n",
                name, seq->max_name_length);
		exit(-1);
	}
    
	strcpy(seq->name, name);
    
}

void sequence_append_name(char *name, Sequence * seq)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (seq == NULL) {
		fprintf(stderr, "[sequence_append_name] Null sequence\n");
		exit(-1);
	}
	if (name == NULL) {
		fprintf(stderr, "[sequence_append_name]Null name\n");
		exit(-1);
	}
#endif
	long len = strlen(name) + strlen(seq->name);
	if (len > seq->max_name_length) {
		fprintf(stderr,
                "[sequence_append_name]The name %s+%s  is to big (Max lenght name:%d\n",
                seq->name, name, seq->max_name_length);
		exit(-1);
	}
    
	strcat(seq->name, name);
    
}

void sequence_to_upper_case(Sequence * seq)
{
	int i;
	for (i = 0; i < seq->length; i++) {
		seq->seq[i] = toupper(seq->seq[i]);
	}
	seq->upper_case = true;
}

void sequence_add_base(char s, char q, Sequence * seq)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (seq == NULL) {
		fprintf(stderr, "[sequence_add_base] The sequence is NULL");
		exit(-1);
	}
#endif
	if (s == '\0')
		return;
    
	seq->seq[seq->length] = s;
	seq->qual[seq->length] = q;
    
	seq->length++;
    
	seq->seq[seq->length] = '\0';
	seq->qual[seq->length] = '\0';
}

/**
 * TODO... write an implentation that reverses seq and qual in the same  step... 
 */
static char *strrev(char *s, int n)
{
	int i = 0;
	while (i < n / 2) {
		*(s + n) = *(s + i);	// uses the null character as the
		// temporary storage.
		*(s + i) = *(s + n - i - 1);
		*(s + n - i - 1) = *(s + n);
		i++;
	}
	*(s + n) = '\0';
    
	return s;
}

static char *nuc_comp(char *s, int n)
{
	int i = 0;
	while (i < n) {
		switch (s[i]) {
            case 'A':
                s[i] = 'T';
                break;
            case 'T':
                s[i] = 'A';
                break;
            case 'C':
                s[i] = 'G';
                break;
            case 'G':
                s[i] = 'C';
                break;
		}
		i++;
	}
	*(s + n) = '\0';
    
	return s;
}

void sequence_reverse(Sequence * seq)
{
	strrev(seq->seq, seq->length);
	strrev(seq->qual, seq->length);
}

void sequence_copy(Sequence * to, Sequence * from)
{
	strcpy(to->seq, from->seq);
	strcpy(to->name, from->name);
	int i;
	for (i=0; i<from->max_length; i++) {
		to->qual[i]=from->qual[i];
	}
	//strcpy(to->qual, from->qual);
	//printf("qual: %s = %s", to->qual, from->qual);
	to->start = from->start;
	to->end = from ->end;
#ifdef ALIGN
	to->count = from ->count;
#endif
	to->length = from ->length;
	to->qual_offset = from->qual_offset;
}

void sequence_complement(Sequence * seq)
{
	nuc_comp(seq->seq, seq->length);
}

void sequence_reverse_complement(Sequence * seq)
{
	sequence_reverse(seq);
	sequence_complement(seq);
}

double base_content(Nucleotide n, SequenceStats * total){
	long long count = 0;
	int i;
	for (i = 0; i<total->max_length; i++){
        count += total->nucleotides[n][i];
	}
	double ret = (((double)count / (double)total->total_bases)*100 ); 
	printf("%lld / %lld = %f\n",count, total->total_bases, ret);
	return (double)ret;
}

void print_stats(FILE * f, SequenceStats * total)
{
	int i, j, k;
	long long count_n, count_d, qual_n, qual_d;
	fprintf(f, "TYPE\tINDEX\tID\tTOTAL\tAVG\n");
	fprintf(f, "reads\t0\t0\t%d\t%5.2f\n", total->reads,
            (double)total->total_bases / total->reads);
	for (j = 0; j < 5; j++) {	// Iterate on the nucleotide
		count_n = 0;
		qual_n = 0;
		for (k = 0; k < 5; k++) {	// Iterate on the dinucleotide
			count_d = 0;
			qual_d = 0;
			for (i = 0; i < total->max_length; i++) {	// iterate on the
				// position 
				if (k == 0) {
					count_n += total->nucleotides[j][i];
					qual_n += total->qual_nucleotides[j][i];
					fprintf(f,
                            "nuc_pos\t%d\t%c\t%lld\t%5.2f\n",
                            i, BASES[j],
                            total->nucleotides[j][i],
                            100 *
                            (double)total->nucleotides[j][i]
                            / total->total_bases);
					fprintf(f,
                            "nuc_qual_pos\t%d\t%c\t%lld\t%5.2f\n",
                            i, BASES[j],
                            total->qual_nucleotides[j][i],
                            (double)
                            total->qual_nucleotides[j][i] /
                            total->nucleotides[j][i]);
				}
				count_d +=
                total->dinculeotides[din_index(j, k)][i];
				qual_d +=
                total->qual_dinculeotides[din_index(j, k)]
                [i];
                
				fprintf(f,
                        "din_pos\t%d\t\"%c.%c\"\t%lld\t%5.2f\n",
                        i, BASES[j], BASES[k],
                        total->dinculeotides[din_index(j, k)]
                        [i],
                        700 *
                        (double)
                        total->dinculeotides[din_index(j, k)][i]
                        / total->total_bases);
				fprintf(f,
                        "din_qual_pos\t%d\t\"%c.%c\"\t%lld\t%5.2f\n",
                        i, BASES[j], BASES[k],
                        total->qual_dinculeotides[din_index
                                                  (j, k)][i],
                        (double)
                        total->qual_dinculeotides[din_index
                                                  (j, k)][i]
                        /
                        total->dinculeotides[din_index(j, k)]
                        [i]);
			}
			fprintf(f, "din\t0\t\"%c.%c\"\t%lld\t%5.2f\n", BASES[j],
                    BASES[k], count_d,
                    (double)100 * count_d / total->total_bases);
			fprintf(f, "din_qual\t0\t\"%c.%c\"\t%lld\t%5.2f\n",
                    BASES[j], BASES[k], qual_d,
                    (double)qual_d / count_d);
		}
		fprintf(f, "nuc\t0\t%c\t%lld\t%5.2f\n", BASES[j], count_n,
                (double)100 * count_n / total->total_bases);
		fprintf(f, "nuc_qual\t0\t%c\t%lld\t%5.2f\n", BASES[j], qual_n,
                (double)qual_n / count_n);
	}
}

char sequence_get_base(int pos, Sequence * seq)
{
	assert(seq  != NULL);
	assert(pos >= 0);
	assert(pos  < seq->length);
	return seq->seq[pos];
}

void sequence_set_base(char base, int pos, Sequence * seq)
{
	assert(seq  != NULL);
	assert(pos >= 0);
	assert(pos  < seq->length);
	seq->seq[pos] = base;
}

char sequence_get_qual(int pos, Sequence * seq)
{
    
	assert(seq  != NULL);
	assert(pos >= 0);
	assert(pos  < seq->length);
	
	return seq->qual[pos];
}

int sequence_get_length(Sequence * seq)
{
	assert(seq != NULL);
	return seq->length;
}

SequenceArray *sequence_array_new(int capacity)
{
	SequenceArray *tmp = calloc(1, sizeof(SequenceArray));
	tmp->sequences = calloc(capacity, sizeof(Sequence *));
	tmp->capacity = capacity;
	tmp->total = 0;
	return tmp;
}

void
sequence_array_add(int max_read_length, int max_name_length, SequenceArray * sa)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (sa == NULL) {
		fprintf(stderr,
                "[sequence_array_add] Can't destroy a null array!\n");
		exit(-1);
	}
#endif
	if (sa->total >= sa->capacity) {
		fprintf(stderr,
                "[sequence_array_add] Not enought memory preallocated for reference sequences, please ensure you had set up the parameter -R with a number at least equal to the number of reference sequences on the file -r \n");
		exit(-1);
	}
	sa->sequences[sa->total] =
    sequence_new(max_read_length, max_name_length, 0);
	sa->total++;
}

Sequence *sequence_array_get_sequence(int pos, SequenceArray * sa)
{
	assert(sa != NULL);
    //	assert(pos < sa->total);
	return sa->sequences[pos];
}

/**
 * Destroys all the sequences in the array. The idea behind
 * is that there shall not be any sequence created outside the 
 * encapsulated methods. 
 */
void sequence_array_destroy(SequenceArray ** sa)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (*sa == NULL || (*sa) == NULL) {
		fprintf(stderr,
                "[sequence_array_destroy] Can't destroy a null array!\n");
		exit(-1);
	}
#endif
	SequenceArray *arr = (*sa);
    
	int i = 0;
#ifndef NO_BOUNDS_CHECK
	for (i = 0; i < arr->capacity; i++) {
		if (i > arr->total && arr->sequences[i] != NULL) {
			fprintf(stderr,
                    "[sequence_array_destroy] There is a pointer of a sequence outside the bounds!\n");
			exit(-1);
		}
	}
#endif
	for (i = 0; i < arr->total; i++) {
		if (arr->sequences[i] != NULL) {
			free_sequence(&(arr->sequences[i]));
		}
	}
	free(arr);
	(*sa) = NULL;
    
}

void sequence_array_clean(SequenceArray * sa)
{
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (sa == NULL) {
		fprintf(stderr,
                "[sequence_array_clean] Can't destroy a null array!\n");
		exit(-1);
	}
#endif
    
	int i;
	for (i = 0; i < sa->total; i++) {
		if (sa->sequences[i] != NULL) {
			sequence_clean(sa->sequences[i]);
		}
	}
}

void sequence_mask(int from, int to, Sequence * seq){
	int i;
	for(i = from;i<to; i++){
		seq->seq[i]='N';	
	}
}

void sequence_print_fasta(FILE * f, Sequence * seq)
{
	assert(seq!=NULL);
    
	if(strlen(seq->seq)>0){
		fprintf(f, ">%s\n%s\n", seq->name, seq->seq);
		fflush(f);
	}
}

void sequence_print_fasta_subseq(FILE * f,int starta, int enda, Sequence * seq)
{
	assert(seq!=NULL);
	assert(starta >= 0);
	//assert(end <= seq->length);
	fprintf(f, ">%s(%d, %d)\n", seq->name, starta, enda);
	int i;
	for (i = 0; i < seq->length; i++) {
		if(i<starta || i > enda){
			fprintf(f, "-");
		}else {
			fprintf(f, "%c",seq->seq[i]);
		}
        
		
		if(i % 81==0){
			fprintf(f, "\n");
		}
	}
	fprintf(f, "\n");
	
}

void sequence_print_fastq(FILE * f, Sequence * seq)
{
    
#ifndef NO_BOUNDS_CHECK		// For debugging purposes, in
	// "production", if coded properly, we can 
	// skip the validation
	if (seq == NULL) {
		fprintf(stderr,
                "[sequence_print_fastq] Can't print a null sequence!\n");
		exit(-1);
	}
#endif
	fprintf(f, "@%s\n%s\n+\n", seq->name, seq->seq);
	int i;
	for (i = 0; i < seq->length; i++) {
		fprintf(f, "%c", seq->qual[i] + seq->qual_offset);
	}
	fprintf(f, "\n");
}

void sequence_iterator(void (*f) (char, int), Sequence * seq)
{
	int i;
	for (i = 0; i < seq->length; i++) {
		f(seq->seq[i], i);
	}
}

int sequence_count_gaps(Sequence * seq, int max)
{
	int i;
	int count = 0;
	for (i = 0; i < max; i++) {
		if (seq->seq[i] == '-') {
			count++;
		}
	}
	return count;
}

int sequence_count_homopolymer(boolean forward, int pos, Sequence * seq)
{
    
	signed int i = 0;
	while (seq->seq[pos] == seq->seq[pos + i]) {
		if (forward)
			i++;
		else
			i--;
	}
	return abs(i);
    
}

void sequence_remove_missing_last_bases(Sequence * seq){
	boolean removed = false;
	char base;
	do{
		removed = false;
		base = sequence_get_base(seq->length - 1, seq);
		if(base == '-'){
			removed = true;
			seq->length--;	
			seq->seq[seq->length]  = '\0';
			seq->qual[seq->length] = '\0';
		}
		
		
	}while(removed);
}

/**
 * Removes the base at position pos
 */
void sequence_remove_base(int pos, Sequence * seq){
	int i;
	char *s = seq->seq;
	char *q = seq->qual;
	for(i = pos; i < seq->length; i++){
		s[i]=s[i+1];
		q[i]=q[i+1];
	} 
	seq->length--;
}

/**
 * Trims the sequence seq to the length l. If
 * the sequence is shorter than the lenght l, no
 * action is done. 
 * 
 * To avoid looping too much, it just sets
 * the last character to null and changes the length
 */
void sequence_trim(int l, Sequence * seq){
	
	if(sequence_get_length(seq) < l) return;
	
	//int i;
	char *s = seq->seq;
	char *q = seq->qual;
	
	
	s[l]='\0';
	q[l]='\0';
	
	seq->length = l;
}


/*
 void sequence_remove_base(int pos, Sequence * seq){
 int i;
 char *s = seq->seq;
 for(i = pos; i < seq->length; i++){
 s[i]=s[i+1];
 } 
 }*/

void sequence_remove_base_up_to_limit(int pos,int limit ,Sequence * seq){
	int i;
	char *s = seq->seq;
	for(i = pos; i < limit; i++){
		s[i]=s[i+1];
	} 
	s[i] = '-';
}

void sequence_insert_base_up_to_limit(char base, int pos,int limit ,Sequence * seq){
	int i;
	char *s = seq->seq;
	char old;
	char new = base;
	for(i = pos; i < limit; i++){
		old = s[i];
		s[i] = new;
		new = old;
	} 
}

int sequence_prev_anchor_base(int pos, Sequence * seq){
	int i;
	char *s = seq->seq;
	for(i = pos; i > 0; i--){
		if(s[i]!='N'){
			return i;
		}
	} 
	return 0;
}

int sequence_next_anchor_base(int pos, Sequence * seq){
	int i;
	char *s = seq->seq;
	for(i = pos; i < seq->length; i++){
		if(s[i]!='N'){
			return i;
		}
	} 
	return 0;
}

/**
 * Returns -1 when no next homopolymer is found before the limit of the search
 */
int sequence_next_hompoplymer(int pos, int max_righ, Sequence *seq){
	int i;
	
	for(i = pos; i > 0 && i > max_righ; i++){
		if(sequence_count_homopolymer(true, i, seq) > 0){
			return i;
		}
	} 
	return -1;
}

/**
 * Returns -1 when no previous homopolymer is found before the limit of the search
 */
int sequence_prev_hompoplymer(int pos, int max_left, Sequence *seq){
	int i;
	
	for(i = pos; i > 0 && i > max_left; i--){
		if(sequence_count_homopolymer(false, i, seq) > 0){
			return i;
		}
	} 
	return -1;
}

/**
 * Returns true if the base is valid on the code. 
 * The code can be any  IUPAC
 */
boolean base_is_valid(char base, char code){
	
	base = toupper(base);
	code = toupper(code);
	
	if(base == code)
		return true;
    
	base = base=='U'?'T':base;
    
	switch (code){
		case 'N':
			return true;
		case 'R':
			return (base == 'G' || base == 'A')?true:false;
		case 'Y':
			return (base == 'T' || base == 'C')?true:false;
		case 'K':
			return (base == 'G' || base == 'T')?true:false;
		case 'M':
			return (base == 'A' || base == 'C')?true:false;		
		case 'S':
			return (base == 'G' || base == 'C')?true:false;
		case 'W':
			return (base == 'A' || base == 'T')?true:false;
		case 'B':
			return (base == 'A')?false: true;
		case 'D':
			return (base == 'C')?false:true;
		case 'H':
			return (base == 'G')?false:true;
		case 'V':
			return (base == 'T')?false:true;
			
	}
	return false;
    
}

boolean base_is_unambiguous(char base){
	switch(base){
		case 'A':
		case 'C':
		case 'T':
		case 'G':
		case 'U':
			return true;
		default:
			return false;
	};
}
/**
 * Compares two sequences alphabetically, however it allows for 
 */ 
int sequence_compare_with_ambiguity(Sequence * seq1, Sequence * seq2){
	assert(seq1 != NULL);
	assert(seq2 != NULL);
	
	int i;
	int len1 = sequence_get_length(seq1);
	int len2 = sequence_get_length(seq2);
	boolean unambig1, unambig2;
	int diff = 0;
	char base1, base2;
	
	for(i = 0; diff == 0 &&i < len2 && i < len1; i++){
		
		base1 = sequence_get_base(i, seq1);
		base2 = sequence_get_base(i, seq2);
		
		unambig1 = base_is_unambiguous(base1);
		unambig2 = base_is_unambiguous(base2);		
		
		if(unambig1 && unambig2){ //simple comparation, both are unambiguous
			diff = base2 - base1;
		}else if(unambig1 && !unambig2){//First one is unambigous, but second one is not 
			if(base_is_valid(base1, base2)){
				diff = 0;
			}else{
				diff = base2 - base1;
			}
		}else if(!unambig1 && unambig2){//First one is unambigous, but second one is not 
			if(base_is_valid(base1, base2)){
				diff = 0;
			}else{
				diff = base2 - base1;
			}	
		}else{//Both are ambiguos just ignore it... diff is the same than before. Would be nice to check if there is a common for both, but seems like at the moment this is not happening. 
			
		}
	}
	if(diff == 0 && len1 != len2){//In case the sequences are the same up to certain base, sort by length
		diff = len2 - len1;
	}
	return diff;
}

/**
 * 
 * Merges two sequences by removing the ambiguity codes. If they are not the same length or 
 * they are not "equal", the function aborts the program. 
 * At the end, both sequences should be the same.
 */ 
void sequence_merge_removing_ambiguity(Sequence * seq1, Sequence * seq2){
	assert(seq1 != NULL);
	assert(seq2 != NULL);
	
	int i;
	int len1 = sequence_get_length(seq1);
	int len2 = sequence_get_length(seq2);
	boolean unambig1, unambig2;
    char base1, base2;
	
	assert(len1 == len2); 
	
	for(i = 0;  i < len1 ; i++){
		
		base1 = sequence_get_base(i, seq1);
		base2 = sequence_get_base(i, seq2);
		
		unambig1 = base_is_unambiguous(base1);
		unambig2 = base_is_unambiguous(base2);		
		
		if(unambig1 && unambig2){ //simple comparation, both are unambiguous
			assert(base1 == base2); //For debugging pruposes, validate the bases are equal. 
		}if(unambig1 && !unambig2){//First one is unambigous, but second one is not 
			assert(base_is_valid(base1, base2));
			sequence_set_base(base1, i, seq2);  //We asume the sequences are already equal, hence we can just asign them. 
		}else if(!unambig1 && unambig2){//First one is unambigous, but second one is not 
			sequence_set_base(base2, i, seq1);
		}
	}
	
	assert(strcmp(seq1->seq, seq2->seq) == 0);//Validate that both sequences are the same after the merging.
	
}


int sequence_differences_with_mask(Sequence * seq, Sequence * ref){
	int i;
	int count = 0;
	for(i = 0; i < sequence_get_length(seq); i++){
		if(!base_is_valid(sequence_get_base(i,seq),sequence_get_base(i,ref))){
			count++;
		}
	} 
	return count;
}

