//
//  seq_io.c
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 12/03/2012.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//


#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include <global.h>
#include <logger.h>
#include <nucleotide.h>
#include <file_format.h>
#include <seq.h>
#include <seq_io.h>



static void print_error_no(void) {
  switch (errno) {
            case EAGAIN:
                fprintf(stderr, "EAGAIN");
                break;
            case EBADF:
                fprintf(stderr, "EBADF");
                break;
            case EACCES:
                fprintf(stderr, "EACCES");
                break;
            case EINVAL:
                fprintf(stderr, "EINVAL");
            case ENOLCK:
                break;
                fprintf(stderr, "ENOLOCK");
                break;
            default:
                break;
        }
        fprintf(stderr, " error: %d\n" ,errno);

}
/*
 * Function that appends a sequence at the end of a file. 
 * It locks the file while it is writing with fcntl. 
 * In the future, I would like to add some handler mechanisms
 * to this, and perhaps a buffered writer, so it writes in chunks. 
 */

void append_sequence(char * filename, Sequence * seq, FileFormat format){
    
    FILE * output = NULL;
    output = fopen(filename, "a");
    if(output == NULL){
        log_and_screen_printf("Unable to open file: %s\n", filename);
        exit(-1);
    }
    
#ifdef LOCKING
    struct flock fl;// = {F_WRLCK, SEEK_SET,   0,      0,     0 };
    fl.l_len = 0;
    fl.l_start = 0;
    fl.l_whence = SEEK_SET;
    fl.l_type = F_WRLCK;
    fl.l_pid = getpid();
   short attemps = 10;
    
    while (fcntl(fileno(output), F_SETLKW, &fl) == -1 && attemps-- > 0) {
        print_error_no();

        //sleep(rand()%10);//Sleep from 0 to 10 seconds...
    }
#endif 
    
    switch (format) {
        case FASTA:
            sequence_print_fasta(output, seq);
            break;
        case FASTQ:
            sequence_print_fastq(output, seq);
            break;
        default:
            fprintf(stderr, "Format not implemented for writing\n");
            assert(false);
            exit(-1);
            break;
    }
#ifdef LOCKING
    fl.l_type=F_UNLCK;
  //  fl.l_whence = SEEK_END;
    if (output != NULL && fcntl(fileno(output), F_SETLKW, &fl) == -1 ) {
        
        log_and_screen_printf("PID:%d ERROR: Unable to unlock!\n", fl.l_pid);
        print_error_no();
    }
#endif
    if (output != NULL) {
        fclose(output);
    }
}


void append_sequence_fh(FILE * output, Sequence * seq, FileFormat format){
    assert(output != NULL);
    switch (format) {
        case FASTA:
            sequence_print_fasta(output, seq);
            break;
        case FASTQ:
            sequence_print_fastq(output, seq);
            break;
        default:
            fprintf(stderr, "Format not implemented for writing\n");
            assert(false);
            exit(-1);
            break;
    }

}



static void parse_casava_1_8_header(Sequence * seq){
 //Test line: @HWI-ST319:185:D0CWKACXX:3:1101:1189:2148 1:N:0:TGATAACA
    assert(seq!=NULL);
    assert(seq->header != NULL);
    casava_sequence_header * header = (casava_sequence_header *) seq->header; 
    
    char filtered = '-';
    sscanf(seq->name, "%[^:]:%d:%[^:]:%d:%d:%d:%d %d:%c:%d:%s", header->instrument, &header->run_number, header->flowcell_id, &header->lane, &header->tile, &header->x_pos, &header->y_pos, &header->read, &filtered, &header->control_number, header->index);
    header->filtered = filtered=='Y';
    
}

static char * get_index_from_casava_1_8(Sequence * seq){
    assert(seq!=NULL);
    assert(seq->header != NULL);
    casava_sequence_header * header = (casava_sequence_header *) seq->header; 
    return header->index;
}



static casava_sequence_header * new_casava_sequence_header_1_8(){
    casava_sequence_header * header = calloc(1, sizeof(casava_sequence_header));
    if (header == NULL) {
        fprintf(stderr, "Unable to allocate header;");
        assert(false);
        exit(-1);
    }
    
    //TODO: make this sizes dynamic? 
    header->instrument = calloc(MAX_FIELD_SIZE, sizeof(char*));
    header->flowcell_id = calloc(MAX_FIELD_SIZE, sizeof(char *));
    header->instrument = calloc(MAX_FIELD_SIZE, sizeof(char *));
    header->index = calloc(MAX_FIELD_SIZE, sizeof(char *));
    header->header_parser = &parse_casava_1_8_header;
    header->get_index = &get_index_from_casava_1_8;
    return header;
}

void * new_sequence_header(sequence_header_type header_type){
    switch (header_type) {
        case CASAVA_1_8:
            return  new_casava_sequence_header_1_8();
            break;
        case UNKNOWN_HEADER:
        case UNKNOWN_HEADER_LAST:
        default:
            log_and_screen_printf("Unknown format. \n");
            assert(false);
            break;
    }
    return NULL;
}
