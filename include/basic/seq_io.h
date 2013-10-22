//
//  seq_io.h
//  Untitled
//
//  Created by Ricardo Ramirez-Gonzalez (TGAC) on 12/03/2012.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#ifndef seq_io_h
#define seq_io_h

#define MAX_FIELD_SIZE 20



typedef struct{
    void (* header_parser)(Sequence * seq); 
    char * (* get_index)(Sequence * seq); 
    char * instrument;
    int run_number;
    char * flowcell_id;
    int lane;
    int tile;
    int x_pos;
    int y_pos;
    int read;
    boolean filtered;
    int control_number;
    char * index;
}casava_sequence_header;

struct{
    char * filename;
    FILE * file;
    SequenceArray * sequences_read_1;
    SequenceArray * sequences_read_2;
    long long written;
    FileFormat format;
    int capacity;
    int size;
    boolean concurrent;
    boolean paird;
}SequenceBufferWriter;

void append_sequence(char * filename, Sequence * seq, FileFormat format);

void append_sequence_fh(FILE * filename, Sequence * seq, FileFormat format);

void * new_sequence_header(sequence_header_type header_type);

#endif
