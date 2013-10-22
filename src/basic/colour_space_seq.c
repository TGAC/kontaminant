#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <binary_kmer.h>
#include <seq.h>
#include <colour_space_seq.h>



//Returns the base in upper case and validates it is a valid base. 
static char validate_base(char c){
    
    if(c == 'n'){
        c = 'N';
    }else if(c == 'a'){
        c = 'A';
    }else if(c == 'c'){
        c = 'C';
    }else if(c == 't'){
        c = 'T';
    }else if(c == 'g'){
        c = 'G';
    }
    assert(c == 'A'  || c == 'T' || c == 'G' || c == 'C' || c == 'N');
    return c;
}


static char base_complement(char base){
    switch (base){
        case 'A' :
        return 'T';
        case 'C' :
        return 'G';
        case 'G' :
        return 'C';
        case 'T' :
        return 'A';
        case 'N':
        return 'N';
    }
    assert(false);//Should not come here... We are not supporting any other code. 
}

static boolean is_complement(char base1, char base2){
    return base_complement(base1) == base2;
}

static char base_transversion(char base){
    switch (base){
        case 'A' :
           return 'C';
        case 'C' :
           return 'A';
        case 'G' :
           return 'T';
        case 'T' :
           return 'G';
        case 'N':
           return 'N';
        default:
            return 0;
       }
}

static boolean is_transversion(char base1, char base2){
    return base_transversion(base1) == base2;
}
//A ↔ G  C ↔ T 
static char base_transition(char base){
    switch (base){
        case 'A' :
            return 'G';
        case 'C' :
            return 'T';
        case 'G' :
            return 'A';
        case 'T' :
            return 'C';
        case 'N':
            return 'N';
        default:
            return 0;
    }
}

static boolean is_transition(char base1, char base2){
    return base_transition(base1) == base2;
}


static char base_to_colour(char first, char second){
    first = validate_base(first);
    second = validate_base(second);
    if(first == 'N' || second == 'N'){
        return '.';
    }else if(first == second){
        return '0';
    }else if(is_transversion(first, second)){
        return '1';
    }else if (is_transition(first, second)){
        return '2';
    }else if(is_complement(first, second)){
        return '3';
    }else{
        assert(false);
        return 0;//Should never ever be here. 
    }
    
}

static char colour_to_base(char base, char colour){
    switch (colour){
        case '0':
            return base;
        case '1':
            return base_transversion(base);
        case '2':
            return base_transition(base);
        case '3':
            return base_complement(base);
        case '.':
            return 'N';
        default:
            assert(false);
            return 0; //Should never ever come here. 
    }
}


/**
* Converts a sequence in colour space to base space. 
*/
Sequence * cs_sequence_to_base_sequence(Sequence * seq_cs, Sequence * seq_base){
    int i;
    int len;
    char current_base;
    assert(seq_cs != NULL);
    assert(seq_base != NULL);
    
    len = sequence_get_length(seq_cs);
    current_base =  sequence_get_base(0, seq_cs);
    sequence_clean(seq_base);
    for(i = 1; i < len; i++  ){
      //  fprintf(stderr, "%c:", current_base);
      //  fprintf(stderr, "%c->", sequence_get_base(i, seq_cs));
        current_base = colour_to_base(current_base, sequence_get_base(i, seq_cs));
    //    fprintf(stderr, "%c\n", current_base);
        sequence_add_base(current_base, sequence_get_qual(i,seq_cs), seq_base);//TODO: AVG the pair base?
    }
    return seq_base;
}

/**
* Converts a sequence in base space to colour space. 
*/
Sequence * base_sequence_to_cs_sequence(Sequence * seq_base, Sequence * seq_cs){
    assert(seq_cs != NULL);
    assert(seq_base != NULL);
    
    int i;
    int len = sequence_get_length(seq_base);
//    printf("len %d\n", len);
    char current_base = 'T'; 
    char next_base;
    char current_colour;
    sequence_clean(seq_cs);
    sequence_add_base(current_base, 40, seq_cs);
    for(i = 0; i < len; i++){
        //fprintf(stderr, "%c:", current_base);
        next_base = sequence_get_base(i, seq_base);
       // fprintf(stderr, "%c->", next_base);
        //current_base = colour_to_base(current_base, sequence_get_b
        current_colour = base_to_colour(current_base, next_base);
    //    fprintf(stderr, "%c\n", current_colour);
        sequence_add_base(current_colour, sequence_get_qual(i,seq_base), seq_cs);
        current_base = next_base;
    }
    return seq_cs;
    
}
//(2) A ↔ G  C ↔ T In genetics, a transition is a point mutation that changes a purine nucleotide to another purine (A ↔ G) or a pyrimidine nucleotide to another pyrimidine (C ↔ T)
//(1) A ↔ C  G ↔ T In molecular biology, transversion refers to the substitution of a purine for a pyrimidine or vice versa.[1
/*
static void start_solid_code(){
    //Make this reentrant by only allowing this to happen once with a lock.
    if(!color_code_started){
        
        solid_code[Adenine][0] = Adenine;
        solid_code[Adenine][1] = Cytosine;
        solid_code[Adenine][2] = Guanine; 
        solid_code[Adenine][3] = Thymine;
        
        solid_code[Cytosine][1] = Adenine;
        solid_code[Cytosine][0] = Cytosine;
        solid_code[Cytosine][3] = Guanine;
        solid_code[Cytosine][2] = Thymine;
        
        solid_code[Guanine][2] = Adenine;
        solid_code[Guanine][3] = Cytosine;
        solid_code[Guanine][0] = Guanine;
        solid_code[Guanine][1] = Thymine;
        
        solid_code[Thymine][3] = Adenine;
        solid_code[Thymine][2] = Cytosine;
        solid_code[Thymine][1] = Guanine;
        solid_code[Thymine][0] = Thymine;
    }
    
}
*/