#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <file_format.h>


    
static char fileFormatStrings[][12] = {"UNSPECIFIED", "FASTA", "FASTQ", "CTX", "ROCHE", "HASH", "CSFASTA",  "KMERS",};

static char sequenceHeaderStrings[][15] = {"UNKNOWN", "CASAVA_1.8"};

static void stringToUpperCase(char *string){
    int i;
    unsigned long  len =  strlen(string);
    for (i = 0; i < len ; i++){
        if (isalpha(string[i])){
            string[i] = toupper(string[i]);
        }

    }
}

char * file_format_to_string(FileFormat ff){
    return fileFormatStrings[ff];
}

//BEWARE! This transforms the format to upper case!
FileFormat string_to_file_format(char * format){
    stringToUpperCase(format);
    int i;
    for(i = 0; i < FILE_FORMAT_LAST; i++)
    {
        if (strcmp(format, fileFormatStrings[i]) == 0) {
            return i;
        }
    }
    return UNSPECIFIED_FORMAT;
}

char * sequence_header_type_to_string(sequence_header_type sht){
    return sequenceHeaderStrings[sht];

}
sequence_header_type string_to_sequence_header_type(char * format){
    return CASAVA_1_8;
}

