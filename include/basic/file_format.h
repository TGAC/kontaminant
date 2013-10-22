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

#ifndef FILE_FORMAT_H_
#define FILE_FORMAT_H_

typedef enum
 {
   UNSPECIFIED_FORMAT   = 0,
   FASTA                = 1,
   FASTQ                = 2,
   CTX                  = 3,
   ROCHE                = 4,
   HASH                 = 5,
   CSFASTA              = 6,
   KMERS                = 7, 
   FILE_FORMAT_LAST     = 8,
 } FileFormat ;

typedef enum {
    UNKNOWN_HEADER      = 0,
    CASAVA_1_8          = 1,
    UNKNOWN_HEADER_LAST = 2
} sequence_header_type;


char * file_format_to_string(FileFormat ff);
FileFormat string_to_file_format(char * format);

char * sequence_header_type_to_string(sequence_header_type sht);
sequence_header_type string_to_sequence_header_type(char * format);

#endif //FILE_FORMAT_H_
