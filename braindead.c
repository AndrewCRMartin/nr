/*************************************************************************

   Program:    nr
   File:       nr.c
   
   Version:    V1.2
   Date:       20.07.00
   Function:   Create a non-redundant sequence data set
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2000
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      a.c.r.martin@reading.ac.uk
               andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   Possible problems:
   ------------------
   It is possible that StoreSequenceFragment() may fail (returning it's
   warning message) in the rather bizarre case that (assuming a 3res 
   fragment) we have stored:
       ABCXXXXXXX : ABC
       BCDXXXXXXX : BCD
       CDEXXXXXXX : CDE
       DEPXXXXXXX : DEP
       EPQXXXXXXX : EPQ
       PQRXXXXXXX : PQR
       QRTXXXXXXX : QRS
       RSTXXXXXXX : RST
   and wish to store:
       ABCDEPQRST
   All the fragment keys have appeared before.
   With a reasonably large fragment size this is very unlikely and
   increasing the fragment size at run time should solve the problem.

   Possible speedups:
   ------------------
   One route to speedup would be to create a linked list of the keys
   in seqhash_temp and work through this list rather than through
   the keys directly so that the keys can be deleted from seqhash_temp
   as we find them. Thus DropSequence() would actually remove an item
   rather than building a list of items to delete and PurgeSequences()
   would no longer be needed.

   Possible changes:
   -----------------
   CompareSequences() could be enhanced to perform non-exact matching to
   allow a limited number of mutations

   Development Time:
   -----------------
   08.06.00-15.06.00     2days
   30.06.00              2hours
   11.07.00              2h
   12.07.00              4h
   19.07.00              2h
   20.07.00              5h

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  15.06.00 Original By: ACRM
   V1.1  30.06.00 Modified to make the sequence hash store just offsets
                  in the file rather than the actual sequence data to
                  reduce memory usage
   V1.2  20.07.00 Rewrote all the deletion logic so it actually works!

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <gdbm.h>
#include <signal.h>
#include <sys/types.h>
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"


/************************************************************************/
/* Defines and macros
*/
#define MAX_KEY_LEN        32
#define MAXBUFF           320
#define HUGEBUFF          800

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
char *GetSequence(FILE *fp, char *key);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   15.06.00 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *fp1, *fp2;
   char id1[MAX_KEY_LEN],
        id2[MAX_KEY_LEN],
        *seq1, *seq2;
   
   
   /* Open the sequence file twice for reading                          */
   fp1 = fopen(argv[1], "r");
   fp2 = fopen(argv[1], "r");
   
   /* for each sequence                                                 */
   while((seq1=GetSequence(fp1, id1))!=NULL)
   {
      /* for each other sequence                                        */
      while((seq2=GetSequence(fp2, id2))!=NULL)
      {
         /* if it's not the same sequence                               */
         if(strcmp(id1,id2))
         {
            if(strstr(seq1, seq2))
            {
               printf("%s is a subsequence of %s\n", id2, id1);
            }
         }
         free(seq2);
      }
      rewind(fp2);
      free(seq1);
   }

   fclose(fp1);
   fclose(fp2);
   
   return(0);
}

char *GetSequence(FILE *fp, char *key)
{
   int  ch,
        size = 0;
   char *sequence, *id,
        ptr[HUGEBUFF];
   sequence = NULL;

   /* Get a line until we have one which starts with a >                */
   while((fgets(ptr,HUGEBUFF,fp))!=NULL)
   {
      if(ptr[0] == '>')
         break;
   }

   /* Find the identifier                                         */
   if((id = strchr(ptr, '|'))!=NULL)
   {
      strncpy(key, id+1, MAX_KEY_LEN-1);
      key[MAX_KEY_LEN-1] = '\0';
      
      /* If the original string (ptr) started with PDB we need to
         take the chain name if specified
      */
      if(!strncmp(ptr+1,"pdb",3))
      {
         TERMINATE(key);
         if((id = strchr(key, '|'))!=NULL)
         {
            if(*(id+1))
            {
               *(id+2) = '\0';
            }
            else
            {
               *id = '\0';
            }
         }
      }
      else  /* Something other than PDB, just take the ID         */
      {
         if((id = strchr(key, '|'))!=NULL)
         {
            *id = '\0';
         }
      }
   }
   else
   {
      strncpy(key,ptr+1,MAX_KEY_LEN);
      key[MAX_KEY_LEN-1] = '\0';
   }


   /* Now keep reading to get the sequence                              */
   for(;;)
   {
      if((ch = fgetc(fp)) != EOF)
      {
         ungetc(ch, fp);                 /* Put the character back      */
         if(ch == '>')                   /* Start of next entry         */
         {
            return(sequence);
         }
         else
         {
            if((fgets(ptr,HUGEBUFF,fp))==NULL)
               return(sequence);
            TERMINATE(ptr);
            size += strlen(ptr);
            if(sequence==NULL)
            {
               sequence = (char *)malloc((1+size)*sizeof(char));
               sequence[0] = '\0';
            }
            else
            {
               sequence = (char *)realloc(sequence, (1+size)*sizeof(char));
            }
           
            strcat(sequence, ptr);
         }
      }
      else
      {
         return(sequence);
      }
   }
   
   
}

      
