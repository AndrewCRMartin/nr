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
#include <sys/types.h>
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"


/************************************************************************/
/* Defines and macros
*/
#define DEFAULT_FRAGSIZE   15
#define MAX_KEY_LEN        32
#define MAXBUFF           320
#define HUGEBUFF          800
#define BLOCK_SIZE       4096
#define MODE             0600

#define TOO_MANY_X_FRAC     (REAL)0.25

#define DEFAULT_SEQHASH     "seqhash"
#define DEFAULT_FRAGHASH    "fraghash"
#define DEFAULT_FRAGTABLE   "fragtablehash"
#define DEFAULT_TEMPSEQHASH "seqhash_temp"
#define DEFAULT_DELETEDHASH "deletedhash"
#define DEFAULT_GDBM_DIR    "/tmp"

#define CREATEDATUM(x,y)                                                 \
   (x).dptr = (y);                                                       \
   (x).dsize = (strlen(y)+1)

#define CLEARHASH(ch_hash, ch_stem)                                      \
do                                                                       \
{  pid_t pid;                                                            \
   char  name[MAXBUFF];                                                  \
   gdbm_close(ch_hash);                                                  \
   pid = getpid();                                                       \
   sprintf(name,"%s/%s.%d",gGDBMDir,ch_stem,pid);                        \
   unlink(name);                                                         \
   if((ch_hash = gdbm_open(name, BLOCK_SIZE,                             \
                            GDBM_WRCREAT|GDBM_FAST,                      \
                            MODE, NULL))==NULL)                          \
   {  fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n",         \
              name);                                                     \
      return(FALSE);                                                     \
   }                                                                     \
} while(0)


/************************************************************************/
/* Globals
*/
int       gVerbose = 0;
GDBM_FILE gDBF_seqdata,
          gDBF_seqdata_temp,
          gDBF_fragdata,
          gDBF_fragtable,
          gDBF_deleted;

char      gGDBMDir[MAXBUFF];


/************************************************************************/
/* Prototypes
*/
int CompareSequences(char *seq1, char *id1, char *seq2, char *id2);
BOOL CreateHashes(void);
BOOL ParseCmdLine(int argc, char **argv, char *outfile, BOOL *FirstIsNR, 
                  int *fragSize, int *firstFile, int *rejectSize);
int main(int argc, char **argv);
void CleanUp(void);
BOOL ReadSequences(FILE *in, char *file, int rejectSize);
BOOL HashSequences(int fragSize, BOOL loadOnly);
void StoreSequenceFragment(char *data, 
                           int fragSize,
                           datum gdbm_seq_seqid, BOOL loadOnly);
BOOL PurgeDeletedSequences(void);
void DropSequence(char *seqid);
BOOL DropRedundancies(int fragSize);
void doDropRedundancy(char *seqid, char *sequence, int fragSize);
BOOL NonRedundantise(char *file, BOOL loadOnly, int fragSize, 
                     int rejectSize);
BOOL MergeSequenceHashes(char *mainhash, char *temphash);
void Usage(void);
char *GetSequence(datum content, BOOL full);
void WriteResults(FILE *out);
char *ThisSequenceRedundant(char *data, int fragSize,
                            datum gdbm_seq_seqid);
BOOL TooManyXs(char *seq);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   15.06.00 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   BOOL FirstIsNR  = FALSE;
   int  fragSize   = DEFAULT_FRAGSIZE,
        rejectSize = 2 * DEFAULT_FRAGSIZE,
        firstFile  = 0,
        i;
   char outfile[MAXBUFF],
        *cptr;
   FILE *out       = stdout;

   /* Set GDBM storage directory.
      Initially use the hard-coded default. Replace this with anything
      specified by an environment variable. This may later be over-ridden
      by someting on the command line
   */
   strcpy(gGDBMDir, DEFAULT_GDBM_DIR);
   if((cptr=getenv("NR_TMPDIR")) != NULL)
   {
      strncpy(gGDBMDir, cptr, MAXBUFF);
      gGDBMDir[MAXBUFF-1] = '\0';
   }
   
   
   if(ParseCmdLine(argc, argv, outfile, &FirstIsNR, &fragSize, 
                   &firstFile, &rejectSize))
   {
      if(firstFile && CreateHashes())
      {
         /* Open a different output file if specified                   */
         if(outfile[0])
         {
            if((out=fopen(outfile,"w"))==NULL)
            {
               fprintf(stderr,"nr: (E001) Can't write %s\n", outfile);
               return(1);
            }
         }
         
         /* Step through each input file                                */
         for(i=firstFile; i<argc; i++)
         {
            NonRedundantise(argv[i], 
                            ((i==firstFile)?FirstIsNR:FALSE),
                            fragSize, rejectSize);
         }
         
         /* Write the NR output                                         */
         WriteResults(out);
         if(out!=stdout) fclose(out);
      }
      CleanUp();
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *outfile, 
                     BOOL *FirstIsNR, int *fragSize, int *firstFile,
                     int *rejectSize)
   -----------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *outfile     Output file (or blank string)
            BOOL   *FirstIsNR   Is first file already non-redundant?
            int    *fragSize    Fragment size for hashing
            int    *firstFile   Offset into argv of first input file
            int    *rejectSize  Reject sequences shorter than this
   Returns: BOOL                Success?

   Parse the command line
   
   09.06.00 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *outfile, BOOL *FirstIsNR, 
                  int *fragSize, int *firstFile, int *rejectSize)
{
   argc--;
   argv++;

   outfile[0] = '\0';
   *firstFile=1;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            *FirstIsNR = TRUE;
            (*firstFile)++;
            break;
         case 'o':
            argc--;
            argv++;
            strncpy(outfile,argv[0],MAXBUFF);
            outfile[MAXBUFF-1] = '\0';
            (*firstFile)+=2;
            break;
         case 'd':
            argc--;
            argv++;
            strncpy(gGDBMDir,argv[0],MAXBUFF);
            gGDBMDir[MAXBUFF-1] = '\0';
            (*firstFile)+=2;
            break;
         case 'f':
            argc--;
            argv++;
            sscanf(argv[0],"%d",fragSize);
            (*firstFile)+=2;
            break;
         case 'r':
            argc--;
            argv++;
            sscanf(argv[0],"%d",rejectSize);
            (*firstFile)+=2;
            break;
         case 'v':
            gVerbose++;
            (*firstFile)++;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         if(argc==0) *firstFile = 0;
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   if(argc==0) *firstFile = 0;
   return(TRUE);
}


/************************************************************************/
BOOL CreateHashes(void)
{
   char  name[MAXBUFF];
   pid_t pid;

   /* Open the hashes we use read/write                                 */
   pid = getpid();
      
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_SEQHASH,pid);
   if((gDBF_seqdata = gdbm_open(name, BLOCK_SIZE, 
                                GDBM_WRCREAT|GDBM_FAST, 
                                MODE, NULL))==NULL)
   {
      fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n", name);
      return(FALSE);
   }
   
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_TEMPSEQHASH,pid);
   if((gDBF_seqdata_temp = gdbm_open(name, BLOCK_SIZE, 
                                     GDBM_WRCREAT|GDBM_FAST, 
                                     MODE, NULL))==NULL)
   {
      fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n", name);
      return(FALSE);
   }
   
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_FRAGHASH,pid);
   if((gDBF_fragdata = gdbm_open(name, BLOCK_SIZE, 
                                 GDBM_WRCREAT|GDBM_FAST,
                                 MODE, NULL))==NULL)
   {
      fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n", 
              name);
      return(FALSE);
   }

   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_FRAGTABLE,pid);
   if((gDBF_fragtable = gdbm_open(name, BLOCK_SIZE, 
                                 GDBM_WRCREAT|GDBM_FAST,
                                 MODE, NULL))==NULL)
   {
      fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n", 
              name);
      return(FALSE);
   }

   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_DELETEDHASH,pid);
   if((gDBF_deleted = gdbm_open(name, BLOCK_SIZE, 
                                 GDBM_WRCREAT|GDBM_FAST,
                                 MODE, NULL))==NULL)
   {
      fprintf(stderr,"E00n: Can't open GDBM hash for r/w: %s\n", 
              name);
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>void CleanUp(void)
   ------------------
   Removes the hash files 

   15.06.00 Original   By: ACRM
*/
void CleanUp(void)
{
   char name[MAXBUFF];
   pid_t pid;

   gdbm_close(gDBF_seqdata_temp);
   gdbm_close(gDBF_seqdata);
   gdbm_close(gDBF_fragdata);
   gdbm_close(gDBF_fragtable);
   gdbm_close(gDBF_deleted);

   pid = getpid();
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_SEQHASH,pid);
   unlink(name);
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_TEMPSEQHASH,pid);
   unlink(name);
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_FRAGHASH,pid);
   unlink(name);
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_FRAGTABLE,pid);
   unlink(name);
   sprintf(name,"%s/%s.%d",gGDBMDir,DEFAULT_DELETEDHASH,pid);
   unlink(name);
}


/************************************************************************/
/*>BOOL ReadSequences(FILE *in, char *file, int rejectSize)
   --------------------------------------------------------
   Input:     FILE   *in      FASTA file to read
              char   *file    Filename
   Returns:   BOOL            Success?

   Read a sequence file into a GDBM hash. Checks for duplicate IDs during
   loading.

   The hash contains IDs as keys and the complete FASTA entry as data
   pointers

   15.06.00 Original   By: ACRM
   30.06.00 Modified such that the data in the hash is an offset into the
            file (from ftell()) rather than the actual data
            Added code to reject seqs <= 2*DEFAULT_FRAGSIZE
   12.07.00 rejectSize passed in as parameter instead of 
            2*DEFAULT_FRAGSIZE
*/
BOOL ReadSequences(FILE *in, char *file, int rejectSize)
{
   char      ptr[HUGEBUFF], *id,
             key[MAX_KEY_LEN],
             *sequence;
   datum     gdbm_seq_seqid,
             gdbm_seq_seqdata;
   long      entryStart = (-1),
             thisEntryStart;
   
   if(gVerbose > 1)
   {
      fprintf(stderr,"Reading Sequences...\n");
   }
   
   /* Read through the FASTA input file using a GDBM hash to store the
      sequence keyed by its identifier
   */
   sequence = NULL;
   key[0]   = '\0';
   while((fgets(ptr,HUGEBUFF,in))!=NULL)
   {
      if(*ptr == '>')             /* Start of new entry                 */
      {
         thisEntryStart = ftell(in) - strlen(ptr);
         
         /* If we have a sequence already then store it                 */
         if(entryStart != (-1) && key[0])
         {
            char entryStartString[80];
            char *sptr;
            sptr = strchr(sequence,'\n');
            if(sptr == NULL)
               sptr = sequence;

            if(strlen(sptr) > rejectSize)
            {
               sprintf(entryStartString,"%s %ld", file, entryStart);
               
               CREATEDATUM(gdbm_seq_seqid,   key);
               CREATEDATUM(gdbm_seq_seqdata, entryStartString);
               
               if(gdbm_store(gDBF_seqdata_temp, gdbm_seq_seqid, 
                             gdbm_seq_seqdata, GDBM_INSERT))
               {
                  fprintf(stderr,"Warning (W001): Duplicate ID: %s\n", 
                          key);
               }
            }
            else if(gVerbose)
            {
               TERMINATE(key);
               fprintf(stderr,"INFO: Sequence %s rejected. Only %d \
residues\n",
                       key, strlen(sequence));
            }
            
            free(sequence);
            sequence = NULL;
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
         

         /* Update the pointer to the start of this entry               */
         entryStart = thisEntryStart;
      }
      else
      {
         TERMINATE(ptr);
      }

      /* Build this line into the sequence string                       */
      sequence = strcatalloc(sequence, ptr);
   }
   
   /* If we have a sequence already then store it                       */
   if((entryStart != (-1)) && key[0])
   {
      char entryStartString[80];
      char *sptr;
      sptr = strchr(sequence,'\n');
      if(sptr == NULL)
         sptr = sequence;
            
      if(strlen(sptr) > rejectSize)
      {
         sprintf(entryStartString,"%s %ld", file, entryStart);
         
         CREATEDATUM(gdbm_seq_seqid,     key);
         CREATEDATUM(gdbm_seq_seqdata, entryStartString);
         
         if(gdbm_store(gDBF_seqdata_temp, gdbm_seq_seqid, 
                       gdbm_seq_seqdata, GDBM_INSERT))
         {
            fprintf(stderr,"Warning (W001): Duplicate ID: %s\n", key);
         }
      }
      else if(gVerbose)
      {
         fprintf(stderr,"Sequence %s rejected. Only %d residues\n",
                 key, strlen(sequence));
      }
      free(sequence);
      sequence = NULL;
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL HashSequences(int fragSize, BOOL loadOnly)
   -----------------------------------------------
   Input:     int    fragSize      Fragment size for hashing
              BOOL   loadOnly      Current set of sequences is already
                                   non-redundant, don't do any checking
   Returns:   BOOL                 Success?

   Create a hash of sequence fragments pointing to sequence IDs

   This the has keys are sequence fragments; the data are the sequence IDs

   15.06.00 Original   By: ACRM
*/
BOOL HashSequences(int fragSize, BOOL loadOnly)
{
   char      *data;
   datum     gdbm_seq_seqid,
             gdbm_seq_seqdata;
   
   if(gVerbose > 1)
   {
      fprintf(stderr,"Hashing Sequence Fragments...\n");
   }
   
   gdbm_seq_seqid = gdbm_firstkey(gDBF_seqdata_temp);
   while(gdbm_seq_seqid.dptr)
   {
      gdbm_seq_seqdata = gdbm_fetch(gDBF_seqdata_temp, gdbm_seq_seqid);
      
      if((data = GetSequence(gdbm_seq_seqdata, FALSE))!=NULL)
      {
         if(TooManyXs(data))
         {
            fprintf(stderr,"W00n: Too many Xs in sequence %s\n", 
                    gdbm_seq_seqid.dptr);
         }
         else
         {
            StoreSequenceFragment(data, fragSize, gdbm_seq_seqid, 
                                  loadOnly);
         }
         free(data);
      }
      if(gdbm_seq_seqdata.dptr!=NULL)
      {
         free(gdbm_seq_seqdata.dptr);
      }
      gdbm_seq_seqid=gdbm_nextkey(gDBF_seqdata_temp,gdbm_seq_seqid);
   }

   return(PurgeDeletedSequences());
}


/************************************************************************/
BOOL TooManyXs(char *seq)
{
   if(strchr(seq, 'X'))
   {
      if(((REAL)countchar(seq,'X') / (REAL)strlen(seq)) > TOO_MANY_X_FRAC)
         return(TRUE);
   }
   return(FALSE);
}


/************************************************************************/
/*>void StoreSequenceFragment(char *data, 
                              int fragSize,
                              datum gdbm_seq_seqid, BOOL loadOnly)
   ---------------------------------------------------------------
   Input:     char        *data            A sequence to store
              int         fragSize         Size of fragment
              datum       gdbm_seq_seqid   GDBM datum of sequence ID
              BOOL        loadOnly         Load only, no checking for
                                           sequence match

   Store a sequence fragment into the GDBM fragment hash. First tries
   the N-terminus as the fragment, if this is already stored, then
   slide along the sequence to find a new fragment and repeat until
   we've managed to store the sequence. If we never find a unique
   fragment, generate a warning message.

   15.06.00 Original By: ACRM 
*/
void StoreSequenceFragment(char *data, 
                           int fragSize,
                           datum gdbm_seq_seqid, BOOL loadOnly)
{
   static char *sFragment=NULL;
   int         maxoffset,
               offset;
   BOOL        done     = FALSE,
               gotX     = FALSE,
               gotNoX   = FALSE;
   datum       gdbm_frag_key,
               gdbm_seq_seqdata;


   if(sFragment==NULL)
   {
      /* Allocate memory for fragment storage                           */
      if((sFragment = (char *)malloc((fragSize+1) * sizeof(char)))==NULL)
      {
         fprintf(stderr,"E00n: No memory for fragment storage\n");
         exit(1);
      }
   }
   
   /* Find max possible offset for a fragment                           */
   maxoffset = strlen(data) - fragSize;
   
   /* Keep trying until we've suceeded in inserting this sequence or
      decided that it is redundant
   */
   for(offset=0; offset<maxoffset; offset++)
   {
      strncpy(sFragment, data+offset, fragSize);
      sFragment[fragSize-1] = '\0';

      CREATEDATUM(gdbm_frag_key, sFragment);

      /* Try to store this fragment in the hash                         */
      if(!gdbm_store(gDBF_fragdata, gdbm_frag_key, gdbm_seq_seqid, 
                     GDBM_INSERT))
      {
         /* Stored OK, store the reverse version as well and break out of
            the loop                            
         */
         gdbm_store(gDBF_fragtable, gdbm_seq_seqid, gdbm_frag_key, 
                    GDBM_INSERT);
         done = TRUE;
         break;
      }
   }

   if(!done)
   {
      if(loadOnly)
      {
         fprintf(stderr,"W00n: Can't find unique fragment. Unable to \
store %s (length=%d)\n", gdbm_seq_seqid.dptr, strlen(data));
         if(gVerbose > 2)
         {
            for(offset=0; offset<maxoffset; offset++)
            {
               strncpy(sFragment, data+offset, fragSize);
               sFragment[fragSize-1] = '\0';
               CREATEDATUM(gdbm_frag_key, sFragment);
               gdbm_seq_seqdata = gdbm_fetch(gDBF_fragdata,
                                             gdbm_frag_key);
               fprintf(stderr,"      Hit with: %s\n",
                       gdbm_seq_seqid.dptr);
               if(gdbm_seq_seqdata.dptr != NULL)
               {
                  free(gdbm_seq_seqdata.dptr);
               }
            }
         }
      }
      else
      {
         char *parent;
         
         /* Run through the fragments again to see if any identified
            match is a parent of this sequence
         */
         if((parent = ThisSequenceRedundant(data, fragSize, 
                                            gdbm_seq_seqid))!=NULL)
         {
            if(gVerbose)
            {
               fprintf(stderr, "INFO: %s superceeds %s\n",
                       parent, gdbm_seq_seqid.dptr);
            }
         }
         else
         {
            fprintf(stderr,"W00n: Can't find unique fragment. \
Unable to store %s (length=%d)\n", gdbm_seq_seqid.dptr, strlen(data));
         }
      }
      
      DropSequence(gdbm_seq_seqid.dptr);
   }
}


/************************************************************************/
char *ThisSequenceRedundant(char *data, int fragSize,
                            datum gdbm_seq_seqid)
{
   static char *sFragment=NULL,
               sID[MAX_KEY_LEN];
   int         maxoffset,
               offset,
               seqnum;
   datum       gdbm_frag_key,
               gdbm_stored_key,
               gdbm_stored_seq;
   char        *frag_sequence;


   if(sFragment==NULL)
   {
      /* Allocate memory for fragment storage                           */
      if((sFragment = (char *)malloc((fragSize+1) * sizeof(char)))==NULL)
      {
         fprintf(stderr,"E00n: No memory for fragment storage\n");
         exit(1);
      }
   }
   
   /* Find max possible offset for a fragment                           */
   maxoffset = strlen(data) - fragSize;
   
   /* We know all fragments are already in the fragment hash. Try each
      in turn to see whether the corresponding stored protein is a
      parent
   */
   for(offset=0; offset<maxoffset; offset++)
   {
      strncpy(sFragment, data+offset, fragSize);
      sFragment[fragSize-1] = '\0';

      if(strchr(sFragment,'X'))
      {
         continue;
      }
   
      CREATEDATUM(gdbm_frag_key, sFragment);

      /* Fetch the identifier for this fragment                         */
      gdbm_stored_key = gdbm_fetch(gDBF_fragdata, gdbm_frag_key);
      /* Fetch the sequence pointer for this identifier                 */
      gdbm_stored_seq = gdbm_fetch(gDBF_seqdata_temp, gdbm_stored_key);
      if(gdbm_stored_seq.dptr == NULL)
      {
         gdbm_stored_seq = gdbm_fetch(gDBF_seqdata, gdbm_stored_key);
      }
      if(gdbm_stored_seq.dptr == NULL)
      {
         return(NULL);
      }

      /* Now fetch the complete sequence for this fragment              */
      if((frag_sequence = GetSequence(gdbm_stored_seq, FALSE))!=NULL)
      {
         if((seqnum = CompareSequences(data, gdbm_seq_seqid.dptr, 
                                       frag_sequence, 
                                       gdbm_stored_key.dptr)))
         {
            strcpy(sID, gdbm_stored_key.dptr);
            free(gdbm_stored_key.dptr);
            free(frag_sequence);
            return(sID);
         }
         free(frag_sequence);
      }
      free(gdbm_stored_key.dptr);
   }

   return(NULL);
}


/************************************************************************/
/* We can't actually delete the seqid->sequence references here since
   that would disrupt the loop through the sequences. Instead we remove
   them from the fragment hashes and simply mark this sequence ID as
   having been deleted

   PurgeDeletedSequences() then does the actual work of deleting the
   entries in the seqid->sequence hash.
*/

void DropSequence(char *seqid)
{
   datum  gdbm_seq_key,
          gdbm_fragment,
          gdbm_deleted;

   /* Delete this sequence from the seqid->sequence hashes              */
   CREATEDATUM(gdbm_seq_key, seqid);
   CREATEDATUM(gdbm_deleted, "1");
   gdbm_store(gDBF_deleted, gdbm_seq_key, gdbm_deleted, GDBM_INSERT);
   
   /* Find the associated fragment and delete that                      */
   gdbm_fragment = gdbm_fetch(gDBF_fragtable, gdbm_seq_key);
   gdbm_delete(gDBF_fragdata, gdbm_fragment);
   if(gdbm_fragment.dptr != NULL)
   {
      free(gdbm_fragment.dptr);
   }

   gdbm_delete(gDBF_fragtable, gdbm_seq_key);
}


/************************************************************************/
BOOL PurgeDeletedSequences(void)
{
   datum     gdbm_seq_seqid;

   /* Loop through the keys of the deleted sequence hash                */
   gdbm_seq_seqid = gdbm_firstkey(gDBF_deleted);
   while(gdbm_seq_seqid.dptr)
   {
      gdbm_delete(gDBF_seqdata, gdbm_seq_seqid);
      gdbm_delete(gDBF_seqdata_temp, gdbm_seq_seqid);

      gdbm_seq_seqid=gdbm_nextkey(gDBF_deleted,gdbm_seq_seqid);
   }

   CLEARHASH(gDBF_deleted, DEFAULT_DELETEDHASH);
   return(TRUE);
}


/************************************************************************/
/*>BOOL NonRedundantise(char *file, BOOL loadOnly, int fragSize,
                        int rejectSize)
   -------------------------------------------------------------
   Input:     char   *file    The file to be processed
              BOOL   loadOnly This file already non-redundant, just load
                              it
              int    fragSize Fragment size
              int    rejectSize Reject fragments up to this length
   Returns:   BOOL            Success?

   Non-redundantise the specified file against those already loaded. If
   loadOnly specified, then this file is already redundant: just load it
   for other files to be processed against.

   15.06.00 Original   By: ACRM
*/
BOOL NonRedundantise(char *file, BOOL loadOnly, int fragSize, 
                     int rejectSize)
{
   FILE *in = NULL;
   BOOL retval=TRUE;

   if(gVerbose > 1)
   {
      fprintf(stderr,"\n\nNON-REDUNDANTISING %s\n\n", 
              ((file)?file:"STDIN"));
   }
   
   /* Open the file                                                     */
   if(file)
   {
      if((in=fopen(file, "r"))==NULL)
      {
         fprintf(stderr,"E00n: Can't read %s\n", file);
         return(FALSE);
      }
   }
   else
   {
      return(FALSE);
   }

   /* Read in the sequence data into a GDBM hash                        */
   if(ReadSequences(in, file, rejectSize))
   {
      if(HashSequences(fragSize,loadOnly))
      {
         if(!loadOnly)
         {
            DropRedundancies(fragSize);
         }

         if(!MergeSequenceHashes(DEFAULT_SEQHASH,DEFAULT_TEMPSEQHASH))
         {
            retval = FALSE;
         }
      }
      else
      {
         retval = FALSE;
      }
   }
   else
   {
      fprintf(stderr,"E00n: Failed to read sequences from %s\n",
              file);
      retval = FALSE;
   }
   
   fclose(in);
   return(retval);
}


/************************************************************************/
/*>BOOL DropRedundancies(int fragSize)
   -----------------------------------
   Input:     int    fragSize       Fragment size
   Returns:   BOOL                  Success?

   Make a second pass through the current temporary sequence data 
   marking redundant sequences for deletion either in the temporary
   hash or in the main hash

   This routine opens the hashes and loops through the temporary file.
   Calls doDropRedundancy() to do the actual work of checking for
   redundancy and marking for deletion

   15.06.00 Original   By: ACRM
*/
BOOL DropRedundancies(int fragSize)
{
   datum     gdbm_seq_seqid,
             gdbm_seq_seqdata,
             gdbm_deleted;
   char      *data;
   
   
   if(gVerbose > 1)
   {
      fprintf(stderr,"Dropping Redundancies...\n");
   }
   
   /* Loop through the keys of the temporary sequence hash              */
   gdbm_seq_seqid = gdbm_firstkey(gDBF_seqdata_temp);
   while(gdbm_seq_seqid.dptr)
   {
      gdbm_deleted = gdbm_fetch(gDBF_deleted, gdbm_seq_seqid);
      if(gdbm_deleted.dptr == NULL)
      {
         /* This sequence hasn't already been marked as deleted         */
         gdbm_seq_seqdata = gdbm_fetch(gDBF_seqdata_temp, gdbm_seq_seqid);
         
         if((data = GetSequence(gdbm_seq_seqdata, FALSE))!=NULL)
         {
            doDropRedundancy(gdbm_seq_seqid.dptr, data, fragSize);
            free(data);
         }
         
         if(gdbm_seq_seqdata.dptr)
         {
            free(gdbm_seq_seqdata.dptr);
         }
      }
      else
      {
         free(gdbm_deleted.dptr);
      }

      gdbm_seq_seqid=gdbm_nextkey(gDBF_seqdata_temp,gdbm_seq_seqid);
   }

   return(PurgeDeletedSequences());
}


/************************************************************************/
/*>void doDropRedundancy(char *seqid, char *sequence, int fragSize)
   -----------------------------------------------------------------------
   Input:     char       *seqid        Sequence identifier to test
              char       *sequence     Sequence to test
              int        fragSize      Fragment size
              
   Does the actual checking of a sequence against the fragment hash and
   looking for redundancy then marking a redundant sequence for deletion

   15.06.00 Original   By: ACRM
*/
void doDropRedundancy(char *seqid, char *sequence, int fragSize)
{
   static char *sFragment=NULL;
   int         maxoffset,
               offset,
               fragnum;
   datum       gdbm_frag_seqid,
               gdbm_frag_key,
               gdbm_seq_data;
   char        *stored_data;
   

   if(sFragment==NULL)
   {
      /* Allocate memory for fragment storage                           */
      if((sFragment = (char *)malloc((fragSize+1) * sizeof(char)))==NULL)
      {
         fprintf(stderr,"E00n: No memory for fragment storage\n");
         exit(1);
      }
   }
   
   /* Find max possible offset for a fragment                           */
   maxoffset = strlen(sequence) - fragSize;

   for(offset=0; offset<maxoffset; offset++)
   {
      strncpy(sFragment, sequence+offset, fragSize);
      sFragment[fragSize-1] = '\0';

      /* Try to fetch a sequence ID for this fragment                   */
      CREATEDATUM(gdbm_frag_key,sFragment);
      gdbm_frag_seqid = gdbm_fetch(gDBF_fragdata, gdbm_frag_key);
      
      /* If we found one                                                */
      if(gdbm_frag_seqid.dptr)
      {
         /* If it isn't the self match                                  */
         if(strcmp(seqid, gdbm_frag_seqid.dptr))
         {
            /* Grab the found sequence
               First try the main sequence hash and then the temp hash
            */
            gdbm_seq_data = gdbm_fetch(gDBF_seqdata, gdbm_frag_seqid);
            if(gdbm_seq_data.dptr == NULL)
            {
               gdbm_seq_data = gdbm_fetch(gDBF_seqdata_temp, 
                                          gdbm_frag_seqid);
            }

            /* If we found it OK                                        */
            if(gdbm_seq_data.dptr)
            {
               /* Compare the sequences                                 */
               if((stored_data = GetSequence(gdbm_seq_data, FALSE))!=NULL)
               {
                  if((fragnum=CompareSequences(sequence, 
                                               seqid,
                                               stored_data,
                                               gdbm_frag_seqid.dptr)))
                  {
                     /* Sequences are the same, if the first one is 
                        longer then replace it in the fragment hash
                     */
                     if(gVerbose)
                     {
                        if(fragnum==1)
                        {
                           fprintf(stderr,"INFO: %s superceeds %s\n",
                                   seqid, 
                                   gdbm_frag_seqid.dptr);
                        }
                        else
                        {
                           fprintf(stderr,"INFO: %s superceeds %s\n",
                                   gdbm_frag_seqid.dptr,
                                   seqid);
                        }
                     }
                  
                     if(fragnum==1)
                     {
                        DropSequence(gdbm_frag_seqid.dptr);
                     }
                     else
                     {
                        DropSequence(seqid);
                     }

                     /* If our probe sequence is declared redundant, then
                        we have finished
                     */
                     if(fragnum == 2)
                     {
                        free(stored_data);
                        break;  /* Out of for(offset...)                */
                     }
                  }  /* if(sequences match)                             */
                  free(stored_data);
               }  /* if(extracted the actual sequence)                  */
               free(gdbm_seq_data.dptr);
            }  /* if(found seq data for this id in the hashes)          */
         }  /* if(not self-match)                                       */
         free(gdbm_frag_seqid.dptr);
      }  /* if(fragment found)                                          */
   }  /* for(offset...)                                                 */
}


/************************************************************************/
/*>BOOL MergeSequenceHashes(char *mainhash, char *temphash)
   --------------------------------------------------------
   Input:     char     *mainhash   Main sequence hash
              char     *temphash   Temporary sequence hash
   Returns:   BOOL                 Success?

   Merge the data from the temporary sequence has (temphash) into the
   main sequence hash (mainhash) and then delete temphash

   15.06.00 Original   By: ACRM
*/
BOOL MergeSequenceHashes(char *mainhash, char *temphash)
{
   datum     key,
             content;

   
   if(gVerbose > 1)
   {
      fprintf(stderr,"Merging Sequence Hashes...\n");
   }
   
   /* Copy everything across from the temp seq hash to the main one     */
   key = gdbm_firstkey(gDBF_seqdata_temp);
   while(key.dptr)
   {
      content = gdbm_fetch(gDBF_seqdata_temp, key);
      
      if(gdbm_store(gDBF_seqdata, key, content, GDBM_INSERT))
      {
         fprintf(stderr,"Warning (W001): Duplicate ID: %s\n", key.dptr);
      }
      if(content.dptr)
      {
         free(content.dptr);
      }
      
      key=gdbm_nextkey(gDBF_seqdata_temp,key);
   }

   CLEARHASH(gDBF_seqdata_temp, DEFAULT_TEMPSEQHASH);
   return(TRUE);
}


/************************************************************************/
/*>char *GetSequence(datum content, BOOL full)
   -------------------------------------------
   Input:     datum   content   GDBM datum structure containing the
                                filename and fseek() pointer
              BOOL    full      Get the header as well as the sequence
   Returns:   char    *         Pointer to sequence data

   The sequence hash contains an fseek() pointer into the sequence data
   file. This routine reads the actual sequence data into a malloc'd
   block of memory.

   15.06.00 Original   By: ACRM
   30.06.00 'content' is now a filename and fseek() pointer into an
            actual sequence file
   11.07.00 Fixed memory leak
*/
char *GetSequence(datum content, BOOL full)
{
   char   *data = NULL,
          ptr[HUGEBUFF],
          filename[MAXBUFF];
   long   offset;
   static FILE *fp = NULL;
   static char lastFilename[MAXBUFF];

   if(content.dptr==NULL)
      return(NULL);

   if(fp==NULL)
      lastFilename[0] = '\0';

   sscanf(content.dptr,"%s %ld", filename, &offset);

   /* If the filename has changed, open the new file                    */
   if(strcmp(filename,lastFilename))
   {
      if(fp!=NULL)        /* Close the old file                         */
         fclose(fp);
      
      fp=fopen(filename,"r");
      strcpy(lastFilename, filename);
   }
   
   if(fp!=NULL)
   {
      if(fseek(fp, offset, SEEK_SET) == (-1))
         return(NULL);
   
      if(!full)
      {
         /* Throw away the first line                                   */
         fgets(ptr, HUGEBUFF, fp);
      }
      
      while(fgets(ptr, HUGEBUFF, fp)!=NULL)
      {
         if((*ptr == '>') &&        /* Start of new entry. Jump out     */
            (data != NULL)) 
         {
            break;
         }
         else
         {
            if(!full)
            {
               TERMINATE(ptr);
            }
            data = strcatalloc(data, ptr);
         }
      }
   }
   
   return(data);
}


/************************************************************************/
/*>void WriteResults(FILE *out)
   ----------------------------
   Input:     FILE   *out      Output file pointer

   Write the non-redundant sequences to the output file

   15.06.00 Original   By: ACRM
   30.06.00 Modified to use GetSequence()
*/
void WriteResults(FILE *out)
{
   char      *seq;
   datum     gdbm_seq_seqid,
             gdbm_seq_seqdata;
   
   if(gVerbose > 1)
   {
      fprintf(stderr,"Writing Results...\n");
   }
   
   gdbm_seq_seqid = gdbm_firstkey(gDBF_seqdata);
   while(gdbm_seq_seqid.dptr)
   {
      gdbm_seq_seqdata = gdbm_fetch(gDBF_seqdata, gdbm_seq_seqid);
      seq = GetSequence(gdbm_seq_seqdata, TRUE);
      if(gdbm_seq_seqdata.dptr)
      {
         free(gdbm_seq_seqdata.dptr);
      }
      fprintf(out,"%s", seq);
      free(seq);
      gdbm_seq_seqid=gdbm_nextkey(gDBF_seqdata,gdbm_seq_seqid);
   }
}


/************************************************************************/
/*>int CompareSequences(char *seq1, char *id1, char *seq2, char *id2)
   ------------------------------------------------------------------
   Input:     char   *seq1   First sequence
              char   *id1    First sequence ID
              char   *seq2   Second sequence
              char   *id2    Second sequence ID
   Returns:   int            0: Sequences are different
                             1: First sequence is longer
                             2: Second sequence is longer

   Compare two sequences, returning 0 if they differ, 1 if the first 
   sequence is a superset of the second or 2 if the second is a superset
   of the first.

   If the sequences are identical, then the identifiers are compared and
   the sequence is chosen with the alphabetically higher identifier. This
   is needed as otherwise it is possible that given 3 sequences of the
   same length, A replaces B; B replaces C; and C replaces A such that
   all 3 sequences have been deleted (unlikely, but who knows). By
   selecting the alphabetically higher ID, B replaces A, C replaces B
   and C replaces A so we end up with C.

   15.06.00 Original   By: ACRM
   14.07.00 Added id parameters
*/
int CompareSequences(char *seq1, char *id1, char *seq2, char *id2)
{
   int len1, 
       len2;
   
   len1 = strlen(seq1);
   len2 = strlen(seq2);
   
   if(len2 < len1)               /* Seq2 is shorter                     */
   {
      if(strstr(seq1, seq2))
      {
         return(1);
      }
   }
   else if(len1 < len2)          /* Seq1 is shorter                     */
   {
      if(strstr(seq2, seq1))
      {
         return(2);
      }
   }
   else                          /* Same length                         */
   {
      if(!strcmp(seq1, seq2))    /* Sequences are identical             */
      {
         /* Compare the identifiers and return the alphabetically higher
            one
         */
         return((strcmp(id1,id2) > 0)?1:2);
      }
   }
   
   
   return(0);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Write a usage message

   15.06.00 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nnr V1.0 (c) 2000 Dr. Andrew C.R. Martin, University \
of Reading.\n");

   fprintf(stderr,"\nUsage: nr [-v] [-o out.faa] [-n] [-f fragsize] \
[-r size] [-d tmpdir] file1.faa [file2.faa ...]\n");
   fprintf(stderr,"       -v  Verbose mode - Report superceeded \
sequences\n");
   fprintf(stderr,"       -o  Specify output file (stdout if not \
specified)\n");
   fprintf(stderr,"       -n  First sequence file is already \
non-redundant\n");
   fprintf(stderr,"       -f  Specify fragment size (default: %d)\n",
DEFAULT_FRAGSIZE);
   fprintf(stderr,"       -r  Reject sequences up to this length \
(default: %d)\n", 2*DEFAULT_FRAGSIZE);
   fprintf(stderr,"       -d  Specify temporary directory \
(default: %s)\n", DEFAULT_GDBM_DIR);
   
   fprintf(stderr,"\nTwo-pass sequence non-reduntantising program with \
minimal memory usage.\n");
   fprintf(stderr,"Selects maximal length versions of any set of \
identical sequences. Thus\n");
   fprintf(stderr,"for the sequences ABCDEF and XABCDEFY only the latter \
sequence will be\n");
   fprintf(stderr,"retained. In the overlapping region the sequences \
must be identical.\n");

   fprintf(stderr,"The hard-coded temporary directory %s is used for \
storing the hash\n", DEFAULT_GDBM_DIR);
   fprintf(stderr,"files. Since these files are large, you probably \
don't want to use /tmp.\n");
   fprintf(stderr,"This may be overridden using the NR_TMPDIR \
environment variable or using \n");
   fprintf(stderr,"the -d switch on the command line. The command line \
switch takes\n");
   fprintf(stderr,"precedence over the environment variable.\n\n");
}

