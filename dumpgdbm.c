/*************************************************************************

   Program:    dumpgdbm
   File:       dumpgdbm.c
   
   Version:    V1.0
   Date:       15.06.00
   Function:   Dump a GDBM database
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 1999
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gdbm.h>
#include <sys/types.h>

/************************************************************************/
/* Defines and macros
*/

#define CREATEDATUM(x,y) \
   (x).dptr = (y);       \
   (x).dsize = strlen(y)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);

/************************************************************************/
int main(int argc, char **argv)
{
   GDBM_FILE dbf;
   int       block_size=4096;
   datum     gdbm_key,
             gdbm_data;

   if((argc==1) || !strncmp(argv[1], "-h", 2))
   {
      Usage();
      return(0);
   }

   if((dbf = gdbm_open(argv[1], block_size, GDBM_READER, 0, NULL))==NULL)
   {
      fprintf(stderr,"Can't open GDBM hash '%s' for reading\n", argv[1]);
      return(1);
   }

   gdbm_key = gdbm_firstkey(dbf);
   while(gdbm_key.dptr)
   {
      gdbm_data = gdbm_fetch(dbf, gdbm_key);
      printf("%s : %s\n", gdbm_key.dptr, gdbm_data.dptr);
      gdbm_key=gdbm_nextkey(dbf,gdbm_key);
   }
   gdbm_close(dbf);
   
   return(0);
}

/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"\ndumpgdbm V1.0 (c) 2000\n\n");
   fprintf(stderr,"Usage: dumpgdbm gdbm_data_file\n");

   fprintf(stderr,"Dumps the contents of a GDBM hash file as key : \
data pairs\n\n");
}
