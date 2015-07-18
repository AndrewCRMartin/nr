

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

      if(strchr(sFragment,'X'))
      {
         gotX = TRUE;
         continue;
      }
      else
      {
         gotNoX = TRUE;
      }
   
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
      if(gotX && !gotNoX)
      {
         fprintf(stderr,"W00n: No fragment found without X characters. \
Unable to store %s\n", gdbm_seq_seqid.dptr);
      }
      else
      {
         if(loadOnly)
         {
            fprintf(stderr,"W00n: Can't find unique fragment. Unable to \
store %s (length=%d)\n", gdbm_seq_seqid.dptr, strlen(data));
            if(gVerbose > 1)
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
            /* Run through the fragments again to see if any identified
               match is a parent of this sequence
            */
         }
      }
      
      DropSequence(gdbm_seq_seqid.dptr);
   }
}
