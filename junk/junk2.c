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
   the N-terminus as the fragment, if this is already stored compare
   the sequences. If they match, then throw away either the new sequence
   or the old one as appropriate. If they don't match, then slide along
   the sequence to find a new fragment and repeat until we've managed to
   store the sequence. If we never find a unique fragment, generate a
   warning message.

   15.06.00 Original   By: ACRM
*/
void StoreSequenceFragment(char *data, 
                           int fragSize,
                           datum gdbm_seq_seqid, BOOL loadOnly)
{
   static char *sFragment=NULL;
   int         maxoffset,
               offset,
               fragnum;
   BOOL        done     = FALSE,
               compared = FALSE,
               gotX     = FALSE,
               gotNoX   = FALSE;
   datum       gdbm_frag_key,
               gdbm_frag_seqid,
               gdbm_seq_seqdata;
   char        *stored_data;

   if(loadOnly) compared = TRUE;
   
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
      if(gdbm_store(gDBF_fragdata, gdbm_frag_key, gdbm_seq_seqid, 
                    GDBM_INSERT))
      {
         /* Store failed - key already in the database.                 

            If we haven't already compared the sequences, then see if
            the new sequence is actually redundant (i.e. they are the
            same sequence with one longer than the other). If it is
            redundant then we either forget about this sequence or 
            replace the entry in the has with this one if it is longer
            than the one we already have stored.

            If the compared flag is set, then we already know that the
            sequence is unique (i.e. this data set is already 
            non-redundant) so we don't bother checking against the other
            sequence. Just let the for() loop continue to find a 
            distinct fragment.
         */
         if(!compared)
         {
            /* Get the sequence ID used previously with this sequence
               fragment as a key 
            */
            gdbm_frag_seqid = gdbm_fetch(gDBF_fragdata, gdbm_frag_key);

            /* Find the complete sequence for the key already stored    */
            gdbm_seq_seqdata = gdbm_fetch(gDBF_seqdata, gdbm_frag_seqid);

            /* Check to see if this one has previously been deleted     */
            if(gdbm_seq_seqdata.dptr == NULL)
            {
               if(gVerbose > 2)
               {
                  fprintf(stderr,"INFO: Fragment %s was already present \
for deleted entry %s\n", gdbm_frag_key.dptr, gdbm_frag_seqid.dptr);
                  fprintf(stderr,"      Replaced by %s\n", 
                           gdbm_seq_seqid.dptr);
               }
               
               /* Yes it had been deleted, so put this sequence in 
                  instead
               */
               gdbm_store(gDBF_fragdata, gdbm_frag_key, gdbm_seq_seqid,
                          GDBM_REPLACE);
               done = TRUE;
               break;
            }

            /* It wasn't previously deleted, get the whole sequence     */
            if((stored_data = GetSequence(gdbm_seq_seqdata, FALSE))!=NULL)
            {
               if((fragnum=CompareSequences(data, 
                                            gdbm_seq_seqid.dptr,
                                            stored_data, 
                                            gdbm_seq_seqdata.dptr)))
               {
                  /* Sequences are the same, if the first one is longer 
                     then replace it in the fragment hash
                  */
                  if(strcmp(gdbm_seq_seqid.dptr, gdbm_frag_seqid.dptr))
                  {
                     if(fragnum==1)
                     {
                        /* New sequence is longer, 
                           drop the old sequence                        
                        */
                        if(gVerbose)
                        {
                           fprintf(stderr,"INFO: %s superceeds %s\n",
                                   gdbm_seq_seqid.dptr, 
                                   gdbm_frag_seqid.dptr);
                        }
                        DropSequence(gdbm_frag_seqid.dptr);

                        /* Insert the new sequence                      */
                        gdbm_store(gDBF_fragdata, gdbm_frag_key, 
                                   gdbm_seq_seqid, GDBM_REPLACE);
                        gdbm_store(gDBF_fragtable, gdbm_seq_seqid, 
                                   gdbm_frag_key,  GDBM_INSERT);
                     }
                     else
                     {
                        /* Old sequence is longer
                           Drop the new sequence                        
                        */
                        if(gVerbose)
                        {
                           fprintf(stderr,"INFO: %s already present as \
%s\n",
                                   gdbm_seq_seqid.dptr, 
                                   gdbm_frag_seqid.dptr);
                        }
                        DropSequence(gdbm_seq_seqid.dptr);
                     }
                  }
                  
                  done = TRUE;
                  break;
               }

               free(stored_data);
               
               /* Carry on through the loop trying to find a unique 
                  fragment.
               */
            }
         }
      }
      else         /* Stored OK, break out of the loop                  */
      {
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
         fprintf(stderr,"W00n: Can't find unique fragment. Unable to \
store %s (length=%d)\n", gdbm_seq_seqid.dptr, strlen(data));
         if(gVerbose > 1)
         {
            for(offset=0; offset<maxoffset; offset++)
            {
               strncpy(sFragment, data+offset, fragSize);
               sFragment[fragSize-1] = '\0';
               CREATEDATUM(gdbm_frag_key, sFragment);
               gdbm_seq_seqdata = gdbm_fetch(gDBF_fragdata,gdbm_frag_key);
               fprintf(stderr,"      Hit with: %s\n",gdbm_seq_seqid.dptr);
            }
         }
      }
      
      DropSequence(gdbm_seq_seqid.dptr);
   }
   
}


