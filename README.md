nr V1.2
=======

(c) 2000 University of Reading, Dr. Andrew C.R. Martin
------------------------------------------------------


`nr` is a program to create a non-redundant sequence database from one
or more FASTA input file(s). The program uses an algorithm designed
for maximal speed with minimal memory usage. 

It selects maximal length versions of any set of identical sequences.
Thus for the sequences ABCDEF and XABCDEFY only the latter sequence
will be retained. In the overlapping region the sequences must be
identical. If two sequences are 100% identical (in sequence and in
length) then the one with the alphabetically higher identifier is
chosen. 

Optionally the first file may already be non-redundant in which case
new non-redundant sequences may be added to it.

```
Usage: nr [-v] [-o out.faa] [-n] [-f fragsize] [-r size] [-d tmpdir]
          file1.faa [file2.faa ...]
          -v  Verbose mode. More information supplied depending on
              the number of times -v appears :
              1 - Report superceeded sequences
              2 - Trace progress through program
              3 - Detail
          -o  Specify output file (stdout if not specified)
          -n  First sequence file is already non-redundant
          -f  Specify fragment size for hashing (default: 15)
          -r  Reject sequences up to this length (default: 30)
          -d  Specify temporary directory (default: /tmp)
              This is used for storing the hash files. Since these 
              files are large, you probably don't want to use /tmp.
              The default may also be overridden using the NR_TMPDIR 
              environment variable. The command line switch takes
              precedence over the environment variable.
```

findequiv.pl
------------

This is a small Perl script to analyse the log file produced by `nr`
(when run with `-v`) which generates a list of the top-level parents
(i.e. those which appear in the final output from `nr`) and all their
descendents.

```
Usage:         findequiv.pl nr.log >parents.lis
        --or--
               cat nr*.log | findequiv.pl >parents.lis
```

Warning and Error Messages
--------------------------

```
W001: Duplicate ID
      The identifier is already used and stored in the sequence hashes

W002: Too many Xs in sequence
      We only allow up to 25% of the sequence length to be the letter
      X. Any more than this is probably meaningless for further
      sequence based work and the chances of finding a unique fragment
      to store in the hashes is reduced.

W003: Can't find unique fragment
      No unique fragment can be found for this sequence so it is
      deleted anyway. Generally this is an indication of the sequence
      being redundant (and there also being lots of other redundant
      sequences in the same group), but in odd and very rare 
      circumstances it can happen that a distinct sequence is lost
      (see stage 2 of the algorithm for an explanation).

E001: Can't write file
      Can't open a file for writing

E002: Can't open GDBM hash for r/w
      Can't open a hash for read/write

E003: No memory for fragment storage
      Out of memory!

E004: Can't read file
      Can't open a file for reading

E005: Failed to read sequences from file
      Error occured while reading the sequence data
```


Algorithm
---------

The algorithm used is as follows:

### 1. Read the data

Read the sequence file creating a hash which gives indexes into the
file keyed by identifier. This allows us quickly to obtain the
sequence for a given identifier. This is loaded as a "temporary"
working hash which is later merged with the main sequence hash.

### 2. Hash the sequences 

The N-terminal fragment (default 15aa) from each sequence is stored in
a hash pointing to its identifier and vice versa, the identifier
pointing to the fragment. If the fragment is already stored, then we
slide the window along to look at the next fragment and so on until we
find a unique fragment for this sequence. If no unique fragment can be
found a warning is given and the sequence is rejected. In most cases
this will be an indication of many similar redundant sequences. It is
however possible that in rather bizarre cases a unique sequence could
be lost. For example (assuming a 3res fragment) we have stored:
```
       ABCXXXXXXX : ABC
       BCDXXXXXXX : BCD
       CDEXXXXXXX : CDE
       DEPXXXXXXX : DEP
       EPQXXXXXXX : EPQ
       PQRXXXXXXX : PQR
       QRTXXXXXXX : QRS
       RSTXXXXXXX : RST
```
and wish to store:
```
       ABCDEPQRST
```
All the fragment keys have appeared before.

With a reasonably large fragment size this is very unlikely and
increasing the fragment size at run time should solve the problem. 

### 3. Drop Redundancies

This stage is skipped if this is the first file and has been flagged
as already non-redundant

Runs through each new sequence in turn (i.e. each key in the temporary
sequence hash). Slides a fragment window along the sequence and looks
to see if it stored in the fragment hash (from stage 2). If a match is
found and it is not a self-match, then we compare the sequences. If
one sequence is a "child" of the other, then it is dropped and a
superceed message is issued.

Dropping a sequence doesn't involve actually removing it from the hash
since this would disturb the loop through the keys. Instead we simply
mark it as deleted and at the end of this stage, all marked sequences
are physically removed from all the hashes.

### 4. Merge the hashes

The temporary working sequence hash is then merged into the main
sequence hash and the temporary hash is deleted.

### 5. Repeat

Stages 1-4 are repeated for any other sequence files specified on the
command line

### 6. Write results

All remaining sequences are written to the output file.



Pseudocode
----------

```
// Pass 1: Create a hash of fragments from the sequences
foreach sequence (ID)
{
   select a fragment (F) from the N-terminus working along till a
          unique one is found;
   store F in a hash (hash on F) with ID;
}

// Pass 2: Test every sequence against the hash
foreach sequence (ID)
{
   foreach overlapping fragment (Fo)
   {
      if(Fo is in the hash)
      {
         IDH = hashedsequenceID;
         if(ID != IDH)
         {
            if(sequence(IDH) is a subset of sequence(ID))
            {
               DeleteFromHash(IDH);
            }
         }
      }
   }
}

// Pass 3: Write the NR data
foreach hashedsequence (ID)
{
   WriteData(ID);
}
```

Notes
-----

Performance improvements might be gained from

### 1. Replacing the hashing routines

Hashing is currently done for ease using GDBM hashes. Hand-coded
hashing is likely to give better performance and would remove the need
for disk storage of the hashes. Example code for hashing is at:
http://www.niksula.cs.hut.fi/~tik76122/dsaa_c2e/files.html
It would also be possible to implement non-unique keys such that
sequences for which no unique fragment can be found could be
processed. 

### 2. Delete-as-we-go

The list of sequence IDs for processing could be generated as a linked
list rather than stepping through all the hash keys. The routine to
drop sequences from the hash could then do so properly rather than
making a list of sequences to drop.

### 3. Sequence storage

The actual sequence data could be stored in memory in the sequence
hashes rather than storing a pointer into the file. This would speed
things considerably (CPU time for processing ~330k sequences is around
13 minutes; elapsed time is ~44 minutes) since much of the elapsed
time is spent doing disk I/O. However this would be at the expense of
much increased memory usage. The current implementation required only
~15M RAM to process ~330k sequences.

### 4. Partial mismatches

With the current method it is not possible to reject partial
mismatches.  While the `CompareSequences()` routine could be modified to
allow a partial mismatch (say 2 differences allowed in the sequences),
this depends on the fragment which is hashed being identical between
the two sequences (in order that they ever get compared).

The only way around that is to hash all non-overlapping fragments from
a sequence instead of just one fragment.


