#!/usr/local/bin/perl
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 1999
#   Author:     Dr. Andrew C. R. Martin
#   Address:    School of Animal and Microbial Sciences,
#               The University of Reading,
#               Whiteknights,
#               P.O. Box 228,
#               Reading RG6 6AJ.
#               England.
#   Phone:      +44 (0)118 987 5123 Extn. 7022
#   Fax:        +44 (0)118 931 0180
#   EMail:      a.c.r.martin@reading.ac.uk
#               andrew@stagleys.demon.co.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
# Read in the superceed list
while(<>)
{
    if(/^INFO:.*superceeds/)
    {
        chomp;
        ($j, $p, $j, $c) = split;
        if(defined($parent{$c}))
        {
            print "$c has more than one parent\n";
        }
        $parent{$c} = $p;          # Store the parent of the child
        push @{$children{$p}}, $c; # Add to array of children for this parent
    }
}

# Build a list of all parents which are not also children
foreach $c (keys %parent)
{
    if(!defined($parent{$parent{$c}}))
    {
        $issuper{$parent{$c}} = 1;
    }
}


# Run through the keys of super-parents and list all children
foreach $super (keys %issuper)
{
    print "$super : ";
    PrintChildren($super);
    print "\n";
}

sub PrintChildren
{
    my($p) = @_;
    my $c;

    foreach $c (@{$children{$p}})
    {
        print "$c ";
        PrintChildren($c);
    }
}
