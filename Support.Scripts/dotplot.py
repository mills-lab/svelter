#!/usr/bin/env python

####
# Modified From: https://github.com/sherbondy/Comparative-Genomics/blob/master/dotplot.py
# 6.047/6.878 - Problem Set 1 - string hashing/dotplots
#
# INSTRUCTIONS FOR USE:
# call program as follows:
#  ./ps1-dotplot.py <FASTA 1> <FASTA 2> <PLOTFILE>
#     e.g. ./ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
#
# Make sure the ps1-dotplot.py is marked as executable:
#     chmod +x ps1-dotplot.py
# or in windows with:
#     python ps1-dotplot.py human-hoxa-region.fa mouse-hoxa-region.fa dotplot.jpg
# once you have put python in your path
#
#
# GNUPLOT
# Gnuplot is used to generate plots for this program.  It is a common plotting
# program installed on most unix systems.  To use gnuplot on athena do the
# following:
#
# athena% add gnu
#
# To test that it works do this:
#
# athena% gnuplot
# gnuplot> plot cos(x)
#
# you should then see a cosine plot appear.
#
# Note: plotting.py and util.py must be in the same directory as this script.
# These files contain the code for generating plots.  You should not
# worry about understanding any of the code contained within these files.
# Much of it is copied from another project and is unrelated.
#



import sys, random
import plotting


def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)


def quality(hits):
    """determines the quality of a list of hits"""

    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    goodhits = []

    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2

        if lower < hit[1] < upper:
            goodhits.append(hit)

    return goodhits


def makeDotplot(filename, hits, lenSeq1, lenSeq2):
    """generate a dotplot from a list of hits
       filename may end in the following file extensions:
         *.ps, *.png, *.jpg
    """
    x, y = zip(* hits)

    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    hits2 = quality(hits)
    #print "%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits)))

    # create plot
    p = plotting.Gnuplot()
    p.enableOutput(False)
    p.plot(x, y, xlab="sequence 2", ylab="sequence 1")
    p.plotfunc(lambda x: slope1 * x + offset1, 0, 1e6, 1e5)
    p.plotfunc(lambda x: slope2 * x + offset2, 0, 1e6, 1e5)

    # set plot labels
    p.set(xmin=0, xmax=lenSeq2, ymin=0, ymax=lenSeq1)
    p.set(main="dotplot (%d hits, %.5f%% hits on diagonal)" %
          (len(hits), 100 * len(hits2) / float(len(hits))))
    p.enableOutput(True)

    # output plot
    p.save(filename)

    return p

invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'}

def subkeys(key, nth_base, inversions):
    subkeys = []
    keylen = len(key)
    	
    # speed tip from:
    # http://wiki.python.org/moin/PythonSpeed/PerformanceTips#String_Concatenation
    if nth_base == 1:
        subkeys = [key]

    elif nth_base != 0:
        for k in range(nth_base):
            substr_list = [key[j] for j in range(keylen) if (j % nth_base == k)]
            subkeys.append("".join(substr_list))

    else:
        # nth_base = 0 is a special case for third base mismatches
        # for every codon, only include the first 2 bases in the hash
        subkeys = ["".join([key[i] for i in range(len(key)) if i % 3 != 2])]

    if inversions:
        for i in range(len(subkeys)):
            subkeys.append("".join([invert_base[c] for c in reversed(subkeys[i])]))

    return subkeys


def kmerhits(seq1, seq2, kmerlen, nth_base=1, inversions=False):
    # hash table for finding hits
    lookup = {}

    # store sequence hashes in hash table
    #print "hashing seq1..."
    seq1len = len(seq1)
    for i in xrange(seq1len - kmerlen + 1):
        key = seq1[i:i+kmerlen]

        for subkey in subkeys(key, nth_base, inversions):
            lookup.setdefault(subkey, []).append(i)

    # match every nth base by 

    # look up hashes in hash table
    #print "hashing seq2..."
    hits = []
    for i in xrange(len(seq2) - kmerlen + 1):
        key = seq2[i:i+kmerlen]

        # only need to specify inversions for one seq
        for subkey in subkeys(key, nth_base, False):
            subhits = lookup.get(subkey, [])
            if subhits != []:
                # store hits to hits list
                for hit in subhits:
                    hits.append((i, hit))
                # break out of loop to avoid doubly counting
                # exact matches
                break

    return hits


def main():
    # parse command-line arguments
    if len(sys.argv) < 4:
        print "you must call program as:  "
        print "   python ps1-dotplot.py <KMER LEN> <FASTA 1> <FASTA 2> <PLOT FILE>"
        print "   PLOT FILE may be *.ps, *.png, *.jpg"
        sys.exit(1)

    kmerlen = int(sys.argv[1])
    file1 = sys.argv[2]
    file2 = sys.argv[3]
    plotfile = sys.argv[4]

    # read sequences
    #print "reading sequences"
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    # match every nth base. Or set to 0 to allow third base mismatches
    nth_base = 1
    inversions = True

    hits = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
    fo=open(plotfile+'.txt','w')
    for x in hits:
	print >>fo,'\t'.join([str(i) for i in x])
    fo.close()
    #
    # hits should be a list of tuples
    # [(index1_in_seq2, index1_in_seq1),
    #  (index2_in_seq2, index2_in_seq1),
    #  ...]
    #

    #print "%d hits found" % len(hits)
    #print "making plot..."
    p = makeDotplot(plotfile, hits, len(seq1), len(seq2))


main()
