# DBGDE
Using pre-made software trying to implement Mutual software algorithm to get differential expression. 


We have:
- One implementation trying to get GFA2 working at: https://github.com/Omig12/Mutual-pepino
- The one I just did is GFA1 with GFApy working and beign read by Bandage. Used segment for both fragment and the segment, the Segment appears first follow by the kmer subsegments.
```
 Format of GFA1 Segment is: S id_of_original_segment contig.
```
```
 Format of GFA1 subSegment is: S  id_of_original_segment:subsegmentid:(A:count,B:count) kmer
```
Where id_of_original_segment is the same as first segment, subsegmentid is this kmer id, a tuple with the count of all the times that kmer appears in the organism.

An example to run the dbg.py is:

```
python dbg.py -k <kmer_size> -A <Organism A file(s)> -B <Organism B file(s)>
```
