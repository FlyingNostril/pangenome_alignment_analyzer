The pangenome alignment analyzer accomplishes a few main goals:
1) Take input regions or genes and produces alignments for each genome included, provided a BLAST match can be found
2) Translate the matching gene sequence from each genome and aligns those as well. 
3) Classify the translated protein alignments and genomic alignments by the number and location of the variants found in the alignment
4) produce Fasta DNA/AA alignments that can be used with alignment viewing software if a more detailed look is desired.


Additional goals to work into the project if I can find the time:
Filter out or retry poor BLAST results and alignments
generate a phylogenetic tree based on the alignment results


Each script in the pipeline has a help associated with it, all you have to do is run:
python3 command.py -h

The scripts do have several dependencies you need to install on your own:
ncbi blast+ 
--needed to run BLAST as well as make the blast databases
muscle version 5 or later
--it uses SUPER and has finer control over the threads and memory used

The only nonstandard python package used is biopython.

The scripts require at least python >=3.7, but I tested everything with python 3.10.12

Happy aligning!
--Samuel Decker
