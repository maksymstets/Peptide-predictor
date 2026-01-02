Peptide-predictor

That is a script to predict peptides obtained after an action of different proteinases. This software was heavily inspired by PeptideCutter tool from Expasy (https://web.expasy.org/peptide_cutter/).

Input/output 

The script can accept FASTA files with one record as an input. The output of the script creates 2 text files: "formatted_sequence.txt", which contains pure aminoacid sequence without FASTA header and file with the results of protein digestion.

Example of usage 

The script provides 3 ways to digest a protein: single digestion, parallel digestion and sequential digestion.

Single digestion mode cleaves protein with only one user-selected enzyme. The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship.

Parallel digestion mode cleaves protein with 2 or more user-selected enzymes and writes results in order of selection. It allows to compare the results of digestion by different enzymes. The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship for both enzymes selected.

Sequential digestion cleaves protein with 2 or more user-selected enzymes one by one: the peptides from first enzymes are substrate for a second one and so on. The outfile displays the peptides obtained as well as sites of cleavage (numeration from an original sequence). The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship for first enzyme. For every next enzyme file contains the name of corrsponding enzyme, the number of peptides produced in this step of digestion, the list of peptides obtained and per-peptide breakdown: "parent peptide", "child peptide", the number of cleavage site from an original sequence.

All modes ask user to input the name of FASTA file ("filename.fasta"), choose the number wich represents the digestion mode and enter the name of enzyme/enzymes wanted.

Requirements: Python 3.10.12 (no tests conducted to assure compatibility with previous or next versions)
