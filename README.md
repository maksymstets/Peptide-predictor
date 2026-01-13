# Peptide predictor
####Video Demo:
That is a script to predict peptides obtained after an action of different proteinases. This software  was heavily inspired by PeptideCutter tool from Expasy (https://web.expasy.org/peptide_cutter/).

##Input/output
The script can accept FASTA files with one or multiple records as an input. The output of the script is PDF file with the results of protein digestion. Each file is autonamed: "Enzyme_chosen_peptides_Digestion_mode_from_Protein_description". That means that user can easily distinguish the different files just by reading the name of the file. It already contains the nme of the enzyme which digests the target protein, digestion mode, and protein description from FASTA header.  Each file is in separate folder:"Enzyme_chosen_Digestion_mode_Protein_description_datetime". All result folders are in folder "results", which creates during the first run of the program.

##Example of usage
The script provides 3 modes to digest a protein: single digestion, parallel digestion and sequential digestion. All modes need user to input the name of FASTA file ("filename.fasta"), choose the number wich represents the digestion mode and enter the number which represents the enzyme/enzymes wanted. User can input via CLI or via interactive prompt. CLI allows for multiple files to analyse without the need for specific Bash script. Example of CLI input: python project.py --input_files file1.fasta file2.fasta --digestion_mode 1 --input_enzyme_number 5

###Single digestion
Single digestion mode cleaves protein with only one user-selected enzyme. The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship.pip install -r requirements.txt

###Parallel digestion
Parallel digestion mode cleaves protein with 2 or more user-selected enzymes and writes results in order of selection. It allows to compare the results of digestion by different enzymes. The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship for both enzymes selected. First page of resulting PDF contains the clickable table of contents to allow user easily orient between pages.

###Sequential digestion
Sequential digestion cleaves protein with 2 or more user-selected enzymes one by one: the peptides from first enzymes are substrate for a second one and so on. The outfile contains the number of peptides obtained, the quantity of cleavage sites, the list of an actual peptides, the cleavages sites list and the site-peptide relationship for first enzyme. For every next enzyme file contains the name of corresponding enzyme, the number of peptides produced in this step of digestion, the list of peptides obtained and per-peptide breakdown: original peptide, produced peptide, the number of cleavage site from the original peptide. If secondary digestion did in fact occur then the produced peptide is blue and cleavage sites are red.  First page of resulting PDF contains the enzyme-selected, number of digestions and clickable table of contents to allow user easily orient between pages.
 
##Requirements:
Python 3.10.12 (no tests conducted to assure compatibility with previous or next versions)
biopython==1.86
fpdf2==2.8.5
pytest==9.0.2




