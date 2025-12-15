import re
def main():
    input_file = input("Enter the input FASTA file name: ")
    
    formatted_fasta_file = "formatted_sequence.txt"
    remove_header_fasta(input_file, formatted_fasta_file)
    enzyme_register = {"Arg-C proteinase": argc_proteinase_cutter, "Trypsin": trypsin_cutter}
    input_enzyme = input("""Choose the enzyme frome the list below. 
    Arg-C proteinase 
    Trypsin 
    Enter the enzyme chosen:""")
    if input_enzyme  in enzyme_register:
        enzyme_selection = enzyme_register[input_enzyme]
        enzyme_selection(formatted_fasta_file, peptides=f"{input_enzyme.replace(' ', '_')}_peptides.txt")
    else:
        print("Please choose a valid enzyme.")
        

def remove_header_fasta(input_file, formatted_fasta_file):
    # This function reads a FASTA file and writes the sequence data to a new file, omitting the header lines.
    try:
        with open(input_file, "r") as infile:
            with open(formatted_fasta_file, "w") as outfile:
                for line in infile:
                    if not line.startswith('>'):
                        outfile.write(line.strip())

    except FileNotFoundError:
        print("The file was not found.")

def trypsin_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site.
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            peptides_list = re.split(r'(?<=[RK])(?!P)', content)
            print(len(peptides_list))
            outfile.write(str(peptides_list))
            return outfile
        
def argc_proteinase_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site.
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            peptides_list = re.split(r'(?<=R)', content)
            print(len(peptides_list))
            outfile.write(str(peptides_list))
            return outfile

if __name__ == "__main__":
    main()