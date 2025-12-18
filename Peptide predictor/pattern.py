
def main():
    input_file = input("Enter the input FASTA file name: ")
    
    formatted_fasta_file = "formatted_sequence.txt"
    remove_header_fasta(input_file, formatted_fasta_file)
    enzyme_register = {"Arg-C proteinase": argc_proteinase_cutter, "Trypsin": trypsin_cutter, "Asp-N endopeptidase": aspn_endopeptidase_cutter, "Chymotrypsin (high specificity)": chymotrypsinh_cutter, "Chymotrypsin (low specificity)": chymotrypsinl_cutter, "Pepsin_pH1.3": pepsin_pH_1_3_cutter,"Pepsin_pH_2":pepsin_pH_2_cutter, "Papain": papain_cutter }
    input_enzyme = input("""Choose the enzyme frome the list below. 
    Arg-C proteinase 
    Asp-N endopeptidase                     
    Trypsin 
    Chymotrypsin (high specificity)
    Chymotrypsin (low specificity)
    Pepsin_pH1.3
    Pepsin_pH_2
    Papain
    Enter the enzyme chosen exactly as shown in the list:""")
    if input_enzyme  in enzyme_register:
        enzyme_selection = enzyme_register[input_enzyme]
        enzyme_selection(formatted_fasta_file, peptides=f"{input_enzyme.replace(' ', '_')}_peptides.txt")
    else:
        print("Please choose a valid enzyme. Write the name exactly as shown in the list.")
        

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
    #This function splits the sequence according to the cleavage site. Mimics trypsin behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P2 = content[site - 1] if site > 0 else None
                P1 = content[site]
                P1_prime = content[site + 1]
                if (P1 == 'K' or P1 == 'R'):
                    cleavage = True
                    if P1_prime == 'P':
                        if not (P1_prime == 'K' and P2 == 'W'):
                            cleavage = False
                        elif not (P1_prime == 'R' and P2 == 'M'):
                            cleavage = False
                    if P1 == 'K' and P1_prime == 'D' and (P2 == 'C' or P2 == 'D'):
                        cleavage = False
                    if P1 == 'K' and (P1_prime == 'H' or P1_prime == 'Y') and P2 == 'C':
                        cleavage = False
                    if P1 == 'R' and P1_prime == 'K' and P2 == 'C':
                        cleavage = False
                    if P1 == 'R' and (P1_prime == 'H' or P1_prime == 'R') and P2 == 'R':
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)
            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")
            print(f"The number of peptides: {len(peptides_list)}")
            
            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile
        
def argc_proteinase_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics Arg-C proteinase behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                if P1 == 'R':
                    cleavage = True
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile

def aspn_endopeptidase_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics Asp-N endopeptidase behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1_prime = content[site + 1]
                if P1_prime == 'D':
                    cleavage = True
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile
        
def chymotrypsinh_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics Chymotrypsin (high specificity) behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                P1_prime = content[site + 1]
                if (P1 == "F" or P1 == "W" or P1 == "Y"):
                    cleavage = True
                    if P1_prime == "P":
                        cleavage = False
                    if P1 == "W" and P1_prime == "M":
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile
        
def chymotrypsinl_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics Chymotrypsin (low specificity) behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                P1_prime = content[site + 1]
                if (P1 == "F" or P1 == "W" or P1 == "Y" or P1 == "L" or P1 == "M" or P1 == "H"):
                    cleavage = True
                    if P1_prime == "P":
                        cleavage = False
                    if P1 == "W" and P1_prime == "M":
                        cleavage = False
                    if P1 == "M" and P1_prime == "Y":
                        cleavage = False
                    if P1 == "H" and (P1_prime == "D" or P1_prime == "M" or P1_prime == "W"):
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile

def pepsin_pH_1_3_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics pepsin specificity under pH 1.3.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                P2 = content[site - 1] if site > 0 else None
                P3 = content[site - 2] if site > 1 else None
                P1_prime = content[site + 1]
                P2_prime = content[site + 2] if site + 2 < len(content) else None
                if (P1 == "F" or P1 == "L") or (P1_prime == "F" or P1_prime == "L"):
                    cleavage = True
                    if P2_prime == "P" or P2 == "P" or (P3 == "H" or P3 == "R" or P3 == "K") or P1 == "R":
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile
        
def pepsin_pH_2_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics pepsin specificity under pH>2.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                P2 = content[site - 1] if site > 0 else None
                P3 = content[site - 2] if site > 1 else None
                P1_prime = content[site + 1]
                P2_prime = content[site + 2] if site + 2 < len(content) else None
                if (P1 == "F" or P1 == "L" or P1 == "W" or P1 == "Y") or (P1_prime == "F" or P1_prime == "L" or P1_prime == "W" or P1_prime == "Y"):
                    cleavage = True
                    if P2_prime == "P" or P2 == "P" or (P3 == "H" or P3 == "R" or P3 == "K") or P1 == "R":
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile

def papain_cutter(enzyme_selected, peptides):
    #This function splits the sequence according to the cleavage site. Mimics papain behaviour.
    cleavage_sites = []
    with open(enzyme_selected, "r") as infile:
        content = infile.read()
        with open(peptides, "w") as outfile:
            for site in range(len(content) - 1):
                cleavage = False
                P1 = content[site]
                P2 = content[site - 1] if site > 0 else None
                P1_prime = content[site + 1]
                if (P1 == "R" or P1 == "K") and (P2 == "A" or P2 == "V" or P2 == "L" or P2 == "I" or P2 == "F" or P2 == "W" or P2 == "Y"):
                    cleavage = True
                    if P1_prime == "V":
                        cleavage = False
                if cleavage:
                    cleavage_sites.append(site + 1)

            positions = [0] + cleavage_sites + [len(content)]
            peptides_list = []
            for i in range(len(positions) - 1):
                peptides_list.append(content[positions[i]:positions[i + 1]])

            print(f"The number of peptides: {len(peptides_list)}")
            print(f"The quantity of cleavage sites: {len(cleavage_sites)}")

            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
            return outfile
if __name__ == "__main__":
    main()