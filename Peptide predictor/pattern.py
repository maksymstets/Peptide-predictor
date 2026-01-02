import sys

def main():
    input_file = input("Enter the input FASTA file name: ")
    formatted_fasta_file = "formatted_sequence.txt"
    remove_header_fasta(input_file, formatted_fasta_file)
    sequence = sequence_maker(formatted_fasta_file)
    enzyme_register = {"Arg-C proteinase": argc_proteinase_logic,
                        "Trypsin": trypsin_logic,
                        "Asp-N endopeptidase": aspn_endopeptidase_logic,
                        "Chymotrypsin (high specificity)": chymotrypsinh_logic,
                        "Chymotrypsin (low specificity)": chymotrypsinl_logic, 
                        "Pepsin_pH1.3": pepsin_pH_1_3_logic,
                        "Pepsin_pH_2":pepsin_pH_2_logic, 
                        "Papain": papain_logic}
    
    digestion_mode = input("""Select digestion mode from the list below.
    1. Single  digestion
    2. Parallel  digestion
    3. Sequential digestion
    Enter the number of digestion mode:""")

    if digestion_mode == "1":
        input_enzyme = input("""Choose the enzyme from the list below. 
    Arg-C proteinase 
    Asp-N endopeptidase                     
    Trypsin 
    Chymotrypsin (high specificity)
    Chymotrypsin (low specificity)
    Pepsin_pH1.3
    Pepsin_pH_2
    Papain
    Enter the enzyme chosen exactly as shown in the list:""")
        if input_enzyme in enzyme_register:
            single_digestion(sequence, input_enzyme, enzyme_register)
        else:
            sys.exit("Please choose a valid enzyme. Write the name exactly as shown in the list.")
            
    elif digestion_mode == "2":
        input_enzyme = input("""Choose the enzyme from the list below. 
    Arg-C proteinase 
    Asp-N endopeptidase                     
    Trypsin 
    Chymotrypsin (high specificity)
    Chymotrypsin (low specificity)
    Pepsin_pH1.3
    Pepsin_pH_2
    Papain
    Enter the enzyme chosen exactly as shown in the list, separated by a comma:""")
        enzyme_names = []
        for enzyme in input_enzyme.split(","):
            enzymes_chosen = enzyme.strip()
            enzyme_names.append(enzymes_chosen)
        
        correct_enzymes = []
        for enzyme in enzyme_names:
            if enzyme in enzyme_register:
                correct_enzymes.append(enzyme)
            else:
                print(f"{enzyme} is not a valid enzyme. Please write the name exactly as shown in the list.")
        if correct_enzymes:
            parallel_digestion(sequence, input_enzyme, correct_enzymes, enzyme_register)
            
        else:
            print("No valid enzymes were selected for parallel digestion.")
        
    elif digestion_mode == "3":
        input_enzyme = input("""Choose the enzyme from the list below. 
    Arg-C proteinase 
    Asp-N endopeptidase                     
    Trypsin 
    Chymotrypsin (high specificity)
    Chymotrypsin (low specificity)
    Pepsin_pH1.3
    Pepsin_pH_2
    Papain
    Enter the enzyme chosen exactly as shown in the list, separated by a comma:""")
        enzyme_names = []
        for enzyme in input_enzyme.split(","):
            enzymes_chosen = enzyme.strip()
            enzyme_names.append(enzymes_chosen)

        correct_enzymes = []
        for enzyme in enzyme_names:
            if enzyme in enzyme_register:
                correct_enzymes.append(enzyme)
            else:
                print(f"{enzyme} is not a valid enzyme. Please write the name exactly as shown in the list.")
        if correct_enzymes:
            _, operations_log = sequential_digestion(sequence, correct_enzymes, enzyme_register)
            output_filename = f'{"_".join(enzyme_names)}_peptides_sequential_digest.txt'
            sequential_output_writer(operations_log, enzyme_names, output_filename)
            
        else:
            print("No valid enzymes were selected for sequential digestion.")
    else:
        print("Please choose a valid digestion mode from the list.")
   

def remove_header_fasta(input_file, formatted_fasta_file):
    # This function reads a FASTA file and writes the sequence data to a new file, omitting the header lines.
    try:
        with open(input_file, "r") as infile:
            with open(formatted_fasta_file, "w") as outfile:
                for line in infile:
                    if not line.startswith('>'):
                        outfile.write(line.strip())

    except FileNotFoundError:
        sys.exit("The file was not found. Check the name of the file.")

def sequence_maker(formatted_fasta_file):
    # This function reads the formatted FASTA file and returns the sequence as a single string.
    try:
        with open(formatted_fasta_file, "r") as infile:
            content = infile.read()
            return str(content)
    except FileNotFoundError:
        sys.exit("The file was not found. Checm the name of the file")

def single_digestion(sequence, input_enzyme, enzyme_register):
    enzyme_selection = enzyme_register[input_enzyme]
    cleavage_sites_finder = enzyme_selection(sequence)
    peptides_list_generator = peptide_generator(cleavage_sites_finder, sequence)
    output_filename = f"{input_enzyme.replace(' ', '_')}_peptides_single_digest.txt"
    output_writer_single(peptides_list_generator, cleavage_sites_finder, output_filename)
        
def parallel_digestion(sequence, input_enzyme, enzyme_names, enzyme_register):
    all_results = []
    for input_enzyme in enzyme_names:
        enzyme_selection = enzyme_register[input_enzyme]
        cleavage_sites_finder = enzyme_selection(sequence)
        peptides_list_generator = peptide_generator(cleavage_sites_finder, sequence)
    
        all_results.append({
            'enzyme': input_enzyme,
            'cleavage_sites': cleavage_sites_finder,
            'peptides': peptides_list_generator
       })
    output_filename = f'{"_".join(enzyme_names)}_peptides_parallel_digest.txt'
    parallel_output_writer(all_results, output_filename)
    
def sequential_digestion(sequence, enzyme_names, enzyme_register):
    current_sequence = [sequence]
    operations_log = []

    for step_num, enzyme_name in enumerate(enzyme_names, 1):
        enzyme_selection = enzyme_register[enzyme_name]
        next_gen_peptides = []

        if step_num == 1:
            # digest the full sequence once
            cleavage_sites = enzyme_selection(sequence)
            next_gen_peptides = peptide_generator(cleavage_sites, sequence)
            operations_log.append({
                'enzyme': enzyme_name,
                'cleavage_sites': cleavage_sites,
                'peptides_list': next_gen_peptides,
                'is_first': True
            })
        else:
            # digest each peptide produced in the previous step and collect all new peptides
            per_peptide_details = []
            for peptide in current_sequence:
                cleavage_sites = enzyme_selection(peptide)
                if cleavage_sites:
                    produced = peptide_generator(cleavage_sites, peptide)
                else:
                    produced = [peptide]
                per_peptide_details.append({
                    'former_peptide': peptide,
                    'cleavage_sites': cleavage_sites,
                    'produced': produced
                })
                next_gen_peptides.extend(produced)

            # store aggregated list of all peptides produced at this step
            operations_log.append({
                'enzyme': enzyme_name,
                'is_first': False,
                'all_peptides': next_gen_peptides,
                'per_peptide': per_peptide_details
            })

        current_sequence = next_gen_peptides

    return current_sequence, operations_log
    
def trypsin_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics trypsin behaviour.
    cleavage_sites = []
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
    return cleavage_sites
         
def argc_proteinase_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics Arg-C proteinase behaviour.
    cleavage_sites = []
    for site in range(len(content) - 1):
        cleavage = False
        P1 = content[site]
        if P1 == 'R':
            cleavage = True
        if cleavage:
            cleavage_sites.append(site + 1)
    return cleavage_sites
            
def aspn_endopeptidase_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics Asp-N endopeptidase behaviour.
    cleavage_sites = []
    for site in range(len(content) - 1):
        cleavage = False
        P1_prime = content[site + 1]
        if P1_prime == 'D':
            cleavage = True
        if cleavage:
            cleavage_sites.append(site + 1)
    return cleavage_sites

def chymotrypsinh_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics Chymotrypsin (high specificity) behaviour.
    cleavage_sites = []
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
    return cleavage_sites
            
def chymotrypsinl_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics Chymotrypsin (low specificity) behaviour.
    cleavage_sites = []
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
    return cleavage_sites

def pepsin_pH_1_3_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics pepsin specificity under pH 1.3.
    cleavage_sites = []
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
    return cleavage_sites
          
def pepsin_pH_2_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics pepsin specificity under pH>2.
    cleavage_sites = []
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
    return cleavage_sites
            
def papain_logic(content):
    #This function splits the sequence according to the cleavage site. Mimics papain behaviour.
    cleavage_sites = []
    
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
    return cleavage_sites
       
def peptide_generator(cleavage_sites, content):
    positions = [0] + cleavage_sites + [len(content)]
    peptides_list = []
    for i in range(len(positions) - 1):
        peptides_list.append(content[positions[i]:positions[i + 1]])
    return peptides_list

def output_writer_single(peptides_list, cleavage_sites, output_filename):
    with open(output_filename, "w") as outfile:
        outfile.write(f"The quantity of cleavage sites: {len(cleavage_sites)}\n\n")
        outfile.write(f"The number of peptides: {len(peptides_list)}\n\n")
        outfile.write(f"The peptides list: {peptides_list}\n\n")
        outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
        outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}")
        print(f"Results written to {output_filename}")
        return outfile

def parallel_output_writer(all_results, output_filename):
    with open(output_filename, "w") as outfile:
        for _, single_result in enumerate(all_results):
            
            enzyme_name = single_result['enzyme']
            peptides_list = single_result['peptides']
            cleavage_sites = single_result['cleavage_sites']
            outfile.write(f"Enzyme: {enzyme_name}\n\n")
            outfile.write(f"The quantity of cleavage sites: {len(cleavage_sites)}\n\n")
            outfile.write(f"The number of peptides: {len(peptides_list)}\n\n")
            outfile.write(f"The peptides list: {peptides_list}\n\n")
            outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
            outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}\n\n")
        print(f"Results written to {output_filename}")
        return outfile

def sequential_output_writer(operations_log, enzyme_names, output_filename):
    
    with open(output_filename, "w") as outfile:

        outfile.write(f"Enzyme cascade: {' → '.join(enzyme_names)}\n")
        outfile.write(f"Number of digestions: {len(operations_log)}\n\n")

        for digestion in operations_log:
            if digestion.get('is_first'):
                cleavage_sites = digestion.get('cleavage_sites', [])
                peptides_list = digestion.get('peptides_list', [])
                outfile.write(f"The quantity of cleavage sites: {len(cleavage_sites)}\n\n")
                outfile.write(f"The number of peptides: {len(peptides_list)}\n\n")
                outfile.write(f"The peptides list: {peptides_list}\n\n")
                outfile.write(f"The cleavage sites list: {cleavage_sites}\n\n")
                outfile.write(f"The site-peptide relationship: {dict(zip(peptides_list, cleavage_sites))}\n\n")
            else:
                enzyme_name = digestion.get('enzyme', 'Unknown')
                all_peptides = digestion.get('all_peptides', [])
                outfile.write(f"Enzyme: {enzyme_name}\n\n")
                outfile.write(f"The number of peptides produced in this step: {len(all_peptides)}\n\n")
                outfile.write(f"All peptides (aggregated) from this step: {all_peptides}\n\n")
                # optional: also print per-peptide breakdown for debugging / detail
                per_peptide = digestion.get('per_peptide', [])
                if per_peptide:
                    outfile.write("Per-peptide breakdown:\n")
                    for entry in per_peptide:
                        outfile.write(f"  From '{entry['former_peptide']}' -> {entry['produced']} (cleavage sites: {entry['cleavage_sites']})\n")
                    outfile.write("\n")
                
        print(f"Results written to {output_filename}")
        return outfile

if __name__ == "__main__":
    main()