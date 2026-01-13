import sys
import argparse
import os
from fpdf import FPDF
from Bio import SeqIO
from datetime import datetime

class PDF(FPDF):
    """
    Custom PDF class extending FPDF with automatic footer generation.
    
    Features:
    - Adds page numbers centered at the bottom
    - Includes clickable "Back to Table of Contents" link
    - Footer can be toggled on/off via show_footer flag
    - Skips footer on page 1 (title page)
    
    Attributes:
        toc_link (int): Link ID for the Table of Contents page
        show_footer (bool): Flag to control footer visibility
    """
    def __init__(self, toc_link = None):
        """
        Initialize custom PDF with optional TOC link.
        
        Args:
            toc_link (int, optional): Link ID for Table of Contents. Defaults to None.
        """
        super().__init__()
        self.toc_link = toc_link
        self.show_footer = False
    def footer(self):
        """
        Automatically called by FPDF at the bottom of each page.
        
        Generates a footer with:
        - Centered page number
        - Right-aligned "Back to Table of Contents" link (if toc_link exists)
        
        Footer is skipped on page 1 and when show_footer is False.
        """
        # Skip footer on first page (title/TOC page)
        if self.page_no() == 1:
            return
        if self.show_footer is False:
            return
        
        # Page number and footer position
        self.set_y(-15)
        self.set_font("Courier", 'I', 8)
        self.set_text_color(0, 0, 0)
        self.cell(0, 5, f"Page {self.page_no()}", align="C", new_x = "RIGHT", new_y = "TOP")

        # "Back to TOC" link, oritentation and color
        if self.toc_link is not None:
            self.set_text_color(0, 0, 255)
            self.cell(0, 5, 'Back to Table of Contents', align='R', link=self.toc_link, new_x="LMARGIN", new_y="NEXT")
            self.set_text_color(0, 0, 0)

def regular_input_storer():
    """
    Interactive mode: Prompt user for input via command line.
    
    Collects:
    - Input FASTA file name
    - Digestion mode (single/parallel/sequential)
    - Enzyme number(s)
    
    Returns:
        tuple: (input_file: str, digestion_mode: str, input_enzyme_number: str)
    """
    input_file = input("Enter the input FASTA file name: ")
    print("\n\n")
    digestion_mode = input("""Select digestion mode from the list below.
    1. Single  digestion
    2. Parallel  digestion
    3. Sequential digestion
    Enter the number of digestion mode: """)
    print("\n\n")
    input_enzyme_number = input("""Choose the enzyme from the list below. 
    1. Arg-C_proteinase 
    2. Asp-N_endopeptidase 
    3. Chymotrypsin_high_specificity
    4. Chymotrypsin_low_specificity                                             
    5. Trypsin 
    6. Papain
    7. Pepsin_pH1.3
    8. Pepsin_pH_2
    
    Enter the enzyme number:""")

    return input_file, digestion_mode, input_enzyme_number

def cli_input_storer():
    """
    CLI mode: Parse command-line arguments using argparse.
    
    Supports batch processing of multiple FASTA files.
    
    Returns:
        tuple: (input_files: list, digestion_mode: str, input_enzyme_number: str)
    
    Example:
        python project.py --input_files file1.fasta file2.fasta --digestion_mode 1 --input_enzyme_number 5
    """
    parser = argparse.ArgumentParser(description="Peptide predictor")
    parser.add_argument("--input_files", nargs = "+", help="Input FASTA file name or filenames")
    parser.add_argument("--digestion_mode",choices=["1", "2", "3"], help="""Choose digestion mode: \
    1. Single  digestion
    2. Parallel  digestion
    3. Sequential digestion""")
    parser.add_argument("--input_enzyme_number", help="""Enter enzyme number or list of enzyme numbers separated by commas:
    1. Arg-C_proteinase 
    2. Asp-N_endopeptidase 
    3. Chymotrypsin_high_specificity
    4. Chymotrypsin_low_specificity                                               
    5. Trypsin 
    6. Papain
    7. Pepsin_pH1.3
    8. Pepsin_pH_2""")

    args = parser.parse_args()

    return args.input_files, args.digestion_mode, args.input_enzyme_number

def input_validator(digestion_mode, input_enzyme_number, sequence, protein_description ):
    """
    Validate user input and route to appropriate digestion function.
    
    This function acts as the main controller:
    1. Maps enzyme numbers to enzyme names
    2. Validates enzyme selections
    3. Routes to single/parallel/sequential digestion based on mode
    
    Args:
        digestion_mode (str): "1" for single, "2" for parallel, "3" for sequential
        input_enzyme_number (str): Comma-separated enzyme numbers (e.g., "3,5")
        sequence (str): Protein amino acid sequence
        protein_description (str): Sanitized protein name for output files
    
    Raises:
        SystemExit: If invalid enzyme numbers or digestion mode provided
    """
    # Map enzyme names to their logic functions
    enzyme_register = {"Arg-C_proteinase": argc_proteinase_logic,
                        "Trypsin": trypsin_logic,
                        "Asp-N_endopeptidase": aspn_endopeptidase_logic,
                        "Chymotrypsin_high_specificity": chymotrypsinh_logic,
                        "Chymotrypsin_low_specificity": chymotrypsinl_logic, 
                        "Pepsin_pH1.3": pepsin_pH_1_3_logic,
                        "Pepsin_pH_2":pepsin_pH_2_logic, 
                        "Papain": papain_logic}
    # Map user input numbers to enzyme names
    enzyme_dict = {"1": "Arg-C_proteinase",
                   "2": "Asp-N_endopeptidase",
                   "3": "Chymotrypsin_high_specificity",
                   "4": "Chymotrypsin_low_specificity",
                   "5": "Trypsin",
                   "6": "Papain",
                   "7": "Pepsin_pH1.3",
                   "8": "Pepsin_pH_2"}

    # Route to appropriate digestion function based on mode (single)
    if digestion_mode == "1":
        
        if input_enzyme_number in enzyme_dict:
            input_enzyme = enzyme_dict[input_enzyme_number]
            if input_enzyme in enzyme_register:
                single_digestion(sequence, input_enzyme, enzyme_register, protein_description)
        else:
            sys.exit("Please choose a valid enzyme number from the list or check if the input is not a sequence of numbers.")
             
    # Route to appropriate digestion function based on mode (parallel)         
    elif digestion_mode == "2":
        # Parse comma-separated enzyme numbers
        enzyme_numbers = []
        for key in input_enzyme_number.split(","):
            keys = key.strip()
            enzyme_numbers.append(keys)
        # Convert numbers to enzyme names     
        enzyme_names = []
        for input_enzyme in enzyme_numbers:
            if input_enzyme in enzyme_dict:
                enzyme_number = enzyme_dict[input_enzyme]
                enzyme_names.append(enzyme_number)
            else:
                sys.exit("Please choose a valid enzyme number from the list.")
         # Validate enzyme names exist in register
        correct_enzymes = []
        for enzyme in enzyme_names:  
            if enzyme in enzyme_register:
                correct_enzymes.append(enzyme)
            else:
                print(f"{enzyme} is not a valid enzyme. Please write the name exactly as shown in the list.")
        if correct_enzymes:
            parallel_digestion(sequence, enzyme_names, enzyme_register, protein_description)
            
        else:
            print("No valid enzymes were selected for parallel digestion.")

    # Route to appropriate digestion function based on mode (sequential)
    elif digestion_mode == "3":
        # Parse comma-separated enzyme numbers
        enzyme_numbers = []
        for key in input_enzyme_number.split(","):
            keys = key.strip()
            enzyme_numbers.append(keys)
        # Convert numbers to enzyme names  
        enzyme_names = []
        for input_enzyme in enzyme_numbers:
            if input_enzyme in enzyme_dict:
                enzyme_number = enzyme_dict[input_enzyme]
                enzyme_names.append(enzyme_number)
            else:
                sys.exit("Please choose a valid enzyme number from the list.")
         # Validate enzyme names exist in register
        correct_enzymes = []
        for enzyme in enzyme_names:
            if enzyme in enzyme_register:
                correct_enzymes.append(enzyme)
            else:
                print(f"{enzyme} is not a valid enzyme. Please write the name exactly as shown in the list.")
        if correct_enzymes:
            _, operations_log = sequential_digestion(sequence, correct_enzymes, enzyme_register)
            
            sequential_output_writer(operations_log, enzyme_names, protein_description)
            
        else:
            print("No valid enzymes were selected for sequential digestion.")
    else:
        print("Please choose a valid digestion mode from the list.")

def parse_multi_fasta(input_file):
    """
    Parse a multi-FASTA file using Biopython and return protein data.
    
    Args:
        input_file (str): Path to FASTA file.
    
    Returns:
        list: List of dictionaries with protein information.
    
    Example:
        [
            {
                'id': 'ABR21772.1',
                'description': 'conglutin beta [Lupinus angustifolius]',
                'sequence': 'MLPML...'
            },
            ...
        ]
    """
    proteins = []
    
    try:
        for record in SeqIO.parse(input_file, "fasta"):
            proteins.append({
                'accession_number': record.id, 
                'description': record.description.replace(record.id, '').strip(),  
                'sequence': str(record.seq)  
            })
        
        return proteins
        
    except FileNotFoundError:
        sys.exit(f"Error: File '{input_file}' not found.")
    except Exception as e:
        sys.exit(f"Error parsing FASTA file: {e}")

def filename_converter(name):
    """
    Sanitize text for use in filenames and folder names.
    
    Removes or replaces characters that are problematic in filesystems:
    - Spaces → underscores
    - Brackets, parentheses → removed
    - Path separators (forward and backslash) → underscores
    - Special characters (: ; , .) → underscores or removed
    
    Args:
        name (str): Original text (e.g., protein description)
    
    Returns:
        str: Filesystem-safe string
    
    Example:
        >>> filename_converter("conglutin beta [Lupinus angustifolius]")
        'conglutin_beta_Lupinus_angustifolius'
    """
    if not name:
        return "unknown"
    
    # Remove or replace problematic characters
    clean = (name.replace(' ', '_')
                  .replace('(', '')
                  .replace(')', '')
                  .replace('[', '')
                  .replace(']', '')
                  .replace('.', '_')
                  .replace(',', '')
                  .replace('/', '_')
                  .replace('\\', '_')
                  .replace(':', '_')
                  .replace(';', '_'))
    if clean:
        return clean

def single_digestion(sequence, input_enzyme, enzyme_register, protein_description):
    """
    Perform single-enzyme digestion on a protein sequence.
    
    Workflow:
    1. Get enzyme cleavage logic function from register
    2. Find cleavage sites in sequence
    3. Generate peptide fragments
    4. Write results to PDF
    
    Args:
        sequence (str): Protein amino acid sequence
        input_enzyme (str): Enzyme name (e.g., "Trypsin")
        enzyme_register (dict): Mapping of enzyme names to logic functions
        protein_description (str): Sanitized protein name for output files
    """
    enzyme_selection = enzyme_register[input_enzyme]
    cleavage_sites_finder = enzyme_selection(sequence)
    peptides_list_generator = peptide_generator(cleavage_sites_finder, sequence)
    output_writer_single(input_enzyme, peptides_list_generator, cleavage_sites_finder, protein_description)
           
def parallel_digestion(sequence, enzyme_names, enzyme_register, protein_description):
    """
    Perform parallel digestion with multiple enzymes.
    
    Each enzyme digests the ORIGINAL sequence independently.
    Results are compiled into a single PDF with separate pages per enzyme.
    
    Workflow:
    1. For each enzyme:
       - Find cleavage sites in original sequence
       - Generate peptide fragments
       - Store results
    2. Compile all results into one PDF with Table of Contents
    
    Args:
        sequence (str): Protein amino acid sequence
        enzyme_names (list): List of enzyme names (e.g., ["Trypsin", "Papain"])
        enzyme_register (dict): Mapping of enzyme names to logic functions
        protein_description (str): Sanitized protein name for output files
    """
    # Digest with each enzyme independently
    all_results = []
    for enzyme_name in enzyme_names:
        enzyme_selection = enzyme_register[enzyme_name]
        cleavage_sites_finder = enzyme_selection(sequence)
        peptides_list_generator = peptide_generator(cleavage_sites_finder, sequence)
    
        all_results.append({
            'enzyme': enzyme_name,
            'cleavage_sites': cleavage_sites_finder,
            'peptides': peptides_list_generator
       })
    parallel_output_writer(all_results, enzyme_names, protein_description)
    
def sequential_digestion(sequence, enzyme_names, enzyme_register):
    """
    Perform sequential digestion (enzyme cascade).
    
    Each enzyme digests the peptides produced by the PREVIOUS enzyme.
    This simulates a multi-step enzymatic process.
    
    Workflow:
    1. First enzyme digests original sequence
    2. Second enzyme digests all peptides from step 1
    3. Third enzyme digests all peptides from step 2
    4. And so on...
    
    Args:
        sequence (str): Protein amino acid sequence
        enzyme_names (list): Ordered list of enzyme names (e.g., ["Trypsin", "Papain"])
        enzyme_register (dict): Mapping of enzyme names to logic functions
    
    Returns:
        tuple: (final_peptides: list, operations_log: list)
            - final_peptides: List of peptide fragments after all digestions
            - operations_log: Detailed record of each digestion step
    """
    current_sequence = [sequence]
    operations_log = []

    for step_num, enzyme_name in enumerate(enzyme_names, 1):
        enzyme_selection = enzyme_register[enzyme_name]
        next_gen_peptides = []
        # digest the full sequence once
        if step_num == 1:
            
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
        # Update current sequence for next iteration
        current_sequence = next_gen_peptides

    return current_sequence, operations_log
    
def trypsin_logic(content):
    """
    Simulate trypsin cleavage specificity.
    
    Trypsin cleaves peptide bonds C-terminal to lysine (K) and arginine (R),
    except when followed by proline (P) or in specific inhibitory contexts.
    
    Cleavage rules:
    - Cleaves after K or R
    - Does NOT cleave if followed by P (except in specific K-P/R-P contexts)
    - Additional exceptions for specific sequence contexts (CKD, CKH, etc.)
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
        
    Reference:
        Based on EXPASY (PeptideCutter) peptidase database trypsin specificity
    """
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
    """
    Simulate Arg-C proteinase cleavage specificity.
    
    Arg-C proteinase cleaves peptide bonds C-terminal to arginine (R) only.
    Simple, high-specificity enzyme.
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate Asp-N endopeptidase cleavage specificity.
    
    Asp-N endopeptidase cleaves peptide bonds N-terminal to aspartic acid (D).
    This is an N-terminal cleavage enzyme (unusual compared to most proteases).
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate chymotrypsin cleavage specificity (high specificity mode).
    
    Chymotrypsin (high specificity) cleaves C-terminal to large hydrophobic
    aromatic amino acids: phenylalanine (F), tryptophan (W), and tyrosine (Y).
    
    Exceptions:
    - Does NOT cleave if followed by proline (P)
    - Does NOT cleave at WM (tryptophan-methionine)
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate chymotrypsin cleavage specificity (low specificity mode).
    
    Chymotrypsin (low specificity) cleaves C-terminal to a broader range
    of hydrophobic amino acids: F, W, Y, L, M, and H.
    
    Exceptions:
    - Does NOT cleave if followed by proline (P)
    - Does NOT cleave at WM (tryptophan-methionine)
    - Does NOT cleave at MY (methionine-tyrosine)
    - Does NOT cleave at HD, HM, or HW (histidine followed by D, M, or W)
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate pepsin cleavage specificity at pH 1.3.
    
    Pepsin at pH 1.3 preferentially cleaves peptide bonds involving
    hydrophobic amino acids, particularly phenylalanine (F) and leucine (L).
    
    Cleavage rules:
    - Cleaves when P1 or P1' is F or L
    - Exceptions: Does NOT cleave if:
      * P2' is proline
      * P2 is proline
      * P3 is H, R, or K (charged residues)
      * P1 is R
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate pepsin cleavage specificity at pH > 2.
    
    Pepsin at pH > 2 has broader specificity, cleaving at hydrophobic
    amino acids: F, L, W, and Y.
    
    Cleavage rules and exceptions are similar to pH 1.3 version,
    but with expanded substrate specificity.
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Simulate papain cleavage specificity.
    
    Papain is a cysteine protease that cleaves peptide bonds C-terminal
    to basic amino acids (R, K) when preceded by hydrophobic residues.
    
    Cleavage rules:
    - P1 must be R or K
    - P2 must be a hydrophobic amino acid (A, V, L, I, F, W, Y)
    - Does NOT cleave if P1' is valine (V)
    
    Args:
        content (str): Protein sequence
    
    Returns:
        list: Zero-indexed positions where cleavage occurs
    """
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
    """
    Generate peptide fragments from sequence and cleavage sites.
    
    Splits the protein sequence at specified cleavage positions to produce
    a list of peptide fragments.
    
    Args:
        cleavage_sites (list): List of cleavage positions (zero-indexed)
        content (str): Protein sequence
    
    Returns:
        list: List of peptide fragments (strings)
    
    Example:
        >>> peptide_generator([2, 5], "ABCDEFG")
        ['AB', 'CDE', 'FG']
    """
    positions = [0] + cleavage_sites + [len(content)]
    peptides_list = []
    for i in range(len(positions) - 1):
        peptides_list.append(content[positions[i]:positions[i + 1]])
    return peptides_list

def output_writer_single(input_enzyme, peptides_list, cleavage_sites, protein_description):
    """
    Generates a PDF report for a single enzyme digestion.

    Creates a timestamped folder in the 'results' directory and saves a PDF 
    containing the enzyme name, statistics on cleavage sites and peptides, 
    and a full list of the resulting peptides.

    Args:
        input_enzyme (str): The name of the enzyme used (e.g., "Trypsin").
        peptides_list (list): A list of peptide strings resulting from digestion.
        cleavage_sites (list): A list of integer indices where cleavage occurred.
        protein_description (str): A cleaned string describing the protein source.

    Returns:
        None: The function saves a file to disk and does not return a value.
    """
     # Build complete path with folder
    timestamp = datetime.now().strftime("%Y_%m_%d-%H_%M")
    folder_name = f"{input_enzyme}_Single_{protein_description}_{timestamp}"
    output_folder = os.path.join("results", folder_name)
    
    # Create folder
    os.makedirs(output_folder, exist_ok=True)
    
    # Build filename
    filename = f"{input_enzyme}_peptides_Single_digest_from_{protein_description}.pdf"
    output_filename = os.path.join(output_folder, filename)
    
    pdf = FPDF()
    pdf.set_font("Courier", 'B', size=12)
    pdf.add_page()
    pdf.cell(0, 5, f"Enzyme: {input_enzyme}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    pdf.set_font("Courier",  size=8)
    pdf.cell(0, 5, f"The quantity of cleavage sites: {len(cleavage_sites)}\n\n", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    pdf.cell(0, 5, f"The number of peptides: {len(peptides_list)}\n\n", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)

    label_peptides_list = "The peptides list:"
    width_peptides_list = pdf.get_string_width(label_peptides_list)
    pdf.cell(width_peptides_list, 5, label_peptides_list,  new_x="RIGHT", new_y="TOP")
    pdf.multi_cell(0, 5, f"{peptides_list}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)

    label_cleavage_list = "The cleavage list:"
    width_cleavage_list = pdf.get_string_width(label_cleavage_list)
    pdf.cell(width_cleavage_list, 5, label_cleavage_list, new_x="RIGHT", new_y="TOP")
    pdf.multi_cell(0, 5, f"{cleavage_sites}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
        
    label_relationship = "The site-peptide relationship:"
    width_relationship = pdf.get_string_width(label_relationship)
    pdf.cell(width_relationship, 5, label_relationship, new_x="RIGHT", new_y="TOP")
    pdf.multi_cell(0, 5, f"{dict(zip(peptides_list, cleavage_sites))}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    print(f"     ✓ Look for results in: results/")
    print(f"        └── {folder_name}/")
    print(f"            └── {filename}\n")
    return pdf.output(output_filename)

def parallel_output_writer(all_results,enzyme_names, protein_description):
    """
    Generates a PDF report for parallel digestion mode.

    Creates a timestamped folder in the 'results' directory and creates a PDF with a Table of Contents linking to separate results for 
    each enzyme. Each section details the cleavage sites and peptides 
    produced by that specific enzyme on the original sequence.

    Args:
        all_results (list): A list of dictionaries, where each dictionary contains:
            - 'enzyme': Name of the enzyme.
            - 'cleavage_sites': List of cut sites.
            - 'peptides': List of resulting peptides.
        enzyme_names (list): A list of names of all enzymes applied.
        protein_description (str): A cleaned string describing the protein source.

    Returns:
        None: The function saves a file to disk.
    """
     # Build complete path with folder
    timestamp = datetime.now().strftime("%Y_%m_%d-%H_%M")
    folder_name = f"{'_'.join(enzyme_names)}_Parallel_{protein_description}_{timestamp}"
    output_folder = os.path.join("results", folder_name)
    
    # Create folder
    os.makedirs(output_folder, exist_ok=True)
    
    # Build filename
    filename = f"{'_'.join(enzyme_names)}_peptides_Parallel_digest_from_{protein_description}.pdf"
    output_filename = os.path.join(output_folder, filename)
    
    toc_link =None
    pdf = PDF(toc_link)
    
    pdf.add_page()
    toc_link = pdf.add_link()
    pdf.set_link(toc_link)
    pdf.toc_link = toc_link
    pdf.set_font("Courier", 'B', size=12)
    pdf.cell(0, 5, "Results for each enzyme are on separate pages below", new_x="LMARGIN", new_y="NEXT" )
    for single_result in all_results:
        single_result['link_id'] = pdf.add_link()

    for single_result in all_results:
        enzyme_name = single_result['enzyme']
        
        pdf.set_font("Courier", size=10)
        pdf.cell(0, 7, f"Jump to Enzyme: {enzyme_name}", link=single_result['link_id'], new_x="LMARGIN", new_y="NEXT")
        
    pdf.show_footer = True

    for single_result in all_results:
        enzyme_name = single_result['enzyme']
        peptides_list = single_result['peptides']
        cleavage_sites = single_result['cleavage_sites']
        
        pdf.add_page()
        pdf.set_link(single_result['link_id'])
        pdf.set_font("Courier", 'B', size=8)
        pdf.cell(0, 5, f"Enzyme: {single_result['enzyme']}\n\n", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        pdf.set_font("Courier", size=8)
        pdf.cell(0, 5, f"The quantity of cleavage sites: {len(cleavage_sites)}\n\n", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        pdf.cell(0, 5, f"The number of peptides: {len(peptides_list)}\n\n", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)

        label_peptides_list = "The peptides list:"
        width_peptides_list = pdf.get_string_width(label_peptides_list)
        pdf.cell(width_peptides_list, 5, label_peptides_list,  new_x="RIGHT", new_y="TOP")
        pdf.multi_cell(0, 5, f"{peptides_list}", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)

        label_cleavage_list = "The cleavage list:"
        width_cleavage_list = pdf.get_string_width(label_cleavage_list)
        pdf.cell(width_cleavage_list, 5, label_cleavage_list, new_x="RIGHT", new_y="TOP")
        pdf.multi_cell(0, 5, f"{cleavage_sites}", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        
        label_relationship = "The site-peptide relationship:"
        width_relationship = pdf.get_string_width(label_relationship)
        pdf.cell(width_relationship, 5, label_relationship, new_x="RIGHT", new_y="TOP")
        pdf.multi_cell(0, 5, f"{dict(zip(peptides_list, cleavage_sites))}", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        pdf.cell(0, 5,"-" * 160, new_x="LMARGIN", new_y="NEXT")
    print(f"     ✓ Look for results in: results/")
    print(f"        └── {folder_name}/")
    print(f"            └── {filename}\n")
    return pdf.output(output_filename)

def sequential_output_writer(operations_log, enzyme_names, protein_description):
    """
    Generates a PDF report for sequential (cascade) digestion mode.

    Creates a timestamped folder in the 'results' directory and creates a detailed PDF tracking the digestion process step-by-step. 
    It includes a summary of the enzyme order and detailed pages for each 
    step, showing how previous peptides were broken down into new ones.

    Args:
        operations_log (list): A list of dictionaries representing each step 
            of the digestion. Contains keys like 'enzyme', 'all_peptides', 
            and 'per_peptide' details.
        enzyme_names (list): A list of enzymes in the order they were applied.
        protein_description (str): A cleaned string describing the protein source.

    Returns:
        None: The function saves a file to disk.
    """
    timestamp = datetime.now().strftime("%Y_%m_%d-%H_%M")
    folder_name = f"{'_'.join(enzyme_names)}_Sequential_{protein_description}_{timestamp}"
    output_folder = os.path.join("results", folder_name)
    
    # Create folder
    os.makedirs(output_folder, exist_ok=True)
    
    # Build filename
    filename = f"{'_'.join(enzyme_names)}_peptides_Sequential_digest_from_{protein_description}.pdf"
    output_filename = os.path.join(output_folder, filename)
    toc_link = None
    pdf = PDF(toc_link)
    pdf.add_page()
    pdf.set_font("Courier", 'B', size=12)
    pdf.cell(0, 5, f"Enzyme cascade: {' -> '.join(enzyme_names)}\n", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    pdf.set_font("Courier", 'B', size=10)
    pdf.cell(0, 5, f"Number of digestions: {len(operations_log)}\n\n", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(5)
    pdf.cell(0, 5, "Results for each enzyme are on separate pages below", new_x="LMARGIN", new_y="NEXT" ) 
    pdf.ln(5)

    for entry in operations_log:
        entry['link_id'] = pdf.add_link()
    toc_link = pdf.add_link()
    pdf.set_link(toc_link)
    pdf.toc_link = toc_link
    pdf.set_font("Courier", "B", 10)
    pdf.cell(0, 10, "Table of Contents", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)
    pdf.set_font("Courier", size=8)

    for entry in operations_log:
        name = entry.get('enzyme', 'Unknown')
    
        pdf.cell(0, 7, f"Jump to Enzyme: {name}", link=entry['link_id'], new_x="LMARGIN", new_y="NEXT")
    pdf.show_footer = True

    gap = 2
    for digestion in operations_log:
        pdf.add_page()
        pdf.set_link(digestion['link_id'])
        pdf.set_font("Courier", 'B', size=12)
        label_enzyme = "Enzyme: "
        width_enzyme = pdf.get_string_width(label_enzyme) + gap
        pdf.cell(width_enzyme, 5, label_enzyme, new_x="RIGHT", new_y="TOP")
    
        pdf.set_font("Courier", size=12)
        pdf.cell(0, 5, f"{digestion.get('enzyme', 'Unknown')}", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)

        if digestion.get('is_first'):
            cleavage_sites = digestion.get('cleavage_sites', [])
            peptides_list = digestion.get('peptides_list', [])

            pdf.set_font("Courier", 'B', size=8)
            label_quantity = "The quantity of cleavage sites: "
            width_quantity = pdf.get_string_width(label_quantity) + gap
            pdf.cell(width_quantity, 5, label_quantity, new_x="RIGHT", new_y="TOP")
            pdf.set_font("Courier", size=8)
            pdf.cell(0, 5, f"{len(cleavage_sites)}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            pdf.set_font("Courier", 'B', size=8)
            label_number = "The number of peptides: "
            width_number = pdf.get_string_width(label_number) + gap
            pdf.cell(width_number, 5, label_number, new_x="RIGHT", new_y="TOP")
            pdf.set_font("Courier", size=8)
            pdf.cell(0, 5, f"{len(peptides_list)}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            pdf.set_font("Courier", 'B', size=8)
            label_peptides_list = "The peptides list: "
            width_peptides_list = pdf.get_string_width(label_peptides_list) + gap
            pdf.cell(width_peptides_list, 5, label_peptides_list, new_x="RIGHT", new_y="TOP")
            pdf.set_font("Courier", size=8)
            pdf.multi_cell(0, 5, f"{peptides_list}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            pdf.set_font("Courier", 'B', size=8)
            label_cleavage_list = "The cleavage sites list:"
            width_cleavage_list = pdf.get_string_width(label_cleavage_list) + gap
            pdf.cell(width_cleavage_list, 5, label_cleavage_list, new_x="RIGHT", new_y="TOP")
            pdf.set_font("Courier", size=8)
            pdf.multi_cell(0, 5, f"{cleavage_sites}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            pdf.set_font("Courier", 'B', size=8)
            label_relationship = "The site-peptide relationship: "
            width_relationship = pdf.get_string_width(label_relationship) + gap
            pdf.cell(width_relationship, 5, label_relationship, new_x="RIGHT", new_y="TOP")
        
            pdf.set_font("Courier", size=8)
            relationship = dict(zip(peptides_list, cleavage_sites))
            pdf.multi_cell(0, 5, f"{relationship}", new_x="LMARGIN", new_y="NEXT")
            pdf.cell(0, 5, "-" * 160, new_x="LMARGIN", new_y="NEXT")    
             
        else:
            all_peptides = digestion.get('all_peptides', [])

            pdf.set_font("Courier", 'B', size=8)
            label_step_number = "The number of peptides produced in this step: "
            width_step_number = pdf.get_string_width(label_step_number) + gap
            pdf.cell(width_step_number, 5, label_step_number, new_x="RIGHT", new_y="TOP")

            pdf.set_font("Courier", size=8)
            pdf.cell(0, 5, f"{len(all_peptides)}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            pdf.set_font("Courier", 'B', size=8)
            label_aggregated = "All peptides (aggregated) from this step: "
            width_aggregated = pdf.get_string_width(label_aggregated) + gap
            pdf.cell(width_aggregated, 5, label_aggregated, new_x="RIGHT", new_y="TOP")

            pdf.set_font("Courier", size=8)
            pdf.multi_cell(0, 5, f"{all_peptides}", new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)

            per_peptide = digestion.get('per_peptide', [])
            if per_peptide:
                pdf.set_font('Courier', 'B', size = 8)
                pdf.cell(0, 5, "Per-peptide breakdown:\n", new_x="LMARGIN", new_y="NEXT")
                pdf.ln(3)

                for entry in per_peptide:
                    has_digested = len(entry.get('cleavage_sites', [])) > 0
                    pdf.set_font('Courier', 'B', size = 8)
                    pdf.set_text_color(0,0,0)
                    pdf.cell(0, 5, f"Original peptide:", new_x="LMARGIN", new_y="NEXT")

                    pdf.set_font('Courier', size = 8)
                    pdf.multi_cell(0, 5, f"{entry['former_peptide']}", new_x="LMARGIN", new_y="NEXT")

                    pdf.set_text_color(0,0,0)
                    pdf.cell(pdf.get_string_width("Produced peptides: "), 5, "Produced peptides:", new_x="RIGHT", new_y="TOP")
                    if has_digested:
                        pdf.set_text_color(0, 0, 139)
                    pdf.multi_cell(0, 5, f"{entry['produced']}", new_x ="LMARGIN", new_y = "NEXT")

                    pdf.set_text_color(0,0,0)
                    pdf.cell(pdf.get_string_width("Cleavage sites: "), 5, "Cleavage sites: ", new_x="RIGHT", new_y="TOP")
                    if has_digested:
                        pdf.set_text_color(255, 0, 0)
                    pdf.multi_cell(0, 5, f"{entry['cleavage_sites']}", new_x="LMARGIN", new_y="NEXT")

                    pdf.set_text_color(0, 0, 0)
                    pdf.ln(3)
                
            pdf.ln(5)
            pdf.cell(0, 5,"-" * 160 , new_x="LMARGIN", new_y="NEXT")
    print(f"     ✓ Look for results in: results/")
    print(f"        └── {folder_name}/")
    print(f"            └── {filename}\n")
    return pdf.output(output_filename)
    
def main():
   # If no CLI arguments are provided, use the interactive mode
    if len(sys.argv) == 1:
        input_file, digestion_mode, input_enzyme_number = regular_input_storer()
        input_files = [input_file] # Wrap in a list to use the same loop below
    else:
        # Get multiple files from CLI
        input_files, digestion_mode, input_enzyme_number = cli_input_storer()
        
        # Guard rail: If CLI was used but files are missing
        if not input_files:
            print("Error: No input files provided.")
            sys.exit(1)

    os.makedirs("results", exist_ok=True)

    # Now the logic handles both 1 file or 100 files identically
    for fasta_path in input_files:
        print(f"\n--- Processing File: {fasta_path} ---")
        
        if not os.path.exists(fasta_path):
            print(f"File not found: {fasta_path}. Skipping.")
            continue

        proteins = parse_multi_fasta(fasta_path)
        
        for idx, protein_data in enumerate(proteins, 1):
            sequence = protein_data['sequence']
            description = protein_data.get('description') or protein_data['accession_number']
            protein_description = filename_converter(description)
            
            try:
                # The validator executes the digestion and PDF creation
                input_validator(digestion_mode, input_enzyme_number, sequence, protein_description)
                print(f"[{idx}/{len(proteins)}] Processed accession number: {protein_data['accession_number']}")
            except Exception as e:
                print(f"Error processing {protein_data['accession_number']}: {e}")

if __name__ == "__main__":
    main()