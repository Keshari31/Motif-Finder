"""
This script is a motif finder and sequence analyzer. It provides functionalities to analyze DNA, RNA, and protein sequences,
including searching for specific elements such as TATA and CCAAT boxes, calculating AT and GC content, searching for motif patterns,
and predicting the presence of a promoter in DNA sequences.

The script allows users to interactively input sequences or read sequences from a FASTA file. It provides a menu-driven interface
to select different analysis options and presents the results in a user-friendly format with colored output.

The code is organized into classes for sequence analysis (Sequence, DNAseq, RNAseq) and provides methods for various analysis tasks.

Usage:
    - Run the script and follow the menu prompts to perform sequence analysis tasks.
"""

from colorama import Fore, Style  # For colored output
from prettytable import PrettyTable  # For table structures
import re  # For regular expressions
class Sequence:
    # Define characters for DNA, RNA, and protein alphabets
    DNA_ALPHABET = {'A', 'T', 'C', 'G', 'N', 'X'}
    RNA_ALPHABET = {'A', 'U', 'C', 'G'}
    PROTEIN_ALPHABET = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'}

    # Constructor to initialize a Sequence object with its sequence and type
    def __init__(self, sequence, sequence_type):
        self.sequence = sequence
        self.sequence_type = sequence_type

    # Static method to split a FASTA file into sequences
    @staticmethod
    def fasta_split(fastafile):
        # Initialize an empty dictionary to store sequences
        sequences = {}
        try:
            with open(fastafile, 'r') as file:
                fasta_content = file.read()
                # Split the FASTA content by > character into individual records
                records = fasta_content.split('>')[1:]
                for record in records:  # Iterate through all the records
                    lines = record.split('\n')  # Split each record into lines
                    header = lines[0]
                    sequence = ''.join(lines[1:])  # Concatenate the remaining lines to form the sequence
                    sequences[header] = [sequence]  # Store the header and sequence in the sequences dictionary
        except FileNotFoundError:
            print(f"{Fore.RED}{Style.BRIGHT}Error: File '{fastafile}' not found.{Style.RESET_ALL}")
        except Exception as e:
            print(f"{Fore.RED}{Style.BRIGHT}An error occurred while reading the file '{fastafile}': {e}{Style.RESET_ALL}")
        return sequences

    # Method to determine the type of sequence (DNA, RNA, protein, or unknown)
    def get_sequence_type(self):

        self.sequence_type = ""
        # Check the sequence against DNA, RNA, and protein alphabets
        if all(char.upper() in self.DNA_ALPHABET for char in self.sequence):
            self.sequence_type = 'DNA'
        elif all(char.upper() in self.RNA_ALPHABET for char in self.sequence):
            self.sequence_type = 'RNA'
        elif all(char.upper() in self.PROTEIN_ALPHABET for char in self.sequence):
            self.sequence_type = 'Amino acid'
        else:
            self.sequence_type = 'unknown'
        return self.sequence_type


    # Method to search for a specific element in the sequence (TATA box, CCAAT box, AT content, GC content)
    def search_element(self, element, max_mismatches):
        if self.sequence_type == 'DNA':
            # If the sequence is DNA, create a DNAseq object for DNA-specific operations
            dna_seq = DNAseq(self.sequence, self.sequence_type)
            # Determine the element type and call the corresponding method in the DNAseq class
            if element.upper() == 'TATA':
                return dna_seq.find_tata_box(int(max_mismatches))
            elif element.upper() == 'CCAAT':
                return dna_seq.find_ccaat_box(int(max_mismatches))
            elif element.upper() == 'AT CONTENT':
                return dna_seq.AT_content()
            elif element.upper() == 'GC CONTENT':
                return dna_seq.GC_content()
            else:
                # Handle invalid element type
                return f"{Fore.RED}Invalid element type.{Style.RESET_ALL}"
        elif self.sequence_type == 'RNA':
            # If the sequence is RNA, create an RNAseq object for RNA-specific operations
            rna_seq = RNAseq(self.sequence, self.sequence_type)
            # Determine the element type and call the corresponding method in the RNAseq class
            if element.upper() == 'AT CONTENT':
                return rna_seq.AT_content()
            elif element.upper() == 'GC CONTENT':
                return rna_seq.GC_content_RNA()
            else:
                # Handle invalid element type
                return f"{Fore.RED}Invalid element type.{Style.RESET_ALL}"
        else:
            # Handle unsupported sequence type
            return f"{Fore.RED}Unsupported sequence type.{Style.RESET_ALL}"

    # Method to search for a specific pattern (motif) in the sequence using regular expressions
    def search_sequence_pattern(self, pattern):
        # Use re.finditer() to find all occurrences of the pattern in the sequence
        matches = re.finditer(pattern.upper(), self.sequence)
        # Initialize empty lists to store motif positions and patterns
        motif_positions = []
        motif_patterns = []
        # Check if any matches are found
        if matches:
            # Iterate through each match found
            for match in matches:
                #If pattern is found add the position and pattern to the list
                position = match.start()
                pattern = match.group()
                motif_positions.append(position)
                motif_patterns.append(pattern)
        # Return the lists containing motif positions and patterns
        return motif_positions, motif_patterns

# Subclass for DNA sequences with additional methods for DNA-specific analysis
class DNAseq(Sequence):
    def __init__(self, sequence, sequence_type):
        super().__init__(sequence, sequence_type)

    # Method to calculate AT content in the DNA sequence
    def AT_content(self):
        # Initialize a counter for counting occurrences of A and T nucleotides
        AT_count = 0
        # Get the length of the sequence
        length = len(self.sequence)
        # Iterate through each element (nucleotide) in the sequence
        for element in self.sequence:
            # Check if the element is either "A" or "T"
            if element == "A" or element == "T":
                # Increment the AT count if the element is "A" or "T"
                AT_count += 1
        # Increment the AT count if the element is "A" or "T"
        self.AT_content = AT_count / length
        return self.AT_content

    # Method to calculate GC content in the DNA sequence
    def GC_content(self):
        GC_count = 0
        length = len(self.sequence)
        for element in self.sequence:
            if element == "G" or element == "C":
                GC_count += 1

        self.GC_content = GC_count / length
        return self.GC_content

    # Method to find the TATA box motif in the DNA sequence
    def find_tata_box(self, max_mismatches=0):
        tata_motif = "TATAAA"  # TATA box motif
        motif_positions = []
        motif = []
        seq_len = len(self.sequence)  # Get the length of the DNA sequence
        motif_len = len(tata_motif)  # Get the length of the TATA box motif

        # Iterate over possible positions in the DNA sequence for the motif
        for i in range(seq_len - motif_len + 1):
            # Count mismatches between the sequence and the TATA box motif at the current position
            mismatches = sum(self.sequence[i + j] != tata_motif[j] for j in range(motif_len))
            # Check if the number of mismatches is within the allowed limit
            if mismatches <= max_mismatches:
                # If so, add the position of the motif to the list of motif positions
                motif_positions.append(i)
                motif.append(self.sequence[i:i + motif_len])
        # Return the list of motif positions and motif patterns
        return motif_positions, motif

    # Method to find the CCAAT box motif in the DNA sequence
    def find_ccaat_box(self, max_mismatches=0):
        ccaat_motif = "GGCCAAT"
        motif_positions=[]
        motif = []
        seq_len = len(self.sequence)
        motif_len = len(ccaat_motif)

        for i in range(seq_len - motif_len + 1):
            mismatches = sum(self.sequence[i + j] != ccaat_motif[j] for j in range(motif_len))
            if mismatches <= max_mismatches:
                motif_positions.append(i)
                motif.append(self.sequence[i:i + motif_len])

        return motif_positions, motif

    # Method to predict the presence of a promoter in the DNA sequence
    def is_promoter_present(self, max_mismatches=1):
        motif_gap = (40, 60)
        gc_threshold = 0.3
        # Find TATA box
        tata_positions, tata_motif = self.find_tata_box(max_mismatches)

        # Find CCAAT box
        ccaat_positions, _ = self.find_ccaat_box(max_mismatches)

        # Find GC content
        gc_content = self.GC_content()
        tata = []
        ccaat = []

        # Check if both motifs are found
        if not tata_positions or not ccaat_positions:
            return False, None, None, None

        for tata_position in tata_positions:
            for ccaat_position in ccaat_positions:
                # Calculate the distance between TATA and CCAAT boxes
                distance = ccaat_position - tata_position
                # Check if the distance falls within the acceptable range
                if motif_gap[0] <= distance <= motif_gap[1]:
                    # If so, store the positions of TATA and CCAAT boxes
                    tata = tata_position
                    ccaat = ccaat_position
                    # Check if the GC content exceeds the threshold
                    if gc_threshold <= gc_content:
                        # If all conditions are met, return True indicating promoter presence
                        return True, tata, ccaat, gc_content
        # If no suitable promoter region is found, return False
        return False, tata, ccaat, gc_content

# Subclass for RNA sequences with additional methods for RNA-specific analysis
class RNAseq(Sequence):
    def __init__(self, sequence, sequence_type):
        super().__init__(sequence, sequence_type)

    # Method to calculate the AU content of the RNA sequence
    def AT_content(self):
        # Initialize a counter for counting occurrences of A and U nucleotides
        AU_count = 0
        length = len(self.sequence)
        for element in self.sequence:
            # Check if the element is either "A" or "U"
            if element == "A" or element == "U":
                # Increment the AU count if the element is "A" or "U"
                AU_count += 1

        self.AU_content = AU_count / length
        return self.AU_content

    # Method to calculate the GC content of the RNA sequence
    def GC_content_RNA(self):
        GC_count = 0
        length = len(self.sequence)
        for element in self.sequence:
            # Check if the element is either "G" or "C"
            if element == "G" or element == "C":
                # Increment the GC count if the element is "G" or "C"
                GC_count += 1

        self.GC_content = GC_count / length
        return self.GC_content

def seq_element():
    # Prompt the user to input the path to the FASTA file
    fasta_file = input("Enter the path to the FASTA file: ")
    # Prompt the user to input the element type to search for
    element_type = input("Enter the element type (TATA, CCAAT, AT CONTENT, GC CONTENT): ")

    # Split the FASTA file into individual sequences
    sequences = Sequence.fasta_split(fasta_file)

    # Initialize a PrettyTable instance with column names
    table = PrettyTable(["Sequence ID", "Position", "Motif"])

    # Check the input element type and perform corresponding operations
    if element_type.upper() == "TATA":
        max_mismatches_tata = None
        while max_mismatches_tata is None:
            max_mismatches_input = input("Enter maximum number of mismatches TATA box exist/default=0: ")
            if max_mismatches_input.isdigit():
                max_mismatches_tata = int(max_mismatches_input)
            else:
                print(f"{Fore.RED}{Style.BRIGHT}Please enter a valid integer for maximum mismatches.{Style.RESET_ALL}")
        print(f"{Fore.CYAN}{Style.BRIGHT}TATA BOX MOTIF{Style.RESET_ALL}")
        for header, sequence in sequences.items():
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            if seq_type == "DNA":
                # Search for the TATA box motif in the sequence
                motif_positions, motif_with_mismatches = seq.search_element(element_type, max_mismatches_tata)
                # If TATA box motif is found, adds its positions to table
                if motif_positions:
                    for i, motif in zip(motif_positions, motif_with_mismatches):
                        table.add_row([header, i, motif])
                else:
                    table.add_row([header, "N/A", "TATA box motif not found"])
            else:
                table.add_row([header, "N/A", f"{Fore.RED}{Style.BRIGHT}Not a DNA sequence{Style.RESET_ALL}"])
        # Print the table
        print(table)

    elif element_type.upper() == "CCAAT":
        max_mismatches_ccaat = None
        while max_mismatches_ccaat is None:
            max_mismatches_input = input("Enter maximum number of mismatches CCAAT box exist/default=0: ")
            if max_mismatches_input.isdigit():
                max_mismatches_ccaat = int(max_mismatches_input)
            else:
                print(f"{Fore.RED}{Style.BRIGHT}Please enter a valid integer for maximum mismatches.{Style.RESET_ALL}")
        print(f"{Fore.CYAN}{Style.BRIGHT}CCAAT BOX MOTIF{Style.RESET_ALL}")
        for header, sequence in sequences.items():
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            if seq_type == "DNA":
                # Search for the CCAAT box motif in the sequence
                motif_positions, motif_with_mismatches = seq.search_element(element_type, max_mismatches_ccaat)
                # If CCAAT box motif is found, add into the table
                if motif_positions:
                    for i, motif in zip(motif_positions, motif_with_mismatches):
                        table.add_row([header, i, motif])
                else:
                    # If CCAAT box motif is not found, print a message
                    table.add_row([header, "N/A", "CCAAT box motif not found"])
            else:
                # If the sequence is not DNA, print a message
                table.add_row([header, "N/A", f"{Fore.RED}{Style.BRIGHT}Not a DNA sequence{Style.RESET_ALL}"])
        # Print the table
        print(table)


    elif element_type.upper() == "AT CONTENT":
        print(f"{Fore.CYAN}{Style.BRIGHT}AT CONTENT{Style.RESET_ALL}")
        # Initialize a new PrettyTable instance with column names
        table2 = PrettyTable(["Sequence ID", "AT CONTENT"])
        for header, sequence in sequences.items():
            # Create a Sequence object with the current sequence
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            # Calculate the AT content of the sequence
            AT_content = seq.search_element(element_type, "")
            table2.add_row([header, AT_content])
        print(table2)

    elif element_type.upper() == "GC CONTENT":
        print(f"{Fore.CYAN}{Style.BRIGHT}GC CONTENT{Style.RESET_ALL}")
        # Initialize a new PrettyTable instance with column names
        table3 = PrettyTable(["Sequence ID", "GC CONTENT"])
        for header, sequence in sequences.items():
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            GC_content = seq.search_element(element_type, "")
            table3.add_row([header, GC_content])
        print(table3)
    else:
        # If the input element type is incorrect, print an error message
        print(f"{Fore.RED}{Style.BRIGHT}Incorrect Element Type!{Style.RESET_ALL}")


    # Prompt the user whether to search for another element
    another_search = input("Do you want to search for another element? (yes/no): ")
    if another_search.lower() == "yes":
        seq_element()  # Recursive call to main()

def motifPattern():
    fasta_file = input("Enter the path to the FASTA file: ")
    pattern = input("Enter the motif pattern to search as a regular expression: ")
    # Split the FASTA file into sequences
    sequences = Sequence.fasta_split(fasta_file)

    print(f"{Fore.CYAN}{Style.BRIGHT}MOTIF PATTERN{Style.RESET_ALL}")
    # Iterate through each sequence in the FASTA file
    for header, sequence in sequences.items():
        # Create a Sequence object with the current sequence
        seq = Sequence(sequence[0], "")
        seq_type = seq.get_sequence_type()
        print(f"{Fore.GREEN}{Style.BRIGHT}Sequence ID: {header}{Style.RESET_ALL}")
        # Search for the motif pattern in the sequence
        motif_positions, motif_pattern = seq.search_sequence_pattern(pattern)
        if motif_positions:
            print("Motif found at positions:", motif_positions)
            # Print each motif with its position
            for i, motif in zip(motif_positions, motif_pattern):
                print(f"Position {i}: {motif}")
        else:
            print(f"{Fore.RED}{Style.BRIGHT}Motif not found in the given {seq_type} sequence.{Style.RESET_ALL}")

    another_search = input("Do you want to search for another motif pattern? (yes/no): ")
    if another_search.lower() == "yes":
        motifPattern()  # Recursive call to motif_pattern() function

def predict_promoter():
    print(f"{Fore.CYAN}{Style.BRIGHT}PROMOTER PREDICTION{Style.RESET_ALL}")

    # Prompt the user to choose between sequence input and FASTA file input
    input_choice = input("Enter '1' for sequence input or '2' for FASTA file input: ")

    if input_choice == '1':
        # Input sequence directly
        sequence = input("Enter DNA sequence for promoter prediction: ")
        seq_type = Sequence(sequence, "").get_sequence_type()
        if seq_type == "DNA":
            dna_seq = DNAseq(sequence, "DNA")
            promoter_presence, tata_positions, ccaat_positions, gc_content = dna_seq.is_promoter_present(max_mismatches=1)
            if promoter_presence:
                print("Promoter Presence: True")
                print(f"TATA positions in sequence: {tata_positions}, CCAAT positions in sequence: {ccaat_positions}, GC content: {gc_content}")
            else:
                print("Promoter Presence: False")
        else:
            print(f"{Fore.RED}{Style.BRIGHT}Please provide a DNA sequence.{Style.RESET_ALL}")
    elif input_choice == '2':
        # Input FASTA file
        fasta_file = input("Enter the path to the FASTA file: ")
        sequences = Sequence.fasta_split(fasta_file)
        for header, sequence in sequences.items():
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            if seq_type == "DNA":
                dna_seq = DNAseq(sequence[0], "DNA")
                promoter_presence, tata_positions, ccaat_positions, gc_content = dna_seq.is_promoter_present(max_mismatches=1)
                if promoter_presence:
                    print(f"{Fore.GREEN}{Style.BRIGHT}Promoter Presence for sequence {header}:{Style.RESET_ALL} True")
                    print(f"TATA positions in sequence: {tata_positions}, CCAAT positions in sequence: {ccaat_positions}, GC content: {gc_content}")
                else:
                    print(f"Promoter Presence for sequence {header}: False")
            else:
                print(f"{Fore.RED}{Style.BRIGHT}Sequence {header} is not a DNA sequence.{Style.RESET_ALL}")
    else:
        print(f"{Fore.RED}{Style.BRIGHT}Invalid input. Please enter '1' for sequence input or '2' for FASTA file input.{Style.RESET_ALL}")

    # Prompt the user whether to predict the promoter of another sequence
    another_search = input("Do you want to predict the promoter of another sequence? (yes/no): ")
    if another_search.lower() == "yes":
        predict_promoter()  # Recursive call to main3()

def sequence_type():
    print(f"{Fore.CYAN}{Style.BRIGHT}SEQUENCE TYPE{Style.RESET_ALL}")
    # Prompt the user to choose between sequence input and FASTA file input
    input_choice = input("Enter '1' for sequence input or '2' for FASTA file input: ")

    if input_choice == '1':
        # If sequence input is chosen, prompt for the sequence
        sequence = input("Enter sequence: ")
        # Create a Sequence object and get its type
        seq = Sequence(sequence, "")
        seq_type = seq.get_sequence_type()
        print("Sequence type:", seq_type)
    elif input_choice == '2':
        # If FASTA file input is chosen, prompt for the file path
        fasta_file = input("Enter the path to the FASTA file: ")
        # Read the FASTA file and extract the sequence
        sequences = Sequence.fasta_split(fasta_file)
        # Iterate over each sequence and print its ID and type
        for header, sequence in sequences.items():
            seq = Sequence(sequence[0], "")
            seq_type = seq.get_sequence_type()
            print(f"Sequence ID: {header:<30}Sequence type: {seq_type}")

    else:
        # If an invalid input is provided, print an error message
        print(f"{Fore.RED}{Style.BRIGHT}Invalid input. Please enter '1' for sequence input or '2' for FASTA file input.{Style.RESET_ALL}")


if __name__ == '__main__':
    # Main loop to display the menu and handle user choices
    while True:
        print(f"{Fore.YELLOW}{Style.BRIGHT}---------------------------------------------------------------{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}{Style.BRIGHT}                          MOTIF FINDER                         {Style.RESET_ALL}")
        print(f"{Fore.YELLOW}{Style.BRIGHT}---------------------------------------------------------------{Style.RESET_ALL}")
        print(f"Please select an option from menu \n{Fore.GREEN}{Style.BRIGHT}1.Sequence Element(TATA/CCAAT Box, AT Content, GC Content)\n2.MOTIF PATTERN \n3.PROMOTER PREDICTION\n4.SEQUENCE TYPE{Style.RESET_ALL}")
        choice = input("Enter the number for your choice (1/2/3/4): ")

        if choice == '1':
            seq_element()
        elif choice == '2':
            motifPattern()
        elif choice == '3':
            predict_promoter()
        elif choice == '4':
            sequence_type()
        else:
            print(f"{Fore.RED}{Style.BRIGHT}Incorrect choice{Style.RESET_ALL}")

        # Ask the user if they want to continue or exit
        another_search = input("Do you want to go back Main Menu? (yes/no): ")
        if another_search.lower() != "yes":
            break  # Exit the loop if the user doesn't want to search again






























#tata{2}
#CGCTGGGGCGCATTATAAAGCACCTCCTCGCCTCCTCGCAGGGCGGTGGGGCGCAGCGGTTCTATCCCCCTCCCCGAGGCGGGGAAGGCCAATAGGCCAATAGGTCTGGAAAGAAACGTGGGTTCGAGGCGGAGAGGAAAAGCGGACCCACCTGCCAGGCTGCGCGGGGAGGCTGGTCCCGGGCTGGGCAGGCGGGCTGGCCTCGCGCCCTCGAGGCACCCGGCGGCGCTGGCTGTGCGGAGGGGCGCCGGCGCGGCCGTATTTGTACCCGCGGGCCCTCACATGGTCTGATCTCTAGATAGCCGCCGCCAAAGAGCTCTTGAAGAATTTTTGCGTCACTTTGAGGCGAATAAACTTAATGCTTCCCCGCGGCCGCGGCTCCGCGCTCCCGCTGGATGGGGTTGCGCTCGCCAGGGAGGGGCCGCGCTACGGGGCGGGGTGCGCGCCCGACCCCAGAGCCAGGAGGGGAGGGACCCCCGACACACACACACGCTCGCAGGGAGGAGCGGAGCGCGGAGCAGCCGGCAGGGCAGGGCCGAGCGAGGAATCCTGACTTTCCTGCTCTTTAACTTTGCGGGAGGGGGGAGGCTGCACAAAGGAGCAGGTGTGCGCCTCCCCTGCCGCCCCGCGCCCACAGGACGGCACACAGGGGTCTGCTCGTGCCGCTTTCTCTTGACCTCGGACACCTCCCCATCGTCCCCATCCCGAAAGCTGCTCTCTCTTTCTTCCTGGGACCTAGCAGTACTCGATACCCATGCGTTTATCTGACTGACCGGCTGGGATGTGAGCTCAGCTTAGAGGCCCTTGTCTAATTCCTCTCCGTCCCCAGAGTCTCACAATAAAGACTTCCACAATGAAGGAATGAAGAACTGGTTCTCATTTTGTTAACCAGCTACAGTTCAAAGCTGCGTCAAAAAGGAGCACAGTGCTTTGGGAGGCTGAGGCGGGAGGATTGCTTGAGCCCAGGAATTCAAGACCAGCCTGGGCAACACAGCAAGACCCCATCTCTAGAAAAAAAATTTTTTTAATTAGCCGGGCGTGGTTGCTCGCGCTTGTGTTCCCAGCTACTCAGGAGGCTGAGGTGGAGGATCAGCTGAGCCCGGGAGGTCGAGGCTGCAGTGAGCCGTGATCGAACCACTGCACTCCAACCTGGTGACAATTTGAATTATCTGAAAGGGCCCAGCAGGAGCCTATTGTTTTGAACCATGTGCATGTAGTTTTGATAATTTTTTCAAAAAGTTTAAAAATTGATTTAATGTTAATATCCTTTTTGCTAAAATTATTTAAAACTTTAAAAAGGGCACATGCAGTACTTGCTTTGGGTTGTCCCTGAGGCGGAAAGGAGAGGGCTCTACCTCTCCCCAAGAGGAGGCCAGGAACGATAGGAAGTAGAAGACCGAAAGGAAATAGCAGTGACAAGTTTGCAGCTCCAGAGAAGCCACCGCCCCTTGTACTTGGAGGAACTGACCCCTGAAAACTGTGCGGCCGGTTGGGCTGAGCGTCTAGAGGGACTGAGCTGGACAACCACGGGCAAGCGAGGGCAGCTCCCAGCGGGTGGAGTCCGCGCGGGATTCTGGTGCCACCTAGACGCCAGGGCGGGGACCGCAAGGTGGGCGGGAGGCTTGGAGGCCGGGATGCGGGGGAATACTGGTAGGGTGCAAGGAGAATGCTGGAGGGGTGCAGGGGGGATGCCGGGGGTGCATGGGGGGATGCTGGGGGGTGCAGGGGGGATACTGCGAGGGGTGCAGGGGGGATAATGGGGGTTGCAGGGGAGATCCTGGGAGAGGTGCAGGGGGATGCTGGAAGGGCTGCAGGGGGGATGCTGGGGGTGCAGGGGAGATGCTGGGGGGGCTGCAGGGGGGATGCTGGGGGTGCAGGGGGGATGCCGCGAGGGGTGAAGGGGGGGATAATGGGGGATGCAGGAGCGATCCTAGGAGGGGTACAGGGAGTTGCTGGGAGGTTGCAGGGGGGATGCTGGGAGGGCTGCAGTAGGGACATTAGGGGTTAGGTGCTGGCATTCAGAGCGGCAACGCCAGGCGGGAGAATCCCCAACAAGCTTG
#AGGAACTATGACCTCGACTATGACTCGGTGCAGCCGTATTTCTACTGCGACGAAGAAGAGAACTTCTACCAGCAGCAGCAGCAGAGCGAGCTGCAGCCGCCGGCGCCCAGCGAAGATATCTGGAAGAAATTCGAGCTGCTGCCCACCCCACCCCTGTCCCCGAGCCGCCGCTCCGGGCTCTGCTCGCCCCCGTATGTTACGGTCGCATCTTTCTCCCCCCCGGGAGATGATGATGGTGGCGGCGGCAGCTTCTCCACAGCCGACCAGCTGGAGATGGTGACCGAGCTGTTGGGAGGAGATATGGTGAACCAGAGCTTCATCTGCGACCCGGACGACGAGACCTTCATCAAAAACATCATCATCCAGGATTGTATGTGGAGCGGCTTCTCGGCTGCTGCCAAGCTCGTCTCAGAGAAGCTGGCCTCGTACCAGGCTGCGCGCAAAGACAGCGGCAGCCCGAGCCCCGCCCGCGGGCACGGCGGCTGCTCCTCCTCCAGCCTCTACCTACAGGACCTGAGCGCTGCCGCGTCGGAGTGCATCGACCCCTCGGTGGTCTTCCCCTACCCGCTCAACGACAGCAGCTCGCCCAAGCCCTGCGCCTCGCCCGACTCCAGCGCCTTCTCTCCATCCTCGGACTCTCTGCTCTCCTCCACGGAGTCCTCCCCGCGGGCCAGCCCTGAGCCCT



