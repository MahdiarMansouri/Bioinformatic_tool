import os
import re
import shutil

from model.DB.db_model import DB


class Combine:
    def create_combined_wgs(self):
        self.db = DB()
        genomes = self.db.search_all_genomes()

        genome_dict = {}
        for genome in genomes:
            genome_dict[genome.name] = genome.file_path

        output_file = "results/combined_wgs.fasta"
        with open(output_file, 'w') as outfile:
            for genome_name, genome_path in genome_dict.items():
                with open(genome_path, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            header = line.strip()
                            # Embed genome name in the query_id
                            new_header = f">{genome_name}|{header[1:]}"
                            outfile.write(new_header + '\n')
                        else:
                            outfile.write(line)

    def combine_files(self, directory_path):
        # Define the output file path
        output_file = os.path.join(directory_path, "combined_fasta.fasta")

        # Initialize a list to store file contents
        fasta_contents = []

        # Loop through all files in the directory
        for filename in os.listdir(directory_path):
            # Check if the file is a FASTA file (assuming they have .fasta extension)
            if filename.endswith(".fasta"):
                file_path = os.path.join(directory_path, filename)
                with open(file_path, 'r') as file:
                    # Add a header with the filename to maintain meaningful names
                    fasta_contents.append(f">{filename}")
                    # Read and append the sequences, ensuring no empty lines or malformed entries
                    sequence_lines = [line.strip() for line in file if line.strip() and not line.startswith('>')]
                    fasta_contents.extend(sequence_lines)

        # Write the combined contents to the output file if at least 2 sequences exist
        if len(fasta_contents) >= 4:  # Minimum 2 sequences with headers and sequences
            with open(output_file, 'w') as output:
                output.write('\n'.join(fasta_contents))
            print(f"FASTA files combined into {output_file}")
        else:
            print("Error: A minimum of 2 sequences is required.")

    def process_files_to_clean_fasta(self, folder_path):
        for idx, file_name in enumerate(os.listdir(folder_path)):
            gene_name = file_name.split(".")[0]
            new_file_name = f"{gene_name}_{idx}.fasta"

            # Build the full paths for old and new names
            old_file_path = os.path.join(folder_path, file_name)
            new_file_path = os.path.join(folder_path, new_file_name)

            # Check if the file exists before renaming
            if os.path.exists(old_file_path):
                print(f"Renaming: {old_file_path} to {new_file_path}")
                os.rename(old_file_path, new_file_path)
            else:
                print(f"File not found: {old_file_path}")

            # Open and read the content of the file
            with open(new_file_path, 'r') as file:
                content = file.read()

            # Define the regex pattern to match sequences of 10 or more uppercase English letters
            pattern = r'[A-Z]{10,}'

            # Find all sequences that match the pattern
            matches = re.findall(pattern, content)

            # Write the matches into a new file
            with open(new_file_path, 'w') as output_file:
                for match in matches:
                    output_file.write(match + '\n')

    def create_results_folder(self):
        os.makedirs("results", exist_ok=True)
        os.makedirs("results/analysis_results", exist_ok=True)
        os.makedirs("results/blast_results", exist_ok=True)

    def create_result_folders_with_seqs(self, gene_name):
        for folder_name in ["cutoff", "duplicate"]:
            seq_paths = self.db.search_seq_paths_by_gene_name_with_cutoff_or_duplicate_1(folder_name, gene_name)
            destination_folder = f"results/blast_results/{gene_name}/{folder_name}"
            os.makedirs(destination_folder, exist_ok=True)

            for seq_path in seq_paths:
                shutil.copy(seq_path, destination_folder)
