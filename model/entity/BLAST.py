import subprocess
import os


class BLAST:
    def create_blast_database(self):
        # Create BLAST database
        os.makedirs('wgs', exist_ok=True)
        result = subprocess.run(["makeblastdb", "-in", "results/combined_wgs.fasta", "-dbtype", "nucl", "-out", "wgs/WGS"],
                                capture_output=True,
                                text=True)
        # print("makeblastdb output:", result.stdout)
        # if result.stderr:
        #     print("makeblastdb error:", result.stderr)

    def blast(self, gene):
        # Perform BLAST search
        os.makedirs(f"results/blast_results/{gene.name}", exist_ok=True)
        result = subprocess.run(
            ["blastn", "-query", f"{gene.file_path}", "-db", "wgs/WGS", "-out", fr"results/blast_results/{gene.name}/{gene.name}.csv", "-outfmt",
             "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qframe sframe qseq sseq"],
            capture_output=True, text=True
        )
        # print("blastn output:", result.stdout)
        # if result.stderr:
            # print("blastn error:", result.stderr)

    def blast_for_duplicate(self, seq1_path, seq2_path):
        blast_command = ["blastn", "-query", seq1_path, "-subject", seq2_path, "-outfmt",
                         "10 pident qlen slen qstart qend sstart send"]
        result = subprocess.run(blast_command, capture_output=True, text=True)
        # print("BLAST output:", result.stdout)
        # if result.stderr:
        #     print("BLAST error:", result.stderr)
        blast_output = result.stdout.strip().replace('\n', ',').split(',')
        return blast_output
