import os
from model.DB.db_model import DB
from model.entity.BLAST import BLAST


class DuplicateCheck:
    def __init__(self):
        self.db = DB()

    def blast_sequences(self, seq1_path, seq2_path):
        # Check if sequence files exist and have content
        if not os.path.isfile(seq1_path) or not os.path.isfile(seq2_path):
            print(f"Sequence file missing: {seq1_path} or {seq2_path}")
            return []

        if os.path.getsize(seq1_path) == 0 or os.path.getsize(seq2_path) == 0:
            print(f"Sequence file is empty: {seq1_path} or {seq2_path}")
            return []

        blast = BLAST()
        return blast.blast_for_duplicate(seq1_path, seq2_path)

    def update_duplicate_column(self, table_name):
        total_blasts = 0
        sequences = self.db.get_sequences(table_name)
        permission_sequences = [[id, sseq_path, 1] for id, sseq_path, identity, alignment_length, query_length in
                                sequences]

        seq_count = len(permission_sequences)
        for i in range(seq_count):
            if permission_sequences[i][2] == 1:
                value = 1
                for j in range(i + 1, seq_count):
                    id1, seq1_path, permission1 = permission_sequences[i]
                    id2, seq2_path, permission2 = permission_sequences[j]
                    blast_output = self.blast_sequences(seq1_path, seq2_path)
                    total_blasts += 1

                    # Ensure the BLAST output is valid
                    if len(blast_output) < 7:
                        print(f"Invalid BLAST output for {seq1_path} vs {seq2_path}: {blast_output}")
                        continue

                    try:
                        identity = float(blast_output[0])
                        qlen = int(blast_output[1])
                        slen = int(blast_output[2])
                        qstart = int(blast_output[3])
                        qend = int(blast_output[4])
                        sstart = int(blast_output[5])
                        send = int(blast_output[6])
                    except ValueError as e:
                        print(f"Error converting BLAST output: {e}")
                        continue

                    q_coverage = ((qend - qstart + 1) / qlen) * 100
                    s_coverage = ((send - sstart + 1) / slen) * 100

                    if identity < 100:
                        self.db.update_duplicate_column(table_name, id1, 1)
                    else:
                        if q_coverage >= s_coverage:
                            self.db.update_duplicate_column(table_name, id2, 1)
                            self.db.update_duplicate_column(table_name, id1, 0)
                            value = 0
                        else:
                            self.db.update_duplicate_column(table_name, id1, 1)
                            self.db.update_duplicate_column(table_name, id2, 0)
                            permission_sequences[j][2] = 0
                    if value == 0:
                        break
        return total_blasts

