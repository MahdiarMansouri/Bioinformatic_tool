import os
import shutil
import pandas as pd
import mysql.connector
from ..entity.Gene import *
from ..entity.Genome import WholeGenome


class DB:
    def __init__(self):
        self.db_info = {
            'host': 'localhost',
            'user': 'root',
            'password': 'root123',
            'database': 'wgs2'
        }

    def connect(self):
        try:
            self.mydb = mysql.connector.connect(
                host=self.db_info['host'],
                user=self.db_info['user'],
                passwd=self.db_info['password'],
                database=self.db_info['database']
            )
            self.cursor = self.mydb.cursor()
            # print("Successfully connected to the database.")
        except mysql.connector.Error as err:
            print(f"Error: {err}")

    def disconnect(self, commit=False):
        if commit:
            self.mydb.commit()
        self.cursor.close()
        self.mydb.close()

    def create_initial_tables(self, folder_paths):
        mydb = mysql.connector.connect(
            host=self.db_info['host'],
            user=self.db_info['user'],
            passwd=self.db_info['password'])
        cursor = mydb.cursor()
        cursor.execute("CREATE DATABASE IF NOT EXISTS WGS2")
        cursor.execute("USE WGS2")

        # Create table query
        # todo: change file_name to gene_name in position 2 instead of column 3

        create_genes_table_query = '''
        CREATE TABLE IF NOT EXISTS gene_files (
            id INT AUTO_INCREMENT PRIMARY KEY,
            name NVARCHAR(255) NOT NULL,
            file_path NVARCHAR(255) NOT NULL,
            file_name NVARCHAR(255) NOT NULL,
            UNIQUE(file_path(255))
        )
        '''

        create_genomes_table_query = '''
                    CREATE TABLE IF NOT EXISTS genome_files (
                        id INT AUTO_INCREMENT PRIMARY KEY,
                        name NVARCHAR(255) NOT NULL,
                        file_path NVARCHAR(255) NOT NULL,
                        file_name NVARCHAR(255) NOT NULL,
                        UNIQUE(file_path(255))
                    )
                    '''

        create_statistical_result_table_query = """
                CREATE TABLE IF NOT EXISTS statistical_result (
                    gene_name VARCHAR(100) PRIMARY KEY,
                    gene_presence_count INT,
                    gene_presence_percentage FLOAT,
                    cutoff_count INT,
                    cutoff_percentage FLOAT,
                    duplicate_count INT,
                    duplicate_percentage FLOAT,
                    diversity_count INT,
                    diversity_percentage FLOAT,
                    distinct_gene_presence_count INT

                )
                """
        cursor.execute(create_genes_table_query)
        cursor.execute(create_genomes_table_query)
        cursor.execute(create_statistical_result_table_query)
        mydb.commit()

        cursor = mydb.cursor()
        for idx, name in enumerate(["gene_files", "genome_files"]):
            # Get existing file paths from database
            cursor.execute(f"SELECT file_path FROM {name}")
            existing_files = set(row[0] for row in cursor.fetchall())

            # Get current file details
            current_files = self.get_files(folder_paths[idx])

            # Insert new files into the table
            new_files = [(name, file_path, file_name) for name, file_path, file_name in current_files if
                         file_path not in existing_files]
            if new_files:
                insert_query = f"INSERT INTO {name} (name, file_path, file_name) VALUES (%s, %s, %s)"
                cursor.executemany(insert_query, new_files)
                mydb.commit()

            # Remove paths from database if file no longer exists
            current_file_paths = set(file_path for name, file_path, file_name in current_files)
            files_to_remove = [file_path for file_path in existing_files if file_path not in current_file_paths]
            if files_to_remove:
                delete_query = f"DELETE FROM {name} WHERE file_path = %s"
                cursor.executemany(delete_query, [(file_path,) for file_path in files_to_remove])
                mydb.commit()
        cursor.close()
        mydb.close()

    def get_files(self, folder_path):
        file_details = []
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                name = file.split('.')[0]
                name = name.replace(" ", "")
                file_path = os.path.join(root, file)
                file_details.append((name, file_path, file))
        return file_details

    def table_exists(self, table_name):
        self.connect()
        self.cursor.execute(f"SHOW TABLES LIKE '{table_name}'")
        result = self.cursor.fetchone()
        self.disconnect()
        return result is not None

    def column_exists(self, table_name, column_name):
        self.connect()
        query = f"SHOW COLUMNS FROM {table_name} LIKE '{column_name}'"
        self.cursor.execute(query)
        result = self.cursor.fetchone()
        self.disconnect()
        return result is not None

    def create_and_insert_blast_results(self, table_name, csv_file):
        self.connect()

        # Check if the table already exists
        if self.table_exists(table_name):
            print(f"Table '{table_name}' already exists. Skipping creation and insertion.")
            self.disconnect()
            return

        # Define table columns and types
        columns = '''
            id INT AUTO_INCREMENT PRIMARY KEY,
            query_id NVARCHAR(100),
            genome_name NVARCHAR(100),
            subject_id VARCHAR(100),
            identity FLOAT,
            alignment_length INT,
            mismatches INT,
            gap_opens INT,
            q_start INT,
            q_end INT,
            s_start INT,
            s_end INT,
            evalue FLOAT,
            bit_score FLOAT,
            query_length INT,
            subject_length INT,
            subject_strand NVARCHAR(20),
            query_frame INT,
            sbjct_frame INT,
            qseq_path NVARCHAR(300),
            sseq_path NVARCHAR(300), 
            cutoff TINYINT,
            duplicate TINYINT
        '''

        create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({columns});"
        insert_query = f"""
            INSERT INTO {table_name} (query_id, genome_name, subject_id, identity, alignment_length,
                                      mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score,
                                      query_length, subject_length, subject_strand, query_frame, sbjct_frame, qseq_path, sseq_path)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        self.connect()
        self.cursor.execute(create_table_query)

        # Read the CSV file
        df = pd.read_csv(f'{csv_file}.csv', header=None)

        base_folder_path = r"results/blast_results"

        folder_path = f"{table_name}_seq_folder"
        os.makedirs(os.path.join(base_folder_path, folder_path), exist_ok=True)

        # Iterate through each row in the DataFrame
        for idx, row in df.iterrows():
            query_id = row[1]

            # Split the query_id by '|'
            parts = query_id.split('|')

            if len(parts) == 2:
                genome_name, original_query_id = parts
            elif len(parts) > 2:
                genome_name = parts[0]
                original_query_id = parts[2]
                # If needed, you can handle additional parts here
            else:
                print(f"Unexpected format in query_id: {query_id}")
                continue

            # Include genome_name in the file names
            qseq_path = os.path.join(base_folder_path, folder_path, f"{table_name}_{genome_name}_qseq_{idx}.fasta")
            sseq_path = os.path.join(base_folder_path, folder_path, f"{table_name}_{genome_name}_sseq_{idx}.fasta")

            with open(qseq_path, 'w') as qf:
                qf.write(row[17])
            with open(sseq_path, 'w') as sf:
                sf.write(row[18])

            row_data = (original_query_id, genome_name) + tuple(row[1:17]) + (qseq_path, sseq_path)
            self.cursor.execute(insert_query, row_data)

        self.disconnect(commit=True)

    def create_result_table(self, table_name):
        self.connect()

        # Check if the table already exists
        if self.table_exists(table_name):
            print(f"Table '{table_name}' already exists. Skipping creation and insertion.")
            self.disconnect()
            return

        # Define table columns and types
        columns = """            
            id INT AUTO_INCREMENT PRIMARY KEY,
            query_id NVARCHAR(100),
            genome_name NVARCHAR(100),
            subject_id VARCHAR(100),
            identity FLOAT,
            alignment_length INT,
            mismatches INT,
            gap_opens INT,
            q_start INT,
            q_end INT,
            s_start INT,
            s_end INT,
            evalue FLOAT,
            bit_score FLOAT,
            query_length INT,
            subject_length INT,
            subject_strand NVARCHAR(20),
            query_frame INT,
            sbjct_frame INT,
            qseq_path NVARCHAR(300),
            sseq_path NVARCHAR(300), 
            cutoff TINYINT,
            duplicate TINYINT
        """

        create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({columns});"

        self.connect()
        self.cursor.execute(create_table_query)

    def insert_blast_result(self, gene_name):
        insert_query = f"""
                    INSERT INTO {gene_name} (query_id, genome_name, subject_id, identity, alignment_length,
                                              mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score,
                                              query_length, subject_length, subject_strand, query_frame, sbjct_frame, qseq_path, sseq_path, cutoff, duplicate)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                """

        base_folder_path = r"results/blast_results"
        folder_path = f"{gene_name}_seq_folder"
        full_path = os.path.join(base_folder_path, folder_path)
        os.makedirs(full_path, exist_ok=True)

        df = pd.read_csv(fr'results/blast_results/{gene_name}.csv', header=None)

        self.connect()
        # Iterate through each row in the DataFrame
        for idx, row in df.iterrows():
            query_id = row[1]

            # Split the query_id by '|'
            parts = query_id.split('|')

            if len(parts) == 2:
                genome_name, original_query_id = parts
            elif len(parts) > 2:
                genome_name = parts[0]
                original_query_id = parts[2]
            else:
                print(f"Unexpected format in query_id: {query_id}")
                continue

            # Include genome_name in the file names
            qseq_path = os.path.join(full_path, f"{gene_name}_{genome_name}_qseq_{idx}.fasta")
            sseq_path = os.path.join(full_path, f"{gene_name}_{genome_name}_sseq_{idx}.fasta")

            # Write sequences to the respective files
            with open(qseq_path, 'w') as qf:
                qf.write(row[17])
            with open(sseq_path, 'w') as sf:
                sf.write(row[18])

            # Prepare row data for insertion
            row_data = (
                original_query_id,  # Original query ID
                genome_name,  # Genome name
                row[1],  # Subject ID (could be genome name as well)
                row[2],  # Identity
                row[3],  # Alignment Length
                row[4],  # Mismatches
                row[5],  # Gap Opens
                row[6],  # q_start
                row[7],  # q_end
                row[8],  # s_start
                row[9],  # s_end
                row[10],  # evalue
                row[11],  # bit_score
                row[12],  # query_length
                row[13],  # subject_length
                row[14],  # subject_strand
                row[15],  # query_frame
                row[16],  # sbjct_frame
                qseq_path,  # Path to qseq file
                sseq_path,  # Path to sseq file
                0,
                0
            )

            # Connect to the database and execute the insert
            self.cursor.execute(insert_query, row_data)

        # Commit and disconnect
        self.disconnect(commit=True)

    def execute_custom_query(self, query, commit=False):
        self.connect()
        self.cursor.execute(query)
        if query.strip().upper().startswith("SELECT"):
            result = self.cursor.fetchall()
            self.disconnect()
            return result
        if commit:
            self.mydb.commit()
        self.disconnect()

    def search_result_table_by_name(self, table_name):
        self.connect()
        self.cursor.execute(f"SELECT * FROM {table_name}")
        rows = self.cursor.fetchall()
        self.disconnect()
        return rows

    def add_row(self, table_name, row_data):
        self.connect()
        insert_query = f"""
                INSERT INTO {table_name} (genome_name, query_id, subject_id, identity, alignment_length, 
                mismatches, gap_opens, q_start,q_end, s_start, s_end, evalue, bit_score, query_length, 
                subject_length, subject_strand, query_frame, sbjct_frame, qseq_path, sseq_path)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )
            """
        self.cursor.execute(insert_query, row_data)
        self.disconnect(commit=True)

    def delete_row_from_result_table_by_condition(self, table_name, condition):
        self.connect()
        self.cursor.execute(f" DELETE FROM {table_name} WHERE {condition}")
        self.disconnect(commit=True)

    def update_result_table_row_by_condition(self, table_name, updates, condition):
        self.connect()
        self.cursor.execute(f"UPDATE {table_name} SET {updates} WHERE {condition}")
        self.disconnect(commit=True)

    def search_distinct_gene_presence_count(self, table_name):
        self.connect()
        distinct_gene_presence_query = f"SELECT COUNT(DISTINCT genome_name) FROM {table_name} WHERE cutoff = 1"
        self.cursor.execute(distinct_gene_presence_query)
        distinct_gene_presence_count = self.cursor.fetchone()
        self.disconnect()
        return distinct_gene_presence_count[0]

    def search_gene_presence_count(self, table_name):
        self.connect()
        cutoff_query = f"SELECT COUNT(*) FROM {table_name} WHERE cutoff = 1"
        self.cursor.execute(cutoff_query)
        gene_presence_count = self.cursor.fetchone()
        self.disconnect()
        return gene_presence_count[0]

    def search_duplicate_gene_count(self, table_name):
        self.connect()
        duplicate_query = f"SELECT COUNT(*) FROM {table_name} WHERE duplicate = 1"
        self.cursor.execute(duplicate_query)
        diversity_count = self.cursor.fetchone()
        self.disconnect()
        return diversity_count[0]

    def search_all_genes(self):
        self.connect()
        self.cursor.execute("SELECT * FROM gene_files")
        genes = self.cursor.fetchall()
        genes_list = []
        for gene in genes:
            gene = Gene(*gene)
            genes_list.append(gene)

        self.disconnect()
        return genes_list

    def search_gene_by_name(self, gene_name):
        self.connect()
        self.cursor.execute("SELECT * FROM gene_files WHERE file_name LIKE %s", [f'%{gene_name}%'])
        genes = self.cursor.fetchall()
        genes_list = []
        for gene in genes:
            gene = Gene(*gene)
            genes_list.append(gene)

        self.disconnect()
        return genes_list

    def search_gene_by_id(self, id):
        self.connect()
        self.cursor.execute("SELECT * FROM gene_files WHERE id=%s", [id])
        gene = self.cursor.fetchone()
        gene = Gene(*gene)
        self.disconnect()
        return gene

    def search_all_genomes(self):
        self.connect()
        self.cursor.execute("SELECT * FROM genome_files")
        genomes = self.cursor.fetchall()
        genomes_list = []
        for genome in genomes:
            genome = WholeGenome(*genome)
            genomes_list.append(genome)

        self.disconnect()
        return genomes_list

    def search_genome_by_name(self, genome_name):
        self.connect()
        self.cursor.execute("SELECT * FROM genome_files WHERE file_name LIKE %s", [f'%{genome_name}%'])
        genomes = self.cursor.fetchall()
        genomes_list = []
        for genome in genomes:
            genome = WholeGenome(*genome)
            genomes_list.append(genome)

        self.disconnect()
        return genomes_list

    def search_genome_by_id(self, id):
        self.connect()
        self.cursor.execute("SELECT * FROM genome_files WHERE id=%s", [id])
        genome = self.cursor.fetchone()
        genome = WholeGenome(*genome)
        self.disconnect()
        return genome

    def search_all_result_table_names(self):
        self.connect()
        self.cursor.execute("SHOW TABLES")
        tables = self.cursor.fetchall()
        self.disconnect()

        result_table_names = []
        for table in tables:
            table_name = table[0]
            if table_name not in ["gene_files", "genome_files", "statistical_result", "genome_gene"]:
                result_table_names.append(table_name)

        return result_table_names

    def search_row_counts(self, table_name):
        self.connect()
        self.cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        rows = self.cursor.fetchall()
        self.disconnect()
        return rows[0]

    def search_row_counts_distinct_name(self, table_name):
        self.connect()
        self.cursor.execute(f"SELECT COUNT(DISTINCT name) FROM genome_files")
        rows = self.cursor.fetchall()
        self.disconnect()
        return rows[0]

    def export_table(self, table_name, output_file_name, file_format, base_location):
        self.connect()
        select_query = f"SELECT * FROM {table_name}"

        # Execute the query
        df = pd.read_sql_query(select_query, self.mydb)
        self.disconnect()

        if file_format == 'csv':
            df.to_csv(f"{base_location}/{output_file_name}.csv", index=False)
        elif file_format == 'excel':
            df.to_excel(f"{base_location}/{output_file_name}.xlsx", index=False)
        elif file_format == 'json':
            df.to_json(f"{base_location}/{output_file_name}.json", orient='records')

    def update_cutoff_column(self, table_name, identity, coverage):
        self.connect()
        update_query = f"""UPDATE {table_name} SET cutoff = CASE
                                    WHEN identity < {identity} OR (alignment_length / query_length) * 100 < {coverage} or evalue > 0.05 THEN 0
                                    ELSE 1
                                END
                            """
        self.cursor.execute(update_query)
        print(f"Column 'cutoff' updated in {table_name} table.")
        self.disconnect(commit=True)

    def insert_statistical_result_table(self, statistical_report):
        insert_query = """
                    INSERT INTO statistical_result (gene_name, cutoff_count, cutoff_percentage, duplicate_count, duplicate_percentage, gene_presence_count, gene_presence_percentage, diversity_count, diversity_percentage, distinct_gene_presence_count)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON DUPLICATE KEY UPDATE
                    cutoff_count = VALUES(cutoff_count),
                    cutoff_percentage = VALUES(cutoff_percentage),
                    duplicate_count = VALUES(duplicate_count),
                    duplicate_percentage = VALUES(duplicate_percentage),
                    gene_presence_count = VALUES(gene_presence_count),
                    gene_presence_percentage = VALUES(gene_presence_percentage),
                    diversity_count = VALUES(diversity_count),
                    diversity_percentage = VALUES(diversity_percentage),
                    distinct_gene_presence_count = VALUES(distinct_gene_presence_count)
                    """

        self.connect()
        self.cursor.execute(insert_query, list(vars(statistical_report).values()))
        self.disconnect(commit=True)

    def get_sequences(self, table_name):
        select_query = f"""
            SELECT id, sseq_path, identity, alignment_length, query_length
            FROM {table_name}
            WHERE cutoff = 1
        """
        self.connect()
        self.cursor.execute(select_query)
        sequences = self.cursor.fetchall()
        self.disconnect()
        return sequences

    def update_duplicate_column(self, table_name, id, value):
        self.connect()
        self.cursor.execute(f"UPDATE {table_name} SET duplicate = {value} WHERE id = {id}")
        self.disconnect(commit=True)

    def show_table_contents(self, table_name):
        rows = self.search_result_table_by_name(table_name)
        # Print column headers
        print("Database Contents:")
        print("-------------------")
        # todo: there is no cursor ...
        columns = [desc[0] for desc in self.cursor.description]
        print("\t".join(columns))
        # Print each row
        for row in rows:
            print("\t".join(str(col) for col in row))

    def create_genome_gene_table(self):
        self.connect()
        # Create the genome_gene table if it doesn't exist
        self.cursor.execute("""
                CREATE TABLE IF NOT EXISTS genome_gene (
                    genome_name VARCHAR(100) PRIMARY KEY
                )
                """)

        # Get all gene tables excluding specific tables
        self.cursor.execute("SHOW TABLES")
        tables = self.search_all_result_table_names()

        for table_name in tables:
            self.connect()
            # Add a column for the gene if it doesn't exist
            self.cursor.execute(f"SHOW COLUMNS FROM genome_gene LIKE '{table_name}'")
            if not self.cursor.fetchone():
                self.cursor.execute(f"ALTER TABLE genome_gene ADD COLUMN {table_name} INT DEFAULT 0")

            # Insert genome names into the genome_gene table if they don't exist
            insert_genome_query = """
                    INSERT INTO genome_gene (genome_name)
                    SELECT DISTINCT name FROM genome_files
                    ON DUPLICATE KEY UPDATE genome_name = VALUES(genome_name)
                    """
            self.cursor.execute(insert_genome_query)

            # Update the genome_gene table for genomes that have the gene (cutoff = 1)
            update_query = f"""
                    UPDATE genome_gene gg
                    JOIN {table_name} g ON gg.genome_name = REPLACE(g.genome_name, '.fas', '')
                    SET gg.{table_name} = 1
                    WHERE g.cutoff = 1
                    """
            self.cursor.execute(update_query)
            self.mydb.commit()
        self.disconnect(commit=True)

