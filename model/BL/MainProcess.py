from model.BL.StatisticalResultProcess import StatisticalResultProcess
from model.entity.Combine import Combine
from model.DB.db_model import DB

from model.BL.DuplicateCheck import *
import time

class MainProcess:
    def __init__(self, gene_path, genome_path, create_drop=True):
        self.create_drop = create_drop
        self.gene_path = gene_path
        self.genome_path = genome_path
        if self.create_drop and os.path.exists("results") and os.path.exists("wgs"):
            os.rmdir("results")
            os.rmdir("wgs")
            db = DB()
            db.execute_custom_query("DROP DATABASE `wgs2`;")


    def process(self):
        start_time = time.time()

        gene_sample_path = self.gene_path
        genome_sample_path = self.genome_path
        folder_paths = [gene_sample_path, genome_sample_path]

        combine = Combine()
        combine.create_results_folder()

        # for folder_path in folder_paths:
        #     combine.process_files_to_clean_fasta(folder_path)

        # combine.process_files_to_clean_fasta(gene_sample_path)

        db = DB()
        db.create_initial_tables(folder_paths)

        combine.create_combined_wgs()

        blast = BLAST()
        blast.create_blast_database()

        dc = DuplicateCheck()

        genes_list = db.search_all_genes()
        identity = 85
        coverage = 90

        for gene in genes_list:
            print(f"------------> Process Starting: {gene.name} <------------")
            s1 = time.time()
            blast.blast(gene)
            db.create_result_table(gene.name)
            time.sleep(1)
            db.insert_blast_result(gene.name)
            print(f"-----------> cutoff started")
            db.update_cutoff_column(gene.name, identity, coverage)
            print(f"-----------> duplicate started")
            total_blasts_duplicate = dc.update_duplicate_column(gene.name)
            print(f"-----------> duplicate finished with {total_blasts_duplicate} blasts")
            print(f"-----------> creating result folders for cutoff and duplicate seqs ")
            combine.create_result_folders_with_seqs(gene.name)
            print('Process Duration: {}'.format(time.time() - s1), end="\n\n")


        print("------------> Analysis Started <------------")
        statistical_analysis = StatisticalResultProcess()
        statistical_analysis.analyze_genes()

        db.create_genome_gene_table()

        db.export_table("statistical_result", "statistical_result", "excel", "results/analysis_results")
        db.export_table("genome_gene", "genome_gene", "excel", "results/analysis_results")
        print("analysis exported in excel files.")

        print('Duration: {}'.format(time.time() - start_time))