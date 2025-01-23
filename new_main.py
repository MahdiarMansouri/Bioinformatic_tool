from model.BL.StatisticalResultProcess import StatisticalResultProcess
from model.entity.Combine import Combine
from model.DB.db_model import DB

from model.BL.DuplicateCheck import *
import time

if __name__ == '__main__':
    start_time = time.time()



    gene_sample_path = r'C:\Users\Mahdiar\Desktop\sample_genes'
    genome_sample_path = r'C:\Users\Mahdiar\Desktop\sample_wgs'
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
        # db.create_and_insert_blast_results(gene.name, gene.name)
        print(f"-----------> cutoff started")
        db.update_cutoff_column(gene.name, identity, coverage)
        print(f"-----------> duplicate started")
        total_blasts_duplicate = dc.update_duplicate_column(gene.name)
        print(f"-----------> duplicate finished with {total_blasts_duplicate} blasts")
        print('Process Duration: {}'.format(time.time() - s1), end="\n\n")

    print("------------> Analysis Started <------------")
    statistical_analysis = StatisticalResultProcess()
    statistical_analysis.analyze_genes()

    db.create_genome_gene_table()

    db.export_table("statistical_result", "statistical_result", "excel", "results/analysis_results")
    db.export_table("genome_gene", "genome_gene", "excel", "results/analysis_results")
    print("analysis exported in excel files.")

    print('Duration: {}'.format(time.time() - start_time))
