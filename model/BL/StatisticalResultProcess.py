import time

from model.DB.db_model import DB
from model.entity.StatisticalReport import StatisticalReport

class StatisticalResultProcess:
    def __init__(self):
        self.db = DB()

    def analyze_genes(self):
        table_names = self.db.search_all_result_table_names()
        print(table_names)

        for table_name in table_names:
            total_blasts = self.db.search_row_counts(table_name)
            total_blasts = total_blasts[0]
            print("total_blasts", total_blasts)
            total_count = self.db.search_row_counts_distinct_name(table_name)
            total_count = total_count[0]
            print("total_count", total_count)

            gene_presence_count = self.db.search_gene_presence_count(table_name)
            print("gene_presence_count", gene_presence_count)
            distinct_gene_presence_count = self.db.search_distinct_gene_presence_count(table_name)
            print("distinct_gene_presence_count", distinct_gene_presence_count)
            cutoff_count = (total_blasts - gene_presence_count)
            if total_blasts == 0:
                continue
            else:
                cutoff_percentage = (cutoff_count / total_blasts) * 100 if total_count else 0
                gene_presence_percentage = (distinct_gene_presence_count / total_count) * 100 if total_count else 0

            diversity_count = self.db.search_duplicate_gene_count(table_name)
            duplicate_count = (gene_presence_count - diversity_count)
            if gene_presence_count == 0:
                continue
            else:
                duplicate_percentage = (duplicate_count / gene_presence_count) * 100 if total_count else 0
                diversity_percentage = (diversity_count / gene_presence_count) * 100 if total_count else 0

            # Insert or update analysis results in gene_analysis table
            statistical_report = StatisticalReport(table_name, cutoff_count, cutoff_percentage, duplicate_count, duplicate_percentage, gene_presence_count,
            gene_presence_percentage, diversity_count, diversity_percentage, distinct_gene_presence_count)

            self.db.insert_statistical_result_table(statistical_report)

        print("Analysis complete.")
