import json

class StatisticalReport:
    def __init__(self, gene_name, gene_presence_count, gene_presence_percentage, cutoff_count, cutoff_percentage,
                 duplicate_count, duplicate_percentage, diversity_count, diversity_percentage,
                 distinct_gene_presence_count):
        self.gene_name = gene_name
        self.gene_presence_count = gene_presence_count
        self.gene_presence_percentage = gene_presence_percentage
        self.cutoff_count = cutoff_count
        self.cutoff_percentage = cutoff_percentage
        self.duplicate_count = duplicate_count
        self.duplicate_percentage = duplicate_percentage
        self.diversity_count = diversity_count
        self.diversity_percentage = diversity_percentage
        self.distinct_gene_presence_count = distinct_gene_presence_count

    def __repr__(self):
        return json.dumps(self.__dict__)
