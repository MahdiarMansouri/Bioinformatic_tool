from model.BL.MainProcess import MainProcess

if __name__ == '__main__':
    gene_sample_path = r'C:\Users\Mahdiar\Desktop\sample_genes'
    genome_sample_path = r'C:\Users\Mahdiar\Desktop\sample_wgs'

    main_process = MainProcess(gene_sample_path, genome_sample_path, create_drop=False)
    main_process.process()
