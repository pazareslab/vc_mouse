import pandas as pd
import sys

# Leer argumentos
input_file = sys.argv[1]
gene_file = sys.argv[2]
output_file = sys.argv[3]

# Cargar archivo de variantes anotadas
vcf_df = pd.read_csv(input_file, sep="\t")

# Cargar genes tumorales desde la columna "Gene" del archivo Intogen.txt
genes_df = pd.read_csv(gene_file, sep="\t")
genes = genes_df["Gene"].str.upper().tolist()

# Filtrar por la columna "SYMBOL"
filtered_df = vcf_df[vcf_df["SYMBOL"].str.upper().isin(genes)]

# Guardar salida
filtered_df.to_csv(output_file, sep="\t", index=False)
