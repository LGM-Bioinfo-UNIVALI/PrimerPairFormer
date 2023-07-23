import pandas as pd

# Leia o arquivo TSV e armazene-o em um DataFrame
# Supondo que o arquivo TSV esteja em 'caminho/do/arquivo.tsv'
df = pd.read_csv('results/results.tsv', sep='\t')

# Criar uma função para calcular a distância entre dois primers
def calcular_distancia(primer1, primer2):
    return abs(primer1 - primer2)

# Criar um dicionário para armazenar as distâncias entre pares de primers para cada organismo
distancias_por_organismo = {}

# Obter a lista de nomes de primers
primers = df['Primer'].tolist()

pares = []
# Iterar sobre os organismos (colunas) e calcular as distâncias entre os pares de primers
for i, primer_i in enumerate(primers):	
    for j, primer_j in enumerate(primers):
    	distancias = []
    	if primer_i != primer_j:
    		for organismo in df.columns[1:]:
    			distancia = calcular_distancia(df[organismo][i], df[organismo][j])
    			
    			# Verificar se a distância está dentro do intervalo de 500 a 800
    			if 500 <= distancia <= 800:
    				distancias.append(distancia)
    		
    		if len(distancias) >= len(df.columns[1:]) / 2:
    			if (primer_j, primer_i) not in pares:
    				pares.append((primer_i, primer_j))
print(len(pares))
    			
    		

