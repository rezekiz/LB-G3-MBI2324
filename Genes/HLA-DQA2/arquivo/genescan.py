# Função de pesquisa bibliográfica

# ATENÇÃO À UTILIZAÇÃO DESTA FUNÇÃO
# DEMASIADAS QUERIES SEGUIDAS PODEM LEVAR A BLOQUEIO POR PARTE DO SERVIDOR
# EM PRINCÍPIO O MÓDULO ENTREZ ESTÁ PREPARADO PARA LIDAR COM ISSO MAS CONVÉM TER SEMPRE CUIDADO

def lit_search(email, term,db = 'pubmed',save = 'N',ext='.txt', amplitude = 10, displ = 'n'):
    import os

    # TODO Adicionar indicações para fazer o print e elencar a lista de resultados
    # TODO Adicionar documentação


    # Importamos o módulo 
    from Bio import Entrez

    # Definimos o email para os fins de pesquisa
    Entrez.email = email
    
    

    # Secção egquery que pode ser skipable

    handle = Entrez.egquery(term=term)
    egq_res = Entrez.read(handle)

    for _ in egq_res['eGQueryResult']:
        if db in _.values():
            print('Resultados em "pubmed":',_['Count'])

    # Sacar as IDs
    handle      = Entrez.esearch(db=db,term=term, retmax = amplitude)
    esearch_res = Entrez.read(handle)

    lista_ids = esearch_res['IdList']
    print('Top {} resultados -> {}'.format(amplitude,lista_ids)) # Converter para uma lista de títulos mais à frente
    
    

    # Sacar os abstracts
    handle = Entrez.efetch(db=db,id=','.join(lista_ids),rettype='abstract',retmode='text')
    fetch_res = handle.read()
    
    # Usar com print para já, para não estar a bagunçar demasiado.
    # Perceber como se pode implementar a escrita de ficheiros para ser mais fácil de mastigar os resultados.

    # NÃO FUNCIONA NO GITHUB / LIVECODING
    # ADICIONAR PARAMETRO save = 'Y' SE SE PRETENDER GUARDAR FICHEIRO 
    
   

    # Verifica a condição de gravação de output como ficheiro de texto
    if save.upper() == 'Y':

        # Verifica se existe pasta para output e cria-a caso não exista
        cwd = os.getcwd()
        outputdir = os.path.join(cwd,'lit-search-output')
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)


        # Muda o working directory para a pasta de output para gerar o ficheiro
        os.chdir(outputdir)

        # Gera o ficheiro           
        with open(term+'on-'+db+ext , 'w', encoding='utf-8') as _:
            _.write(fetch_res)

        # Retorna ao working directory anterior
        os.chdir(cwd)

    handle.close()

    if displ.upper() == 'Y':
       print(fetch_res,sep='\n')
       return fetch_res
    else:
        return fetch_res

# Exemplo lit_search(email,'hla-dqa2',save='y',amplitude=1,db='nucleotide')



# Criação de ficheiros para análise

def guardar_ficheiro(email, basedados, id_gene, nome_ficheiro):
    
    import os
    from Bio import SeqIO
    from Bio import Entrez

    Entrez.email = email
    
    filename = nome_ficheiro

    if not os.path.isfile(filename):
        net_handle = Entrez.efetch( db = basedados, id = id_gene, rettype="gb", retmode="text")
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        diretoria = os.getcwd()
        print("O ficheiro foi guardado com sucesso na seguinte diretoria:", diretoria)

    return 

# Exemplo guardar_ficheiro("oliveira1mariana@hotmail.com", "nucleotide","NC_000006", "NC_000006.gbk")

# Análise do ficheiro Genbank

def parsing(nome_ficheiro):

    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")

    print("ID do gene:", record.id)

    print("Nome do gene:", record.name)

    print("Descrição do gene", record.description)

    print("Comprimento da sequência:", len(record.seq), "bp")

    return

# Exemplo parsing("HLA-DQA1.gbk")

# Verificação de anotações correspondentes aos genes de interesse

def anot(nome_ficheiro):
    
    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")

    print("Quantidade de anotações:", len(record.annotations))

    print()
    
    print("Lista de anotações:")
    
    for anotacao in record.annotations:
        print(anotacao, "->", record.annotations[anotacao])

    return   

# Exemplo anot("HLA-DQA1.gbk")

# Informação complementar de features e qualifiers

def features_qualifiers(nome_ficheiro):
    
    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")
    
    print("Quantidade de features:", len(record.features))

    for feature in record.features:
        print(feature)
    return

# Exemplo features_qualifiers("HLA-DQA1.gbk")


def analise_ncbi(email, db, id_gene):

    filename = id_gene+'.gbk'

    guardar_ficheiro(email,db,id_gene,filename)
    parsing(filename)
    anot(filename)
    features_qualifiers(filename)




