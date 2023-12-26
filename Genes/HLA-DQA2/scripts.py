# ATENÇÃO À UTILIZAÇÃO DESTA FUNÇÃO
# DEMASIADAS QUERIES SEGUIDAS PODEM LEVAR A BLOQUEIO POR PARTE DO SERVIDOR
# EM PRINCÍPIO O MÓDULO ENTREZ ESTÁ PREPARADO PARA LIDAR COM ISSO MAS CONVÉM TER SEMPRE CUIDADO

def pesquisa_ncbi(email, term, db = 'pubmed', retmax = 10, rettype = 'abstract', retmode = 'text', display = 'Y' , save = 'N', ext='txt'):
    '''
    NOTA - PARA USAR RETMODE XML TÊM DE SER DECLARADA UMA NOVA VARIÁVEL POIS A FUNÇÃO AINDA NÃO CRIA FICHEIROS XML
    
    Função para a pesquisa de bibliografia na base de dados NCBI utilizando BioPython
    Esta função contempla um conjunto de parâmetros para seleção de base de dados e outputs desejado adaptando-se
    ao objetivo da pesquisa. 

    Por defeito a função pesquisa a base de dados pubmed e capta os 10 primeiros resultados para uma pesquisa
    generalizada. 

    Tem a opção de gerar um ficheiro de texto com os resultados da pesquisa de artigos e respetivos abstracts.

    No entanto a função pode ser especializada para uma pesquisa em base de dados de genes para captar ficheiros
    GenBank e FASTA para posterior análise com as funções específicas que estamos a desenvolver.

    Parâmetros
    ----------

    email : str
        email a ser utilizado pelo módulo Entrez para  pesquisa (e.g. exempl@mail.com)

    term : str
        termo a ser utilizado na pesquisa (e.g. 'gene-123 functions')

    db  : str
        base de dados onde irá ser realizada a pesquisa. Por defeito é utilizada a 'pubmed'
        Para efeitos deste trabalho iremos também utilizar: 'gene', 'nuccore' e 'protein'
        Outras opções disponíveis podem ser encontradas aqui: https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

    retmax : int
        define o número de resultados que queremos retirar para o output final (10 por defeito)

    rettype : str
        define o tipo de resultado que pretendemos retirar da base de dados ('abstract' por defeito)
        pode ser mudado para 'gb' para base de dados 'gene' ou 'fasta' para base de dados de nucleotidos ou proteina (sequencias)
        
    retmode : str
        define o retmode para a captura de resultados da pesquisa na base de dados. Por defeito 
        o retmode é text mas pode ser adaptado para 'xml' quando é pretendido obter um ficheiro xml (e.g base de dados nucleotide, ou genbank)

    display : str
        define se é pretendido o output da função legível (útil para quando se quer refinar a pesquisa antes/sem gravar o ficheiro)
        por defeito apresenta o resultado em formato legível
    
    save : str
        define se é pretendido gravar ou não o ficheiro de resultados na pasta 'lit-search-output'
        por defeito este parâmetro está inactivo

            
    
    '''
    # TODO criar ficheiros XML
    


    # Importação de módulos
    from Bio import Entrez
    import os

    # Definimos o email para os fins de pesquisa
    Entrez.email = email
        
    # Secção egquery que pode ser skipable

    handle = Entrez.egquery(term=term)
    egq_res = Entrez.read(handle)

    for _ in egq_res['eGQueryResult']:
        if db in _.values():
            if retmax > 1:
                print('Encontrados {} resultados em {}. Irão ser processados {} resultados.\n'.format(_['Count'],db,retmax))
            else:
                print('Encontrados {} resultados em {}. Irá ser processado 1 resultado.\n'.format(_['Count'],db))

    # Sacar as IDs
    handle      = Entrez.esearch(db=db,term=term, retmode=retmode, retmax=retmax)
    esearch_res = Entrez.read(handle)

    lista_ids = esearch_res['IdList']
    # print('Top 10 artigos ->',lista_ids) # Converter para uma lista de títulos mais à frente'''
    
    
    # Transfere a informação
    handle = Entrez.efetch(db=db,id=','.join(lista_ids),rettype=rettype,retmode=retmode)
    fetch_res = handle.read()
    
    # Usar com print para já, para não estar a bagunçar demasiado.
    # Perceber como se pode implementar a escrita de ficheiros para ser mais fácil de mastigar os resultados.

    # NÃO FUNCIONA NO GITHUB / LIVECODING
    # ADICIONAR PARAMETRO save = 'Y' SE SE PRETENDER GUARDAR FICHEIRO 
    
    if display.upper() == 'Y':
        print(fetch_res,sep='\n')

    # Verifica a condição de gravação de output como ficheiro de texto
    
    if retmode == 'xml': save = 'n'
    
    if save.upper() == 'Y':  

        # Definimos a extensão em função da base de dados e retmode
        if rettype == 'gb':
            ext = '.gbk'
        
        elif rettype == 'fasta':
            ext = '.fasta'
        
        elif retmode == 'xml':
            ext = '.xml'

        else: ext = '.txt'

        # Verifica se existe pasta para output e cria-a caso não exista
        cwd = os.getcwd()
        outputdir = os.path.join(cwd,'gene_search() output')
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)


        # Muda o working directory para a pasta de output para gerar o ficheiro
        os.chdir(outputdir)

        # Gera o ficheiro           
        with open(term+ext , 'w', encoding='utf-8') as _:
            _.write(fetch_res)

        # Retorna ao working directory anterior
        os.chdir(cwd)


    handle.close()

    return fetch_res

### Ferramentas de análise de features

def parsing(nome_ficheiro):

    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")

    print("ID do gene:", record.id)

    print("Nome do gene:", record.name)

    print("Descrição do gene", record.description)

    print("Comprimento da sequência:", len(record.seq), "bp")

    return

def anot(nome_ficheiro):
    
    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")

    print("Quantidade de anotações:", len(record.annotations))

    print()
    
    print("Lista de anotações:")
    
    for anotacao in record.annotations:
        print(anotacao, "->", record.annotations[anotacao])

    return   


def features_qualifiers(nome_ficheiro):
    
    from Bio import SeqIO
    
    record = SeqIO.read(nome_ficheiro, "genbank")
    
    print("Quantidade de features:", len(record.features))

    for feature in record.features:
        print(feature)
    return


