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

    
    if rettype == 'fasta' or rettype == 'gb': retmax = 1

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

        filename = term+ext

        # Verifica se existe pasta para output e cria-a caso não exista
        cwd = os.getcwd()
        outputdir = os.path.join(cwd,'output')
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)


        # Muda o working directory para a pasta de output para gerar o ficheiro
        os.chdir(outputdir)

        if os.path.isfile(filename): print('Ficheiro já existe. Não será criado novo ficheiro.')



        # Gera o ficheiro           
        with open(term+ext , 'w', encoding='utf-8') as _:
            _.write(fetch_res)

        # Retorna ao working directory anterior
        os.chdir(cwd)

        print('Ficheiro gravado.')


    handle.close()

    print('Concluído.')

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


##############################
# Secção de Análise Proteína #
##############################

# Carrega a sequência de AA
def carregar_sequencia(filename):
    from Bio import SeqIO

    try:    
        # Lemos o ficheiro .fasta
        seq_proteina = SeqIO.read(open(filename),format='fasta')

        sp = seq_proteina.seq # Isolamos a sequência
        return sp
    
    except FileNotFoundError:
        print('Ficheiro FASTA não existe. Apresentar um ficheiro FASTA válido.')

# Contagem de AA

def conta_aa(fastafile):
    # TODO colocar um dicionário bonitinho
    from collections import Counter
    seq_prot = carregar_sequencia(fastafile)
    contagem_aa = dict(Counter(seq_prot))
    
    return contagem_aa


# Parsing de Blast Records 

def blast_parse(filename, lim = None, ethreshold = 0.05):
    
    from Bio.Blast import NCBIXML

    handle = open(filename)
    blast_res = NCBIXML.parse(handle)
    blast_records = list(blast_res)

    num = 1
    for record in blast_records:

        if lim == None: lim = len(record.alignments)    

        for aligment in record.alignments[0:lim]:
            #print(aligment)
            
            for hsp in aligment.hsps:
                #print(hsp)
                
                if hsp.expect < ethreshold:
                    print('****Alinhamento {} ****'.format(num))
                    print('Sequência:',aligment.title)
                    print('Tamanho da Sequência:',aligment.length)
                    print('e-value:',hsp.expect)
                    print(hsp.query[0:40]+'...')
                    print(hsp.match[0:40]+'...')
                    print(hsp.sbjct[0:40]+'...')
                    print()
                    num += 1    

# Analise na BD SwissProt
# TODO integrar esta função na função de análise

def swiss_prot_scan(swiss_id):
    from Bio import SeqIO
    from Bio import ExPASy
    # Efetuamos uma pesquisa na base de dados SwissProt 
    handle = ExPASy.get_sprot_raw(swiss_id) # Encontramos o ID ao ler o ficheiro fasta
    sr = SeqIO.read(handle, "swiss")
    print(
        f'ID {sr.id}',
        f'Sequência: {sr.seq}',
        f'Tamanho da sequência: {len(sr.seq)} bp',
        f'Nome: {sr.name}',
        f'Descrição: {sr.description}',
        f'Taxonomia: {sr.annotations["taxonomy"]}',
        f'Organismo: {sr.annotations["organism"]}',
        f'Keywords: {sr.annotations["keywords"]}',
        sep = '\n')
    

# Função que executa análise Blast caso o ficheiro blast não exista para fazer parsing
    
def analise_blast(email, fastafile, xmlfile, lim = None, ethreshold = 0.05):
    '''
    Função que carrega uma sequência de aminoácidos a partir de ficheiro fasta e realiza
    uma análise blast caso esta ainda não esteja feita ou então faz parse de blast records
    até ao limite definido.

    email : str
        parâmetro de email para utilização do serviço NCBIWWW

    fastafile : str
        caminho para ficheiro fasta com a sequência da proteína

    xmlfile : str
        caminho para o ficheiro XML com blastrecords caso já esteja previamente realizada

    lim : int
        limite de resultados que queremos visualizar dos resultados do blast
        por defeito não está definido

    e-threshold : float
        valor-limite de e-value a que queremos delimitar a nossa análise blast
        por defeito está definido 0.05
    
    '''
    import os
    from Bio.Blast import NCBIWWW
    
    NCBIWWW.email = email

    if not fastafile.endswith('.fasta'): raise ValueError('Ficheiro .fasta inválido')

    # Muda para a pasta da função para as operações
        
    cwd = os.getcwd()
    outputdir = os.path.join(cwd,'output')
    os.chdir(outputdir)
    
                

    try:
        sequencia_proteina = carregar_sequencia(fastafile)
    except FileNotFoundError: 
        print('Ficheiro FASTA não existe. Apresentar um ficheiro FASTA válido.')
        print('Sugestão: Realizar um pesquisa_ncbi() com o ID da proteína na base de dados protein')
    
    # Verifica se já existe um ficheiro XML para os resultados do BLAST. Ignora a pesquisa caso já exista.
    # Para poupar tempo.
    # Se o ficheiro já existir executa o parsing do ficheiro.

    if not xmlfile.endswith('.xml'):
        xmlfile = xmlfile + '.xml'

    # Verifica a existência do ficheiro

    if os.path.isfile(xmlfile): 
       
        
        blast_parse(xmlfile,lim, ethreshold)
        os.chdir(cwd)
        
    
    else: 
    # Executa o BLAST da Proteína
        handle = NCBIWWW.qblast('blastp','swissprot',sequencia_proteina)

        # TODO Converter o savefile numa sub-função
        # Verifica se existe pasta para output e cria-a caso não exista
        cwd = os.getcwd()
        outputdir = os.path.join(cwd,'output')
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)

        # Muda o working directory para a pasta de output para gerar o ficheiro
        os.chdir(outputdir)

        with open(xmlfile, 'w') as _:
            _.write(handle.read())

        blast_parse(xmlfile,lim,ethreshold)

        # Retorna ao working directory anterior
        os.chdir(cwd)
        

        

        print(f'Ficheiro {xmlfile} gravado. Para realizar um blast parse ir buscar o ficheiro à pasta output')
        
        handle.close()