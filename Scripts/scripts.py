####################################
# Secção de Pesquisa Bibliográfica #
####################################

def pesquisa_ncbi(email, term, db = 'pubmed', retmax = 10, rettype = 'abstract', retmode = 'text', display = 'Y' , save = 'N', ext='txt'):
    '''

    Desenvolvido por Rui Sousa

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

    if rettype == 'fasta' or rettype == 'gb': retmax = 1

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
    
    # Transfere a informação
    handle = Entrez.efetch(db=db,id=','.join(lista_ids),rettype=rettype,retmode=retmode)
    fetch_res = handle.read()
    
    if display.upper() == 'Y':
        print(fetch_res,sep='\n')


    # Verifica a condição de gravação de output como ficheiro de texto
    
    if retmode == 'xml': save = 'n'
    
    if save.upper() == 'Y':  

        # Definimos a extensão em função da base de dados e retmode
        

        filename = term + type_handler(rettype,retmode)

        if os.path.isfile(filename): 
            print('Ficheiro já existe. Não será criado novo ficheiro.')
            return fetch_res

        # Gera o ficheiro           
        guardar_ficheiro(filename, fetch_res)


    handle.close()

    return fetch_res


# SUBFUNÇÕES

def guardar_ficheiro(ficheiro, content):

    '''
    Função que recebe um nome de ficheiro (ficheiro) e carrega o conteúdo 
    da variável content para esse ficheiro.
    '''    
    with open(ficheiro, 'w', encoding = 'utf-8') as _:
        _.write(content)

    print(f'{ficheiro} guardado com sucesso.')

def type_handler(rettype,retmode):

    '''
    Função auxiliar da função pesquisa NCBI que permite
    controlar o tipo de ficheiro de output da função
    '''

    if rettype == 'gb':
        ext = '.gbk'
        
    elif rettype == 'fasta':
        ext = '.fasta'
        
    elif retmode == 'xml':
        ext = '.xml'

    else: ext = '.txt'

    return ext

# PESQUISA KEGG

def get_kegg(organismo,idgene,displ = 1):
    '''
    Desenvolvido por Rui Sousa

    Função que pesquisa na base de dados KEGG.

    Na primeira versão vai ser preciso introduzir o ID do gene (consultar a base de dados gene, a função 
    pesquisa_ncbi pode ajudar aqui)

    Referencial de IDs, etc: https://www.kegg.jp/kegg/rest/keggapi.html
    
    Numa segunda versão será implementada uma pesquisa com prompt com termo de pesquisa e confirmação se 
    é o gene pretendido para depois transferir o ficheiro relativo ao gene

    Parametros
    ----------

    organismo : str
        identificador do organismo de acordo com o referencial de IDs da API KEGG

    idgene : str
        identificador do gene de acordo com a base de dados GENE da NCBI

    displ : int
        trigger para visualização de informação (ligado por defeito)

    Devolve
    -------

    Ficheiro de texto captado na página do KEGG referente ao gene em questão.

    '''
    import os
    from Bio.KEGG import REST

    if displ != 1 or 0:
        raise ValueError('Display tem de ser 1 ou 0.')

    request = REST.kegg_get(organismo+':'+idgene)

    ficheiro = organismo+idgene+'.kegg'

    if not os.path.isfile(ficheiro):
        with open(ficheiro,'w',encoding='utf-8') as _:
            _.write(request.read())
    
    if displ == 1:
        read(ficheiro)
    
    return request.read()

def read(ficheiro):
    handle = open(ficheiro)
    for line in handle:
        print(line)
    handle.close()


#################################
# Secção de análise de features #
#################################

# Desenvolvido por Mariana Oliveira

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

'''
(DEPRECATED)
Funções desenvolvidas por Rui Sousa.

Conjunto de funções para processamento básico de sequências proteicas (contagem de aminoácidos, analisar dados da base de dados swiss prot).
'''
def carregar_sequencia(filename):
    from Bio import SeqIO

    try:    
        # Lemos o ficheiro .fasta
        seq_proteina = SeqIO.read(open(filename),format='fasta')

        sp = seq_proteina.seq # Isolamos a sequência
        return sp
    
    except FileNotFoundError:
        print('Ficheiro FASTA não existe. Apresentar um ficheiro FASTA válido.')

def conta_aa(fastafile):
    # TODO colocar um dicionário bonitinho
    from collections import Counter
    seq_prot = carregar_sequencia(fastafile)
    contagem_aa = dict(Counter(seq_prot))
    
    return contagem_aa

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


    
#############################################
# Secção de Análise de Homologias por BLAST #
#############################################

'''
(DEPRECATED)
Funções desenvolvidas por Rui Sousa

Conjunto de funções para instanciação de blast remoto via package BioPython.
'''
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
        
    
    else: 
    # Executa o BLAST da Proteína
        handle = NCBIWWW.qblast('blastp','swissprot',sequencia_proteina)
        with open(xmlfile, 'w') as _:
            _.write(handle.read())

        blast_parse(xmlfile,lim,ethreshold)

        # Retorna ao working directory anterior
        

        

        print(f'Ficheiro {xmlfile} gravado. Para realizar um blast parse ir buscar o ficheiro à pasta output')
        
        handle.close()


# Funções desenvolvidas por Carlos Gomes, usadas ao longo do trabalho

def blast_remoto(arquivo_entrada, formato = 'fasta', banco_dados = 'nr', tipo_blast = 'blastp', save = 'N'):

    

    '''
    Função para a realização do Blast Remoto

    Parâmetros:
    arquivo_entrada : str
        Caminho para o arquivo de entrada contendo a sequência que será usada como consulta

    formato : str
        Formato do arquivo de entrada
    
    banco_dados : str
        Banco de dados usado para a busca BLAST ('nr' para banco de dados não redundante do NCBI)
    
    tipo_blast : str
        Tipo de busca BLAST ('blastp' para proteínas)

    save : str
        Parâmetro para salvar ou não o ficheiro recorrendo à função save_blast()
    
    
    Retorna:
    resultado : handle
        Objeto que contém os resultados do BLAST
    '''

    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    

    # Lê a sequência do arquivo
    record = SeqIO.read(open(arquivo_entrada), format = formato)

    # Executa o Blast Remoto
    resultado = NCBIWWW.qblast(tipo_blast, banco_dados, record.format(formato))

    # Guarda os resultados num ficheiro XML se pretendido

    if save.upper() == 'Y':
        save_blast(resultado)
    

    return resultado


# Guardar Blast
def save_blast(resultado, nome_arquivo = 'blast.xml'):

    '''
    Função para salvar os resultados do BLAST num arquivo XML

    Parâmetros:
    resultado : handle
        Objeto que contém os resultados do BLAST
    
    nome_arquivo : str
        Nome do arquivo onde os resultados vão ser guardados 
    
    O 'with' fecha automaticamente o arquivo 
    '''
    
    with open(nome_arquivo,'w') as save_file:
        save_file.write(resultado.read())

def resultado_blast(nome_arquivo='blast.xml'):
    '''
    Função para ler e processar os resultados do BLAST armazenados num arquivo XML

    Parâmetros:
    nome_arquivo : str
        Nome do arquivo XML que contêm os resultados do BLAST 

    Retorna:
    blast_record
        Objeto que representa os resultados do BLAST
    
    Raises:
    FileNotFoundError:
        Se o arquivo especificado pelo 'nome_arquivo' não for encontrado

    IOError:
        Se ocorrer um erro durante a leitura do arquivo.
    '''
    from Bio.Blast import NCBIXML
    
    try:
        # Abre o arquivo XML com os resultados do BLAST
        with open(nome_arquivo) as resultado:
            # Lê o conteúdo do arquivo e retorna o Blast record
            blast_record = NCBIXML.read(resultado)
    except FileNotFoundError as e:
        # Captura e relança a exceção para fornecer uma mensagem de erro mais informativa
        raise FileNotFoundError(f'O arquivo {nome_arquivo} não foi encontrado.') from e
    except IOError as e:
        # Captura e relança a exceção para fornecer uma mensagem de erro mais informativa
        raise IOError(f'Erro ao ler ou processar o arquivo {nome_arquivo}.') from e
    
    return blast_record

def parametros(blast_record):
    
    '''
    Imprime os parâmetros do BLAST

    Parâmetros:
    blast_record 
        Objeto que representa os resultados do BLAST

    Prints:
    Database - str
        Nome do banco de dados utilizado no BLAST
    
    Matrix - str
        Nome da matriz de substituição utilizada
    
    Gap penalties - tuple
        Tuple contendo os gap penalties
    
    '''
    print('PARAMETROS:')
    print('Database - ', blast_record.database)
    print('Matrix - ', blast_record.matrix)
    print('Gap penalties - ', blast_record.gap_penalties)


def numero_de_hits(blast_record):
    
    '''
    Imprime o número de hits no resultado do BLAST

    Parâmetros:
    blast_record 
        Objeto que representa os resultados do Blast

    Prints:
    Número de Hits - int
        Número total de alinhamentos encontrados durante a busca BLAST.
    '''

    print('Número de Hits: ', len(blast_record.alignments))



def primeiro_alinhamento(blast_record):
    
    '''
    Imprime informações sobre o primeiro alinhamento do resultado do BLAST.

    Parâmetros:
    blast_record 
        Objeto que representa os resultados do BLAST

    Prints:
    Primeiro Alinhamento:
    
    Acession - str
        Número de acesso associado
    
    Hit Id - str
        Identificador único associado 
    
    Definição - str
        Definição 
    
    Length - int
        Comprimento 
    
    HSPs - int
        Número de High Scoring Pairs (HSPs)
    
    E-value - float
        Valor de E associado ao primeiro HSP
    
    Score - float
        Pontuação do primeiro HSP
    
    Length - int
        Comprimento do alinhamento do primeiro HSP
    '''
    
    
    # Obtém o primeiro alinhamento do BLAST
    first_alignment = blast_record.alignments[0]

    # Imprime informações sobre o primeiro alinhamento
    print('Primeiro Alinhamento:')
    print('Acession: ', first_alignment.accession)
    print('Hit Id: ', first_alignment.hit_id)
    print('Definição: ', first_alignment.hit_def)
    print('Length: ', first_alignment.length)
    print('HSPs: ', len(first_alignment.hsps))

    # Obtém o primeiro HSP (High Scoring Pair) do primeiro alinhamento
    first_hsp = first_alignment.hsps[0]

    # Imprime informações sobre o primeiro HSP
    print('E-value: ', first_hsp.expect)
    print('Score: ', first_hsp.score)
    print('Length: ', first_hsp.align_length)


    
def alinhamentos(blast_record, num_alinhamentos=5):
    '''
    Imprime informações sobre os alinhamentos no resultado da busca BLAST.

    Parâmetros:
    blast_record 
        Objeto que representa os resultados do BLAST

    num_alinhamentos : int
        Número máximo de alinhamentos

    Prints:
    Informações sobre os alinhamentos:
            Sequence: str
                Título
            Accession: str
                Número de acesso
            Definition: str
                Definição
            E-Value: float
                Valor de E associado 
    '''
    for i in range(min(num_alinhamentos, len(blast_record.alignments))):
        alignment = blast_record.alignments[i]
        print(f'Alinhamento {i + 1}:')
        print('Sequence: ', alignment.title)
        print('Accession: ', alignment.accession)
        print('Definition: ', alignment.hit_def)
        for hsp in alignment.hsps:
            print('E-Value: ', hsp.expect)
        print('\n')

##############################
#  Análise de estruturas 3D  #
##############################

def PDB_viewer(id,ficheiro):

    '''
    Desenvolvido por Rui Sousa

    Função que instancia o módulo PDBParser do package BioPython para a visualização 
    de ficheiros PDB

    Roadmap: Implementar uma função para pesquisa para transferir ficheiros diretamente da base de dados PDB
    a partir do terminal usando o package rcsbsearchapi

    Parametros
    ----------

    id : str
        ID da proteína na base de dados PDB

    ficheiro : str
        Caminho para o ficheiro PDB

    Funcionalidade
    --------------

    Gera um modelo 3D interativo da estrutura da proteína
    '''


    from Bio.PDB import PDBParser
    import nglview as nv

    # Carrega o ficheiro e constroi a estrutura 3D para visualizar
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(id,ficheiro)

    # Visualiza a estrutura criada usando o nglviewer

    return nv.show_biopython(structure)


#######################
#  Análise de motifs  #
#######################

def motif_viewer(ficheiro, disp_logo = 1):
    '''
    Desenvolvido por Rui Sousa

    Função que recebe como input um ficheiro de resultados MEME em XML e devolve a lista de motifs
    e os logos dos diferentes motifs.

    Roadmap: Adicionar uma interação com a instanciação do MEME usando Docker
    '''

    from Bio.motifs import meme
    from logomaker import Logo
    import pandas as pd

    with open(ficheiro) as _:
        record = meme.read(_)

    for motif in record:
        for instance in motif.instances:
            print(instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue)
    
    if disp_logo == 1:

        for motif in record:
            mdf = pd.DataFrame(motif.counts)
            logo = Logo(mdf)

