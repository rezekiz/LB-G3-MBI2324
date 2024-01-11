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


