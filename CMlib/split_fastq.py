import os
from PyQt5.QtWidgets import QHeaderView, QPushButton,QProgressDialog,QProgressBar,QDialog
from PyQt5.QtCore import Qt,QBasicTimer

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])



def split_fastq(indexnow,df,fastq,output):
    '''

    :param df:  barcode csv
    :param fastq:  split fastqfile
    :param output: output directory
    :return:
    '''
    all_fastq = set()
    lev_fastq = set()
    barcodes = dict()
    barcodes_reverse = dict()


    sample = df.loc[indexnow]['Sample']
    left_code = df.loc[indexnow]['Barcode_L'].upper()
    right_code_raw = df.loc[indexnow]['Barcode_R'].upper()
    right_code = reverse_complement(right_code_raw)
    left_length = len(left_code)
    right_length = len(right_code)
    left_code_reverse = reverse_complement(left_code)
    right_code_reverse = reverse_complement(left_code)
    dirname = os.path.join(output,sample)
    mkdircmd = ' '.join(['mkdir',dirname])
    print(mkdircmd)
    os.system(mkdircmd)
    fileq = open(os.path.join(output,sample,sample + '.extendedFrags.fastq'),'w')

    flag = 0
    with open(fastq) as seqfile:
        for i in seqfile:
            inf = i.rstrip()
            flag += 1
            if flag == 1:
                name = inf
            if flag == 2:
                seq = inf
            if flag == 4:
                quality = inf
                value = name + '\n' + seq + '\n' + '+' + '\n' + quality
                flag = 0
                if seq[:left_length] == left_code and seq[-right_length:] == right_code:
                    seq1 = seq[left_length:-right_length]
                    quality1 = quality[left_length:-right_length]
                    print(name,seq1,'+',quality1,sep='\n',file=fileq)
                if seq[:right_length] == right_code_raw and seq[-left_length:] == left_code_reverse:
                    seq1 = seq[right_length:-left_length]
                    quality1 = quality[right_length:-left_length]
                    print(name,seq1,'+',quality1,sep='\n',file=fileq)
    fileq.close()
