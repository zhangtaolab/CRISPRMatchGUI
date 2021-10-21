import os
import pysam
from pyfasta import Fasta
import matplotlib
from scipy import stats
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from os import path
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def deletion_len(samfilename, genename, output):
    outfile_del = os.path.join(output, samfilename + '_del_aln.fa')
    mut_file = os.path.join(output,'mut_rate.all.txt')
    fw = open(outfile_del,'r')
    lenth = dict()
    for i in range(1,21):
        lenth[i] = 0
    lenth['>21'] = 0
    for line in fw:
        a = line.strip()
        if a.startswith('>'):
            count_tmp = a.split('>')[1]
        else:
            delete = re.findall(r"--*",a)
            if delete:
                for i in delete:
                    size = len(i)
                    if size > 20:
                        lenth['>21'] += int(count_tmp)
                    else:
                        lenth[len(i)] +=int(count_tmp)
            else:
                if re.findall(r"\w-\w|^-\w|\w-$",a):
                    lenth['1'] += int(count_tmp)
    data = pd.read_table(mut_file,index_col=0)
    deletelentpd=pd.DataFrame.from_dict(lenth, orient='index')
    deletion_sum=deletelentpd.iloc[:, 0].sum() ##统计所有deletion
    # deletion_sum = data.loc[samfilename].total_read_count
    deletelentpd['ratio'] = deletelentpd[0]/deletion_sum

    deletelentpd_p = pd.DataFrame.from_dict(lenth, orient='index')
    deletion_sum_p=deletelentpd.iloc[:, 0].sum() ##统计所有deletion
    # deletion_sum_p = data.loc[samfilename].total_read_count  ##统计所有deletion
    deletelentpd_p['ratio'] = deletelentpd_p[0]/deletion_sum_p*100
    return(deletelentpd,deletelentpd_p)

def deletionbarplot(regmean,stdrr,namenow,deletelentpdall,pdfname):
    fig, ax = plt.subplots()
    x = np.array(range(1, 22))
    ymajorFormatter = FormatStrFormatter('%1.1f') ## 设置坐标轴格式
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.bar(x, regmean, color='red')
    # add errorbar, elinewidth：errorbar line with; capsize/capthick:上下横线长短／粗细，ls:linestyle='None'去掉连接线。 ecolor: errorbar line color
    ax.errorbar(x, regmean, yerr=stdrr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None',ecolor='black')
    ax.set_title(namenow, size=15, fontdict={'family': 'Times New Roman'})
    ax.set_ylabel('Deletion Size (%)', size=15, fontdict={'family': 'Times New Roman'})
    ax.set_xticks(x)
    ax.set_xticklabels(deletelentpdall.index, rotation=35, fontdict={'family': 'Arial'}, size=5)
    #ax.legend((bar1[0], bar2[0], bar3[0]), ('rep1', 'rep2', 'control'))
    plt.savefig(pdfname)
    plt.close(fig)
    print("group", namenow, "deletion size finished!")

def barchart_filter(groupinfo, output, bamdir):

    groupinfor = pd.read_csv(groupinfo)
    #print(groupinfor)
    #groupinfor = groupinfor.dropna(axis=0, how='any',thresh=7)
    groupinfor = groupinfor.fillna("UNKNOWN")  ##填充表格中NaN处
    #print(groupinfor)

    for idx in groupinfor.index:

        repbam1 = os.path.join(bamdir, groupinfor.loc[idx]['rep1'] + '.bam')
        repbam2 = os.path.join(bamdir, groupinfor.loc[idx]['rep2'] + '.bam')
        repbam3 = os.path.join(bamdir, groupinfor.loc[idx]['rep3'] + '.bam')

        repdel1 = groupinfor.loc[idx]['rep1']
        repdel2 = groupinfor.loc[idx]['rep2']
        repdel3 = groupinfor.loc[idx]['rep3']
        #ckbam = os.path.join(bamdir, groupinfor.loc[idx]['control'] + '.bam')
        genename = groupinfor.loc[idx]['gene']
        namenow = groupinfor.loc[idx]['group']
        pdfname = os.path.join(output, namenow + '_deletion_size.pdf')
        csvname = os.path.join(output, namenow + '_deletion_size.csv')

        #if (path.exists(repbam1) and path.exists(repbam2)) and path.exists(ckbam):
        if (path.exists(repbam1) and path.exists(repbam2) and path.exists(repbam3)):
            (deletelentpd1,deletelentpd1_p) = deletion_len(samfilename=repdel1, genename=genename, output=output)
            (deletelentpd2,deletelentpd2_p) = deletion_len(samfilename=repdel2, genename=genename, output=output)
            (deletelentpd3, deletelentpd3_p) = deletion_len(samfilename=repdel3, genename=genename, output=output)
            deletelentpdall = deletelentpd1
            #deletelentpdCK = deletion_len(samfilename=ckbam, genename=genename)
            #deletelentpd1 = pd.DataFrame.from_dict(pd1, orient='index')
            col1 = deletelentpd1_p.iloc[:, 1] ##提取pandas表格的第三列数值
            y1 = col1.values
            #deletelentpd2 = pd.DataFrame.from_dict(pd2, orient='index')
            col2 = deletelentpd2_p.iloc[:, 1]
            y2 = col2.values
            col3 = deletelentpd3_p.iloc[:, 1]
            y3 = col3.values
            #deletelentpdCK = pd.DataFrame.from_dict(pdCK, orient='index')
            #colCK = deletelentpdCK.iloc[:, 1]
            #yCK = colCK.values
            reg = pd.concat([col1, col2, col3], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd1, deletelentpd2, deletelentpd3], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep1_count','rep1_Ratio','rep2_count','rep2_Ratio','rep3_count','rep3_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif (path.exists(repbam1) and path.exists(repbam2)):
            print("Rep 3 is missing")
            (deletelentpd1,deletelentpd1_p) = deletion_len(samfilename=repdel1, genename=genename, output=output)
            (deletelentpd2,deletelentpd2_p) = deletion_len(samfilename=repdel2, genename=genename, output=output)
            deletelentpdall = deletelentpd1
            #deletelentpdCK = deletion_len(samfilename=ckbam, genename=genename)
            #deletelentpd1 = pd.DataFrame.from_dict(pd1, orient='index')
            col1 = deletelentpd1_p.iloc[:, 1] ##提取pandas表格的第三列数值
            y1 = col1.values
            #deletelentpd2 = pd.DataFrame.from_dict(pd2, orient='index')
            col2 = deletelentpd2_p.iloc[:, 1]
            y2 = col2.values
            #deletelentpdCK = pd.DataFrame.from_dict(pdCK, orient='index')
            #colCK = deletelentpdCK.iloc[:, 1]
            #yCK = colCK.values
            reg = pd.concat([col1, col2], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd1, deletelentpd2], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep1_count','rep1_Ratio','rep2_count','rep2_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif (path.exists(repbam1) and path.exists(repbam3)):
            print("Rep 2 is missing")
            (deletelentpd1,deletelentpd1_p) = deletion_len(samfilename=repdel1, genename=genename, output=output)
            (deletelentpd3,deletelentpd3_p) = deletion_len(samfilename=repdel3, genename=genename, output=output)
            deletelentpdall = deletelentpd1
            #deletelentpdCK = deletion_len(samfilename=ckbam, genename=genename)
            #deletelentpd1 = pd.DataFrame.from_dict(pd1, orient='index')
            col1 = deletelentpd1_p.iloc[:, 1] ##提取pandas表格的第三列数值
            y1 = col1.values
            #deletelentpd2 = pd.DataFrame.from_dict(pd2, orient='index')
            col3 = deletelentpd3_p.iloc[:, 1]
            y2 = col3.values
            #deletelentpdCK = pd.DataFrame.from_dict(pdCK, orient='index')
            #colCK = deletelentpdCK.iloc[:, 1]
            #yCK = colCK.values
            reg = pd.concat([col1, col3], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd1, deletelentpd3], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep1_count','rep1_Ratio','rep3_count','rep3_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif (path.exists(repbam2) and path.exists(repbam3)):
            print("Rep 1 is missing")
            (deletelentpd3,deletelentpd3_p) = deletion_len(samfilename=repdel3, genename=genename, output=output)
            (deletelentpd2,deletelentpd2_p) = deletion_len(samfilename=repdel2, genename=genename, output=output)
            deletelentpdall = deletelentpd3
            #deletelentpdCK = deletion_len(samfilename=ckbam, genename=genename)
            #deletelentpd1 = pd.DataFrame.from_dict(pd1, orient='index')
            col3 = deletelentpd3_p.iloc[:, 1] ##提取pandas表格的第三列数值
            y3 = col3.values
            #deletelentpd2 = pd.DataFrame.from_dict(pd2, orient='index')
            col2 = deletelentpd2_p.iloc[:, 1]
            y2 = col2.values
            #deletelentpdCK = pd.DataFrame.from_dict(pdCK, orient='index')
            #colCK = deletelentpdCK.iloc[:, 1]
            #yCK = colCK.values
            reg = pd.concat([col3, col2], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd3, deletelentpd2], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep2_count','rep2_Ratio','rep3_count','rep3_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif path.exists(repbam1):
            print("Rep2 and Rep3 are missing")
            (deletelentpd1,deletelentpd1_p) = deletion_len(samfilename=repdel1, genename=genename, output=output)
            deletelentpdall = deletelentpd1
            col1 = deletelentpd1_p.iloc[:, 1] ##提取pandas表格的第三列数值

            reg = pd.concat([col1], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd1], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep1_count','rep1_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif path.exists(repbam2):
            print("Rep1 and Rep3 are missing")
            (deletelentpd2,deletelentpd2_p) = deletion_len(samfilename=repdel2, genename=genename, output=output)
            deletelentpdall = deletelentpd2
            col2 = deletelentpd2_p.iloc[:, 1] ##提取pandas表格的第三列数值

            reg = pd.concat([col2], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd2], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep2_count','rep2_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)

        elif path.exists(repbam3):
            print("Rep1 and Rep2 are missing")
            (deletelentpd3,deletelentpd3_p) = deletion_len(samfilename=repdel3, genename=genename, output=output)
            deletelentpdall = deletelentpd3
            col3 = deletelentpd3_p.iloc[:, 1] ##提取pandas表格的第三列数值
            reg = pd.concat([col3], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)
            len_raw_data = pd.concat([deletelentpd3], axis=1)
            len_raw_data['ratio_mean'] = regmean
            len_raw_data['ratio_stdrr'] = stdrr
            len_raw_data.to_csv(csvname, index=True, index_label='Index', header=['rep3_count','rep3_Ratio','ratio_mean','ratio_stdrr'], sep=',',encoding='utf-8')
            print("Output",csvname, "finished!" )

            ##plot figures
            deletionbarplot(regmean, stdrr, namenow, deletelentpdall, pdfname)


