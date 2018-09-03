# CRISPRMatchGUI
## Brief introduction
An automatic calculation and visualization tool for high-throughput CRISPR genome-editing data analysis
## I. <u>Requirements</u>
Anaconda</br>
python3</br>
bwa</br>
samtools</br>
FLASH</br>
pyqt5</br>

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **<font color=red>Note:</font>** Using `Anaconda` to Install all packages (`bwa,samtools,picard,FLASH`)

## II. <u>Manually Install</u>
CentOS Linux release 7.3.1611 (terminal)
1. Install Anaconda</br>
```
$ yum install wget git
$ mkdir /home/software
$ cd /home/software
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
$ bash Anaconda3-5.0.1-Linux-x86_64.sh
```
2. Install required packages  
```
$ conda install bwa \  
                samtools \  
                pyqt=5.6 \  
                flash \ 
                matplotlib \  
                pysam \  
                pandas \  
                argparse \  
                numpy \
```
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **<font color=red>Note:</font>** To ensure the tool working, please using `Anaconda` to install all packages (`bwa,samtools,pyqt,FLASH ...`)

3. Download CRISPRMatchGUI and test
```
$ cd /home/software
$ git clone https://github.com/zhangtaolab/CRISPRMatchGUI.git
$ python3 /home/software/CRISPRMatchGUI/start.py
  
```
## III. <u>Start running</u>
1. Files for mutation calculation  
- **File1**: Genome-editing target sequences  
[Fasta format example](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/Samples_gene.fa)
- **File2**: NGS samples information  
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) <font color=red>*note*:</font>   
For CRISPR-Cas9 system, the `'Note'` must contain `'gRNA'` label.  
For CRISPR-Cpf1 system, the `'Note'` must contain `'crRNA'` label.  
*example*:  
[sample information](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/sample_infor.csv)  
- **File3**: NGS group information  
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) <font color=red>*note*:</font> At present, two repeats are supported<br>
*example*:</br>
[group information](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/group_info.csv)  
- **Note**: the information files `File1`, `File2` and `File3` are required!  
</br>
2. Video manual:</br>

- Link: https://v.youku.com/v_show/id_XMzgwODc4ODQ2NA==.html?spm=a2h3j.8428770.3416059.1


3. Merge paired-end reads:</br>
```
$ cd /home/software/CRISPRMatchGUI/
$ python3 /home/software/CRISPRMatchGUI/merge.py
```
- *example*:<br/>
[paired-end reads](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/merge_sample/)  
