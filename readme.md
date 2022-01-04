# CRISPRMatchGUI
## Brief introduction
The Graphical User Interface(GUI) for CRISPRMatch--An automatic calculation and visualization tool for high-throughput CRISPR genome-editing data analysis
## I. <u>Requirements(软件所需依赖包)</u>
Anaconda</br>
python3</br>
bwa</br>
samtools</br>
FLASH</br>
pyqt5</br>

[**<font color="#FF0000">Note:</font>**] Using `Anaconda` to Install all packages (`bwa,samtools,picard,FLASH`) ##应用conda统一安装即可，方便快捷

## II. <u>Manually Install（手动安装，非虚拟机）</u>
CentOS Linux release 7.3.1611 (terminal)
1. Install Anaconda</br>
```
$ yum install wget git  ##安装git和wget程序
$ mkdir /home/software  ##创建下载软件的文件夹，这里以software为例
$ cd /home/software     ##进入到software文件夹下
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh  ##下载linux对应的conda版本，注意尽量选择3.5版本的！
$ bash Anaconda3-5.0.1-Linux-x86_64.sh  ##用bash命令安装conda
```
2. Install required packages  ##利用conda安装所有依赖包（建议用清华软件源镜像替换，缩短安装时间）
```
$ conda install bwa \             ##用\符号和回车分隔多个软件
                samtools \  
                pyqt=5.6 \  
                flash \ 
                matplotlib \  
                pysam \  
                pandas \  
                argparse \  
                numpy \
```
**<font color="#FF0000">Note:</font>** To ensure the tool working, please using `Anaconda` to install all packages (`bwa,samtools,pyqt,FLASH ...`)

3. Download CRISPRMatchGUI and test  ##下载本软件的软件包
```
$ cd /home/software   ##进入到software文件夹下
$ git clone https://github.com/zhangtaolab/CRISPRMatchGUI.git  ##利用git方式下载本软件包
$ cd /home/software/CRISPRMatchGUI/   ##进入本软件文件夹
$ python3 /home/software/CRISPRMatchGUI/start.py  ##使用python3打开软件包中的start.py程序，即可实现软件运行
  
```
## III. <u>Start running（运行方法）</u>
1. Video manual(用户手册)</br>

>(1)CRISPRMatch虚拟机使用教程
- Link: https://v.youku.com/v_show/id_XMzgwODc4ODQ2NA==.html?spm=a2h3j.8428770.3416059.1

>(2)双端测序数据合并教程
- Link: https://v.youku.com/v_show/id_XMzkzMTY5NTEwOA==.html?scm=20140719.manual.114461.video_XMzkzMTY5NTEwOA==

>(3)拆分混池测序结果(带有barcode信息)
- Link: https://v.youku.com/v_show/id_XMzkzMTY5MzY4NA==.html?scm=20140719.manual.114461.video_XMzkzMTY5MzY4NA==

>(4)虚拟机读取usb设备(改方法可实现大数据集计算)
- Link: https://v.youku.com/v_show/id_XMzk0MDgyMjA2MA==.html?scm=20140719.manual.114461.video_XMzk0MDgyMjA2MA==

2. Mirroring file for Windows (虚拟机下载地址)</br>
- Link: https://pan.baidu.com/s/1L8KPij9SP2Mp9v7RYgS5_w  code: CPF1

3. Files for mutation calculation(编辑计算所需三个信息文件)</br>
- **File1**: Genome-editing target sequences  
[Fasta format example](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/Samples_gene.fa)
- **File2**: NGS samples information
<font color="#FF0000">*note*:</font>
For CRISPR-Cas9 system, the `'Note'` must contain `'gRNA'` label.  
For CRISPR-Cpf1 system, the `'Note'` must contain `'crRNA'` label.  
*example*:  
[sample information](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/sample_infor.csv)  
- **File3**: NGS group information  
<font color="#FF0000">*note*:</font> At present, two repeats are supported<br>
*example*:</br>
[group information](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/sample_test/group_info.csv)  
- **Note**: the information files `File1`, `File2` and `File3` are required!  
</br>

4. Merge paired-end reads(运行双端测序数据合并程序)</br>
```
$ cd /home/software/CRISPRMatchGUI/               ##进入本软件文件夹
$ python3 /home/software/CRISPRMatchGUI/merge.py  ##运行双端测序数据合并
```
- *example*:<br/>
[paired-end reads](https://github.com/zhangtaolab/CRISPRMatchGUI/tree/master/merge_sample/)  
</br>

5. Split sequencing file(运行拆分混池测序结果)</br>
```
$ cd /home/software/CRISPRMatchGUI/               ##进入本软件文件夹 
$ python3 /home/software/CRISPRMatchGUI/split.py  ##运行拆分混池测序程序
```

