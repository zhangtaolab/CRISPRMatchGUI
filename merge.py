import sys
import pandas as pd
import matplotlib.pyplot as plt

from PyQt5 import uic,QtWidgets
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QHeaderView
import os

from subprocess import Popen, PIPE

path = os.getcwd()
qtCreatorFile = os.path.join(path,'CMlib/flash_merge.ui')

Ui_showtable, QtBaseClass = uic.loadUiType(qtCreatorFile)

class showtable(QtWidgets.QDialog, Ui_showtable):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        Ui_showtable.__init__(self)
        self.setupUi(self)

        self.setWindowTitle('Merge FastQ')
        self.leftbtn.clicked.connect(lambda: self.getfastq("left"))
        self.left.setReadOnly(True)  ##设置不可输入
        self.rightbtn.clicked.connect(lambda: self.getfastq("right"))
        self.right.setReadOnly(True)
        self.outputbtn.clicked.connect(self.outputdir)
        self.output.setReadOnly(True)
        self.pushButton.clicked.connect(self.merge)




    def getfastq(self,file):
        fastqPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', path)
        if fastqPath != "":
            if file == "left":
                print("fastq Direction", fastqPath)
                self.left.setText(fastqPath)
                self.leftfastq = fastqPath
            if file == "right":
                print("fastq Direction", fastqPath)
                self.right.setText(fastqPath)
                self.rightfastq = fastqPath

    def outputdir(self):
        outputdirpath = QtWidgets.QFileDialog.getExistingDirectory(self, 'open directory', path)
        if outputdirpath != "":
            print("Direction", outputdirpath)
            self.outputdirpath = outputdirpath
            self.output.setText(outputdirpath)
            self.outputfiledir = outputdirpath

    def merge(self):
        outname = self.name.text().rstrip()
        threadnumber = self.spinBox.value()
        if outname != "" and self.leftfastq != "" and self.rightfastq != "" and self.outputfiledir != "":
            flashpath = self.which('flash')
            if flashpath:
                flashversion = self.flash('flash')
                if flashversion == 'None':
                    self.showMessageBox("warning","Please input flash directory")
                else:
                    flashbin=flashpath[0]
                    flashcmd = ' '.join([flashbin, '-o', outname, '-t', str(threadnumber), '-d', self.outputfiledir, self.leftfastq, self.rightfastq, '2>&1 | tee', os.path.join(self.outputfiledir, outname + '_flash.log')])
                    print(flashcmd)
                    runflash = Popen(flashcmd, shell=True)
                    runflash.communicate()


                    msgBox = QtWidgets.QMessageBox()
                    msgBox.setWindowTitle("Information")
                    msgBox.setIcon(QtWidgets.QMessageBox.Information)
                    msgBox.setText("Project Done!")
                    msgBox.setDetailedText(''.join(['File ',outname,'.extendedFrags.fastq ','is located in ',self.outputfiledir,'/']))
                    msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
                    msgBox.exec_()

        else:
            self.showMessageBox("warning","Please set output name!")

    def flash(self,filename):
        """
        :param filename:
        :return: flash version
        """
        flashpath = self.which(filename)
        flashcmd = ' '.join([flashpath[0], '--version'])
        # location= samtoolspath[0]
        flashrun = Popen(flashcmd, stdout=PIPE, stderr=PIPE, shell=True)
        i = flashrun.stdout.readlines()[0]
        version = i.decode('utf-8').rstrip('\n')
        flashrun.communicate()
        return version

    def which(self,filename):
        """docstring for which"""
        locations = os.environ.get("PATH").split(os.pathsep)
        candidates = []
        for location in locations:
            candidate = os.path.join(location, filename)
            if os.path.isfile(candidate):
                candidates.append(candidate)
        return candidates


    # ############## warning message #########
    def showMessageBox(self, title, message):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Warning)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.exec_()
    ##################################################

    #     self.checkbtn.clicked.connect(self.sampleEdit)
    #     print("open group table")
    #
    # def setuptable(self,pd):
    #
    #     self.df=pd
    #     rown = len(self.df.index)
    #     coln = len(self.df.columns)
    #     self.model = QStandardItemModel(rown, 8)
    #     # labels = list(self.df.columns.values)
    #     # rown=len(self.df.index)
    #     # self.model = QStandardItemModel(rown,9)
    #     labels=['group','rep1','rep2','control','gene','strand','start','end']
    #     self.model.setHorizontalHeaderLabels(labels)
    #     # self.tableView.resize(500,300)
    #     #下面代码让表格100填满窗口
    #     self.tableView.horizontalHeader().setStretchLastSection(True)
    #     self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    #
    #     for row in range(rown):
    #         #print(self.df.loc[row].Sample)
    #         for column in range(8):
    #             item = QStandardItem(str(self.df.loc[row][labels[column]]))
    #             self.model.setItem(row, column, item)
    #
    #     self.tableView.setModel(self.model)
    #
    # def sampleEdit(self):
    #     self.close() ## 关闭窗口




if __name__ == "__main__":
    app =  QtWidgets.QApplication(sys.argv)
    window = showtable()
    window.show()
    sys.exit(app.exec_())