import sys
import pandas as pd
import matplotlib.pyplot as plt

from PyQt5 import uic,QtWidgets
from CMlib.show_barcodestable import showtable as showbarcode
from multiprocessing import Pool
from functools import  partial
from CMlib.split_fastq import split_fastq
from PyQt5.QtWidgets import QHeaderView, QPushButton,QProgressDialog
from PyQt5.QtCore import Qt, QBasicTimer
import os

from subprocess import Popen, PIPE

path = os.getcwd()
qtCreatorFile = os.path.join(path,'CMlib/split_lanes.ui')

Ui_showtable, QtBaseClass = uic.loadUiType(qtCreatorFile)

class showtable(QtWidgets.QDialog, Ui_showtable):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        Ui_showtable.__init__(self)
        self.setupUi(self)

        self.setWindowTitle('Split FastQ')
        self.resize(500,400)
        # self.fastqbtn.clicked.connect(lambda: self.getfastq("left"))
        self.fastqbtn.clicked.connect(self.getfastq)
        self.fastqline.setReadOnly(True)  ##设置不可输入
        # self.rightbtn.clicked.connect(lambda: self.getfastq("right"))
        # self.right.setReadOnly(True)
        self.barcodebtn.clicked.connect(self.barcodeinfo)
        self.barcodeline.setReadOnly(True)
        self.outputbtn.clicked.connect(self.outputdir)
        self.outputline.setReadOnly(True)

        self.showbtn.clicked.connect(self.showtable)
        self.splitbtn.clicked.connect(self.split)
        self.resetbtn.clicked.connect(self.reset)

        self.path1 = ""
        self.path2 = ""
        self.path1check = ""
        self.edit = ""


    def getfastq(self):
        fastqPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', path)
        if fastqPath != "":
            print("fastq Direction", fastqPath)
            self.fastqline.setText(fastqPath)
            # self.fastq = fastqPath
            self.path2 = fastqPath




    def barcodeinfo(self):
        barcodepath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', path)
        if barcodepath != "":
            print("Direction", barcodepath)
            self.barcodeline.setText(barcodepath)
            self.dfbarcode = pd.read_csv(str(barcodepath))
            self.path1 = barcodepath

    def outputdir(self):
        outputdirpath = QtWidgets.QFileDialog.getExistingDirectory(self, 'open directory', path)
        if outputdirpath != "":
            print("Direction", outputdirpath)
            self.outputdirpath = outputdirpath
            self.outputline.setText(outputdirpath)
            self.outputfiledir = outputdirpath


    def showtable(self):
        if self.path1 !="":
            self.ui = showbarcode()  ##打开showtable新窗口
            result = self.ui.setuptable(self.dfbarcode)  ##传递倒入sample csv

            if result == "yes":
                self.ui.show()  ##显示窗
                # self.newdfbarcode, self.edit = self.ui.sampleEdit()

                # print(newdf)
                self.path1check = self.path1
            else:
                self.path1check = ""
                self.path1 = ""
        else:
            self.showMessageBox('Warning', 'Please load Sample information Table first')
            self.path1check=""

    def reset(self):
        self.path1 = ""
        self.path2 = ""
        self.path1check = ""
        self.edit = ""


    def split(self):
        if self.path2 != "":
            if self.path1check != "":
                self.edit, self.newdfbarcode = self.ui.resulttest()  ##check table has fixed
                if self.edit == "yes":
                    self.showbarprocess("Prepare for splitting...")
                    # self.figures = Example()
                    # self.figures.initUI()
                    # self.figures.show()

                    pool = Pool(4)
                    pool.map(partial(split_fastq, df=self.newdfbarcode, fastq=self.path2, output=self.outputfiledir),
                             list(self.newdfbarcode.index))
                    # pool.map(partial(split_fastq,df=self.dfbarcode,fastq=self.path2,output=self.outputfiledir),list(self.dfbarcode.index))
                    print('done')
                    # self.showMessageBox('Warning', 'Please click show buttons for information checking')
                    self.showfinishBox('Information','The project have been done!')
                else:
                    self.showMessageBox('Warning', 'Please click "Confirm" button for barcodes checking')
                    self.edit = ""


            else:
                self.showMessageBox('Warning', 'Please click show buttons for information checking and confirming')
                self.path1check = ""
        else:
            self.showMessageBox('Warning', 'Please load fastq file first')
            self.path2 = ""



    # ############## warning message #########
    def showMessageBox(self, title, message):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Warning)
        msgBox.setWindowTitle(title)
        msgBox.setText(message)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.exec_()
    ##################################################

    def showfinishBox(self, title, message):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setWindowTitle(title)
        msgBox.setIcon(QtWidgets.QMessageBox.Information)
        msgBox.setText(message)
        msgBox.setDetailedText("The project has finished, please check the result!")
        msgBox.exec_()

    ################################################

    def showbarprocess(self,content):

        num = int(100000)
        progress = QProgressDialog(parent=self)
        progress.setWindowTitle("Start Processing ...")
        progress.setLabelText(content)
        # progress.setCancelButtonText("0")
        progress.setCancelButton(None) ##不显示cancel button
        progress.setMinimumDuration(5)
        progress.setWindowModality(Qt.WindowModal)
        progress.setRange(0, num)

        for i in range(num):
            progress.setValue(i)
        else:
            progress.setValue(num)

        progress.cancel()  ##直接关闭


if __name__ == "__main__":
    app =  QtWidgets.QApplication(sys.argv)
    window = showtable()
    window.show()
    sys.exit(app.exec_())