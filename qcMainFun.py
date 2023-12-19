from rna_seq_tools.FileDealer import *
from rna_seq_tools.QualityController import *
import os
import multiprocessing
import sys
import argparse

def fastqc(expNum, f, outDir):
    data_path = f
    qcResult_path = os.path.join(outDir,"qcTest")
    qcControl = QualityController()
    qcControl.FastQC(data_path, qcResult_path)

def fastqcTrimmed(expNum, f, outDir):
    #view quality of trimmed seq data 
    
    data_path = f
    qcResult_path = os.path.join(outDir,"qcTestTrimmed")
    qcControl = QualityController()
    qcControl.FastQC(data_path, qcResult_path)

def getFastqcReport(expNum):
    #view quality of original seq data 
    data_root = '/data'
    fastqc_root = '/data/qcTest'
    seqFiles = getSeqFiles(data_root, expNum)
    fastqcReport(fastqc_root, expNum, seqFiles) 

def trimAdaptor(expNum, file, outDir):
    # trim adpator 
    outroot = os.path.join(outDir, "trimmed")
    encoding = "33"
    infile = file
    trimQC = SingleEndQualityController(outroot, encoding, infile)
    trimQC.TrimmomaticAdaptor()

def trimLowScore(expNum, file, outDir):
    # trim low score base 
    outroot = os.path.join(outDir,"trimmed")
    encoding = "33"
    infile = file
    trimQC = SingleEndQualityController(outroot, encoding, infile)
    trimQC.TrimmomaticLowScore()

def trimLowScoreAdpator(expNum, file, outDir):
    # trim low Score and adaptor
    outroot = os.path.join(outDir, "trimmed")
    encoding = "33"
    infile = file
    trimQC = SingleEndQualityController(outroot, encoding, infile)
    trimQC.TrimmomaticLowScoreAdaptor()

def trimAdaptorTruSeq2(expNum, file, outDir):
    # trim adpator 
    outroot = os.path.join(outDir, "trimmed")
    encoding = "33"
    infile = file
    trimQC = SingleEndQualityController(outroot, encoding, infile)
    trimQC.TrimmomaticAdaptorTruSeq2()    

def filterrRNA(expNum, f, outDir):
    #to filter out rRNA from raw reads
    qc = QualityController()
    ref_root = "/data/reference/rRNA"
    refs = ["rfam_5s_id98.fasta", "rfam_5.8s_id98.fasta", "silva_euk_28s_id98.fasta", "silva_euk_18s_id95.fasta"]
    refObjs = []
    for r in refs:
        file = FileObj(ref_root, r)
        refObjs.append(file)
      
    index_root = "/data/reference/rRNA/sortmeRNA_index"
    indexs = ["rfam_5s_id98_db", "rfam_5.8s_id98_db", "silva_euk_28s_id98_db", "silva_euk_18s_id95_db"]
    indexObjs = []
    
    for i in indexs:
        file = FileObj(index_root, i)
        indexObjs.append(file)
    
    fArr = f.split("/")
    fName = fArr[-1]
    readObj = FileObj(outDir, fName)
    print(fName, outDir)
    print(readObj.path)
    print(readObj.fileName)

    qc.sortmerna(refObjs, indexObjs, readObj)
 
def processFun(expNum, fs, outDir, fun):
    if fun == "fastqc":
        [fastqc(expNum, f, outDir) for f in fs]
    elif fun == "fastqctrimmed":
        [fastqcTrimmed(expNum, f, outDir) for f in fs]
    elif fun == "sortrna":
        [filterrRNA(expNum, f, outDir) for f in fs]
    elif fun == "trimadaptor":
        [trimAdaptor(expNum, f, outDir) for f in fs]
    elif fun == "trimlowscore":
        [trimLowScore(expNum, f, outDir) for f in fs]
    elif fun == "trimlowscoreadaptor":
        [trimLowScoreAdaptor(expNum, f, outDir) for f in fs]
    elif fun == "trimadaptor2":
        [trimAdaptorTruSeq2(expNum, f, outDir) for f in fs]
        

if __name__ == "__main__" :
    #create option for qcMainFun Tool
    parser = argparse.ArgumentParser()
    parser.add_argument("-exp", "--experiment", type=str, help="experiment name", required = True)
    parser.add_argument("-p", "--process", type=int, help="number of process used during running", default = 1)
    parser.add_argument("-dir", "--directory", type=str, help="directory of seq raw data stored", default = "/data")
    parser.add_argument("-fun", "--function", type=str, help="toolname, fastqc, fastqctrimmed, sortrna, trimadaptor, trimlowscore, trimlowscoreadaptor, trimadaptor2", required = True)
    
    try:
        args = parser.parse_args()
        expNum = args.experiment
        n = args.process
        dataDir = args.directory
        fun = args.function

        #trimmed files  
        files = os.listdir(dataDir)
        files = [os.path.join(dataDir, f) for f in files]

        file_list = {}
        
        for i in range(n):
            file_list[i] = []
        
        #list input files 
        print("You input files are:......")
        
        for x in files:
            print(x)
        
        #check input files are corrected
        isCorrect = input("Are files corrected?(Y/N)")

        if isCorrect.upper() == "Y":      
            count = 0
            for f in files:
                if ".gz" in f:
                    if count < n-1:
                        file_list[count].append(f)
                        count += 1
                    else:
                        file_list[count].append(f)
                        count = 0
            for k in file_list:
                p = multiprocessing.Process(target = processFun, args=(expNum, file_list[k],dataDir, fun))    
                p.start()
        else:
            print("function is stopped")
    except:
        print("Fail to run quality controller")
