from rna_seq_tools.FileDealer import *
from rna_seq_tools.Aligner import *
import os
import multiprocessing
import sys
from numpy import integer
import argparse

def starAlign(expNum, f,spe, outDir, readLength, ercc):
    aligner = Aligners()
    
    if spe == "human":
        index_dir = "/data/reference/hg38/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Homo_sapiens.GRCh38.89.gtf")
        
    elif spe == "rat":
        index_dir = "/data/reference/rat6/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Rattus_norvegicus.Rnor_6.0.90.gtf")
    elif spe == "mouse":
        index_dir = "/data/reference/mouse38/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Mus_musculus.GRCm38.90.gtf")

    if ercc.lower() == "yes":
        index_dir = index_dir + "_ercc"
        gtf.fileName = gtf.fileName + "_ercc"
        gtf._FileObj__getPath()
    
    fArr = f.split("/")
    fName = fArr[-1] 
    readObj = FileObj(outDir, fName)
    base_file = fName.replace(".fastq.gz", "")
    dir_name = base_file + "_aligned"
    out_path = os.path.join(outDir, "alignment/genome/", dir_name)
    out_dir = FileObj(out_path, base_file)
    aligner.STAR(index_dir, readObj, gtf, out_dir, readLength)
    print(readObj.path)
    print(out_dir.path)

def starPEAlign(expNum, seq_file, spe, outDir, readLength, ercc):
    aligner = Aligners()
    
    if spe == "human":
        index_dir = "/data/reference/hg38/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Homo_sapiens.GRCh38.89.gtf")
        
    elif spe == "rat":
        index_dir = "/data/reference/rat6/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Rattus_norvegicus.Rnor_6.0.90.gtf")
    
    elif spe == "mouse":
        index_dir = "/data/reference/mouse38/STAR_index"+ "_" + str(readLength)
        gtf = FileObj("/data/reference/annotation", "Mus_musculus.GRCm38.90.gtf")

    if ercc.lower() == "yes":
        index_dir = index_dir + "_ercc"
        gtf.fileName = gtf.fileName + "_ercc"
        gtf._FileObj__getPath()
    
    fArr1 = seq_file[0].split("/")
    fName1 =  fArr1[-1]
    
    fArr2 = seq_file[1].split("/")
    fName2 = fArr2[-1]
    
    readObj1 = FileObj(outDir, fName1)
    readObj2 = FileObj(outDir, fName2)
    
    base_file = commonSubString(fName1, fName2)
    dir_name = base_file + "aligned"
    out_path = os.path.join(outDir, "alignment/genome/", dir_name)
    out_dir = FileObj(out_path, base_file)
    aligner.STARPE(index_dir, readObj1,readObj2, gtf, out_dir, readLength)

def commonSubString(seq_file1, seq_file2):
    n = len(seq_file1)
    out = ""
    for i in range(n):
        if seq_file1[i] == seq_file2[i]:
            out += seq_file1[i]
        else:
            return(out)
    

def processFun(expNum, fs, spe, outDir, fun, pairend, readLength, ercc):
    if fun == "star" and pairend == "yes":
        [starPEAlign(expNum, f,spe, outDir, readLength, ercc) for f in fs]
    elif fun=="star" and pairend == "no":
        [starAlign(expNum, f,spe, outDir, readLength, ercc) for f in fs]


if __name__ == "__main__":
    #create option for qcMainFun Tool
    parser = argparse.ArgumentParser()
    parser.add_argument("-exp", "--experiment", type=str, help="experiment name", required = True)
    parser.add_argument("-spe", "--species", type=str, help="species: human, mouse, rat", required = True)
    parser.add_argument("-p", "--process", type=int, help="number of process used during running", default = 1)
    parser.add_argument("-dir", "--directory", type=str, help="directory of seq raw data stored", default = "/data")
    parser.add_argument("-fun", "--function", type=str, help="toolname, star", required = True)
    parser.add_argument("-pe", "--pairend", type=str, help="is pairend sequence",default="no")
    parser.add_argument("-l", "--readlength", type=int, help="maximum read length", required = True, default = 75)
    parser.add_argument("-e", "--ercc", type = str, help = "ercc genome", required = True, default = "no")

    try:
        # get args from commind line, first arg is experiment number
        # and the second args is number of process
        args = parser.parse_args()
        expNum = args.experiment
        n = args.process
        spe = args.species
        dataDir = args.directory
        fun = args.function
        pairend = args.pairend
        readLength = args.readlength
        ercc = args.ercc

        files = os.listdir(dataDir)
        files = [os.path.join(dataDir, f) for f in files]
        files = list(filter(lambda x: ".fastq.gz" in x, files))
        file_list = {}

        for i in range(n):
            file_list[i] = []

        #list input files
        print("You input files are:......")
        
        # 02/04/19 Edit. Make new line after each file name 
        for x in sorted(files):
            print(x)
        # print(sorted(files))
        #check input files are corrected
        isCorrect = input("Are files corrected?(Y/N)")

        if isCorrect.upper() == "Y":
            nfile = len(files)
            count = 0
            if pairend.upper() == "YES":
                files = sorted(files)
                for i in range(0,nfile,2):
                    tmp = [files[i], files[i+1]]
                    if count < n-1:
                        file_list[count].append(tmp)
                        count += 1
                    else:
                        file_list[count].append(tmp)
                        count = 0
            elif pairend.upper() == "NO":
                for f in files:
                    print(f)
                    if count < n-1:
                        file_list[count].append(f)
                        count += 1
                    else:
                        file_list[count].append(f)
                        count = 0
            for k in file_list:
                print(file_list[k])
                p = multiprocessing.Process(target=processFun, args = (expNum, file_list[k], spe, dataDir, fun, pairend, readLength, ercc))
                p.start()
        else:
            print("function is stopped")

    except:
        print("Fail to run align function")

