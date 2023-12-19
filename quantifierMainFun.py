from rna_seq_tools.Quantifier import *
from rna_seq_tools.FileDealer import *
import multiprocessing
import os
import argparse

def countFun(b, gtfFile, exp):
    q = Quantifier()
    name = b.fileName.replace("Aligned.sortedByCoord.out.bam", "")
    name += "_count.txt"
    directory = '/data' + exp + '/htseq/'
    outFile = FileObj(directory=directory, name=name)
    q.countRead(b, gtfFile, outFile)

def countPEFun(b, gtfFile, exp):
    q = Quantifier()
    name = b.fileName.replace("Aligned.sortedByCoord.out.bam", "")
    name += "_count.txt"
    directory = '/data' + exp + '/htseq/'
    outFile = FileObj(directory=directory, name=name)
    q.countReadPE(b, gtfFile, outFile)

def processFun(fs, gtfFile, exp):
    for b in fs:
        print(b.fileName)
        #countFun(b, gtfFile, exp)
        countPEFun(b, gtfFile, exp)

if __name__ == "__main__" :
    #create option for qcMainFun Tool
    parser = argparse.ArgumentParser()
    parser.add_argument("-exp", "--experiment", type=str, help="experiment name", required = True)
    parser.add_argument("-p", "--process", type=int, help="number of process used during running", default = 1)
    parser.add_argument("-s", "--species", type=str, help="species sequenced", default = "human")
    
    try:
        args = parser.parse_args()
        exp = args.experiment
        n = args.process
        spe = args.species
    except:
        print("Failed to run quantifier")

    dirs = os.listdir(os.path.join("/data" ,exp, "/alignment/genome/"))
    
    bamroot = os.path.join( "/data", exp, "/alignment/genome/")
    bamFiles = []
    
    if spe == "human":
        gtfFile = FileObj(directory="/data/reference/annotation/", name = "Homo_sapiens.GRCh38.89.gtf")
    
    elif spe == "rat":
        gtfFile = FileObj(directory="/data/reference/annotation/", name = "Rattus_norvegicus.Rnor_6.0.90.gtf")
    
    elif spe == "mouse":
        gtfFile = FileObj(directory="/data/reference/annotation/", name = "Mus_musculus.GRCm38.90.gtf")
        
    for d in dirs:
        path = bamroot + d
        f = d.replace("_aligned", "")
        fileName = f + "Aligned.sortedByCoord.out.bam"
        fileObj = FileObj(path, fileName)
        bamFiles.append(fileObj)

    file_list = {}
    for i in range(n):
        file_list[i] = []
        
    count = 0
    for b in bamFiles:
        if count <n-1:
            file_list[count].append(b)
            count += 1
        else:
            file_list[count].append(b)
            count = 0
    
    for fs in file_list:
        p = multiprocessing.Process(target = processFun, args = (file_list[fs], gtfFile, exp))
        p.start()
