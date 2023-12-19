
from rna_seq_tools.SeqFileParser import *
import os
import multiprocessing
import sys
import argparse
import re

def processFun(fs, outDir):
    for f in fs:
        #print(f.fileLink)
        downloadObj = DownLoadFileTool()
        downloadObj.downloadFileWget(f, outDir)
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-exp", "--experiment", type=str, help="experiment name", required = True)
    parser.add_argument("-f", "--filename", type=str, help="file only contain download link", required = True)
    parser.add_argument("-dir", "--directory", type=str, help="file download directory", default = "/data")
    parser.add_argument("-p", "--process", type=int, help="number of process used during running", default = 1)
    
    try:
        # get args from command line, first arg is experiment number
        # and the second args is number of process
        args = parser.parse_args()
        expNum = args.experiment
        n = args.process
        dataDir = args.directory
        fileName = args.filename
        
        with open(fileName, 'r') as f:
            links = f.readlines()
        
        downloadFiles = []   
        for l in links:
            l= l.strip()
            f = File()
            f.fileLink = l
            f.expNum = expNum
            m = re.search(r'^https.+Samples\/(.+)?\?.+', l)
            f.fileName = m.group(1)
            downloadFiles.append(f)
        
        
        #list input files 
        print("You input files are:......")
        print([f.fileName for f in downloadFiles])
        
        #check input files are corrected
        isCorrect = input("Are files corrected?(Y/N)")
        
        if isCorrect.upper() == "Y":  
            file_list = {}
            for i in range(n):
                file_list[i] = []
            count = 0
            for f in downloadFiles:
                if count < n-1:
                    file_list[count].append(f)
                    count += 1
                else:
                    file_list[count].append(f)
                    count = 0
            for k in file_list:
                p = multiprocessing.Process(target=processFun, args = (file_list[k], dataDir))
                p.start()
                
                
        else:
            print("function is stopped")
    
    except:
        print("Fail to run align function")

    

