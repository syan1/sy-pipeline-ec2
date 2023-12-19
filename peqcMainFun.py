from rna_seq_tools.FileDealer import FileObj
from rna_seq_tools.QCpe import PairEndQualityController
import os
import multiprocessing
import sys
from numpy import integer
import argparse

def commonSubString(seq_file1, seq_file2):
    n = len(seq_file1)
    out = ""
    for i in range(n):
        if seq_file1[i] == seq_file2[i]:
            out += seq_file1[i]
        else:
            return(out)

def trimTruSeq2PE(expNum, seq_file, outDir):
    trimmer = PairEndQualityController()
    fArr1 = seq_file[0].split("/")
    fName1 =  fArr1[-1]
    
    fArr2 = seq_file[1].split("/")
    fName2 = fArr2[-1]
    
    readObj1 = FileObj(outDir, fName1)
    readObj2 = FileObj(outDir, fName2)

    base_file = commonSubString(fName1, fName2)
    base_file = base_file.strip("_")
    dir_name = base_file + "trimmed"
    out_path = os.path.join(outDir, dir_name)
    out_dir = FileObj(out_path, base_file)
    trimmer.TrimTruSeq2PE(read1 = readObj1, read2 = readObj2, bname = base_file, out_dir = out_dir)

def processFun(expNum, fs, outDir, fun):
    if fun == "trimpe":
        [trimTruSeq2PE(expNum, f, outDir) for f in fs]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-exp", "--experiment", type=str, help="experiment name", required = True)
    parser.add_argument("-p", "--process", type=int, help="number of process used during running", default = 1)
    parser.add_argument("-dir", "--directory", type=str, help="directory of seq raw data stored", default = "/data/cofactor_genomics/")
    parser.add_argument("-fun", "--function", type=str, help="toolname, star", required = True)

    try:
        args = parser.parse_args()
        expNum = args.experiment
        n = args.process
        dataDir = args.directory
        fun = args.function

        files = os.listdir(dataDir)
        files = [os.path.join(dataDir, f) for f in files]
        files = list(filter(lambda x: ".fastq.gz" in x, files))
        file_list = {}

        for i in range(n):
            file_list[i] = []

        print("Your input files are: ...")
        for x in sorted(files):
            print(x)

        isCorrect = input("Are files corrected? (Y/N)")
        if isCorrect.upper() == "Y":
            nfile = len(files)
            count = 0
            files = sorted(files)
            for i in range(0, nfile, 2):
                tmp = [files[i], files[i+1]]
                if count < n-1:
                    file_list[count].append(tmp)
                    count += 1
                else:
                    file_list[count].append(tmp)
                    count = 0
        for x in file_list:
            print("Process: %s uses these files: %s" % (x, file_list[x]) )
        isCorrect2 = input("Are the files assigned to process correctly? (Y/N)")
        if isCorrect2.upper() == "Y":
            for k in file_list:
                p = multiprocessing.Process(target = processFun, args = (expNum, file_list[k], dataDir, fun))
                p.start()
    except:
        print("there is an error")
        # Printer 1 for args
        try:
            args
        except NameError:
            print("args don't exist")
        else:
            print("these are the args:")
            print(args)
        # Printer 2 for file names
        try:
            files
            for x in files:
                print(x)
        except NameError:
            print("files don't exist")
        else:
            print("these are the files")
            for x in files:
                print(x)
        # Printer 3 for processed files
        print("these are the file_list items...")
        try:
            file_list
            for x in file_list:
                print(x)
        except NameError:
            print("file list doesn't exist")
        else:
            print("these are the file_list:")
            for x in file_list:
                print(file_list[x])
