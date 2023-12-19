import os
import logging
import re
import numpy as np
import pandas as pd
import tempfile
import shutil
import errno
from rna_seq_tools.FileDealer import FileObj
import pdb

# pdb.set_trace()
class BasicStat(object):
    def __init__(self):
        self.fileName = None
        self.fileType = None
        self.encoding = None
        self.totalSeq = 0
        self.flaggedPoor = 0
        self.seqLen = 0
        self.gcContent = 0
        
class OverRepresentSeq(object):
    def __init__(self):
        self.data = pd.DataFrame({'Sequence':[], 
                                  'Count':[],
                                  'Percentage':[],
                                  'Possible Source':[]})
   
class QualityController(object):
    def __init__(self, file=None):
        self.file = file
        self.passed = []
        self.fails = []
        self.itemDict = {'Basic Statistics': 'basic_statistics', 
                         'Per base sequence quality': 'per_base_quality', 
                         'Per tile sequence quality': 'per_tile_quality', 
                         'Per sequence quality scores': 'per_sequence_quality',
                         'Per base sequence content': 'per_base_sequence_content', 
                         'Per sequence GC content': 'per_sequence_gc_content',
                         'Per base N content': 'per_base_n_content',
                         'Sequence Length Distribution': 'sequence_length_distribution',
                         'Sequence Duplication Levels': 'duplication_levels',
                         'Overrepresented sequences': 'overrepresented_sequences', 
                         'Adapter Content': 'adapter_content', 
                         'Kmer Content': 'kmer_profile'}
        
    def FastQC(self, data , outLocation = None):
        if outLocation:
            if not os.path.exists(outLocation):
                os.system("sudo mkdir -p " + outLocation)
                #os.makedirs(outLocation)
            outLocation = "--outdir=" + outLocation
        
        # 02/04/19 Edit: corrected path to be more general, don't use sudo + specific/path     
        cmd = "fastqc --extract " + data + " " + outLocation  
        
        try:
            os.system(cmd)
        except:
            mss = self.file + " failed to executive"
            logging.warning(mss)
            
    def parseFastQC(self, directory):
        summaryLoc = os.path.join(directory, "summary.txt")
        #get summary infomation
        with open(summaryLoc) as file:
            summary = file.readlines()
            
        summary = [x.strip("\n").split("\t") for x in summary]

        outDir = os.path.join(directory, self.file+".parsedOut")
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        for item in summary:
            if item[0] == "PASS":
                self.passed.append(item)
            else:
                self.fails.append(item)
            
            if item[1] == 'Basic Statistics':
                basicStatObj = BasicStat()
                self._fillInStat(directory, basicStatObj)
                
    def __fillInStat(self, directory, basicStatObj):
        dataLoc = os.path.join(directory, "fastqc_data")
        with open(dataLoc, "r") as r:
            content = r.read()
            find = re.findall(r">>Basic\sStatistics(.*?)>>END", content, re.S)[0]
    
            lines = find.split("\n")[2:]
            lines = [x.split("\t") for x in lines]
                
        #create two files to record the pass items and fails items 
#         self.__writePass(outDir)
#         self.__writeFails(outDir)
        
    def __writePass(self, outDir):
        fileName = os.path.join(outDir, self.file + ".passed")  
        with open(fileName, "a") as out:
            [out.write(i) for i in self.passed]
            out.write("total passed: " + str(len(self.passed)))
    
    def __writeFails(self, outDir, fails):
        fileName = os.path.join(outDir, self.file + ".passed")  
        with open(fileName, "a") as out:
            [out.write(i) for i in self.fails]
            out.write("total fail: " + str(len(self.fails)))
    
    #from here use fileObj, code above need to modified
    def sortmerna(self, refObjs, indexObjs, readObj):
        try:
            tmp_dir = tempfile.mkdtemp(dir='/data/tmp/')  # create dir
            # decompress fastq.gz file and save it into tempDir 
            temp_name = readObj.fileName.strip(".gz")
            tempFile = FileObj(tmp_dir, temp_name)
            unzip_cmd = "gunzip -c " + readObj.path+ " >" + tempFile.path
            
            ref = ""
            for (r, i) in zip(refObjs, indexObjs):
                ref += r.path + "," + i.path + ":"
            ref = ref.strip(":") 
            
            rRNA_dir = os.path.join(readObj.dir, "rRNA_filtered")
            if not os.path.exists(rRNA_dir):
                os.system("mkdir " + rRNA_dir )
            if ".fastq.gz" in readObj.fileName:
                rrna_name = readObj.fileName.strip(".fastq.gz")
                non_rrna_name = readObj.fileName.strip(".fastq.gz")
            elif ".txt.gz" in readObj.fileName:
                rrna_name = readObj.fileName.strip(".txt.gz")
                non_rrna_name = readObj.fileName.strip(".txt.gz")
                
            rrna_name += "_rRNA"
            rrnaObj = FileObj(rRNA_dir, rrna_name)
            non_rrna_name += "_non_rRNA"
            non_rrnaObj = FileObj(rRNA_dir, non_rrna_name)
            
            # 02/04/19 Edit: corrected path to be more general, don't use sudo + specific/path
            # Changed to increase thread -a 4 instead of -a 2.
            # 02/06/19 Edit: changed thread back to 2 because more threads did not run properly. ie shorter files
            # 02/06/19 Edit: Also added .fastq file extension to sortmerna as currently files are written without and that breaks zip_cmd     
            # 02/08/19 Edit: Reverted .fastq addition.
            # 02/08/19 Edit: Increased thread count and memory count. CAUTION: Can only run with p1 at once since it takes 24 for 1 processes now.
            cmd = "sortmerna --ref " + ref \
                  + " --reads " + tempFile.path \
                  + " --num_alignments 1 -a 12 -m 12288 --fastx " \
                  + "--aligned " + rrnaObj.path \
                  + " --other " + non_rrnaObj.path + " --log -v" 
                  
            zip_cmd = "gzip " + rrnaObj.path + ".fastq " + non_rrnaObj.path + ".fastq" 
        
            try:
                os.system(unzip_cmd)
                os.system(cmd)
                os.system(zip_cmd)
            except:
                mss = readObj.fileName + " failed to executive"
                logging.warning(mss)
     
        finally:
            try:
                shutil.rmtree(tmp_dir)  # delete directory
            except OSError as exc:
                if exc.errno != errno.ENOENT:  # ENOENT - no such file or directory
                    raise  # re-raise exception



class AlignmentControl(object):
    def coverageandGC(self, samFiles, out_dir, expNum):
        try:
            tmp_dir = tempfile.mkdtemp()  # create dir
            sortedBamFiles = []
            for samfile in samFiles:
                self.__sam2bam(samfile, tmp_dir, sortedBamFiles)
            
            out_prefix = os.path.join(out_dir, expNum)
            i = ",".join([x.path for x in sortedBamFiles])
            coverage_cmd = "geneBody_coverage.py" + \
                            " -r /data/cofactor_genomics/reference/annotation/Homo_sapiens.GRCh38.89.bed" + \
                            " -i " + i + \
                            " -o " + out_prefix
                           
            try:
                os.system(coverage_cmd)
            except:
                mss ="coverage failed to executive"
                logging.warning(mss)
            
            for bam in sortedBamFiles:
                outGC = os.path.join(out_dir, bam.fileName) 
                gc_cmd = "read_GC.py -i " + bam.path + " -o " + outGC
            
        finally:
            try:
                shutil.rmtree(tmp_dir)  # delete directory
            except OSError as exc:
                if exc.errno != errno.ENOENT:  # ENOENT - no such file or directory
                    raise  # re-raise exception


        #end here, below code is also need to modified using FileObj
               

class SingleEndQualityController(QualityController):
    def __init__(self, outroot, encoding, infile):
        if ".fastq.gz" in infile:
            filename = infile.strip('.fastq.gz')
        elif ".txt.gz" in infile:
            filename = infile.strip('.txt.gz')
        else:
            filename = infile
        outfile = filename + "_trimmed.fastq.gz"
        self.outfile = os.path.join(outroot, outfile)
        
        logfile = filename + "_log"
        self.logfile = os.path.join(outroot, logfile) 
        self.encoding = encoding
        self.infile = infile

    def TrimmomaticAdaptor(self):         
        cmd = "trimmomatic SE " + self.infile + " " \
              + self.outfile + " -phred" + self.encoding \
              + ' -trimlog ' + self.logfile \
              + " ILLUMINACLIP:/data/anaconda3/envs/QC/share/trimmomatic-0.36-6/adapters/TruSeq3-SE.fa:2:30:10"
        try:
            os.system(cmd)
        except:
            mss = self.infile + " failed to executive"
            logging.warning(mss)
        
    def TrimmomaticLowScore(self):
        cmd = "trimmomatic SE " + self.infile + " " \
              + self.outfile + " -phred" + self.encoding \
              + ' -trimlog ' + self.logfile \
              + ' SLIDINGWINDOW:5:20' \
              + ' SLIDINGWINDOW:5:20' 
        try:
            os.system(cmd)
        except:
            mss = self.infile + " failed to executive"
            logging.warning(mss)

    def TrimmomaticLowScoreAdaptor(self): 
        cmd = "trimmomatic SE " + self.infile + " " \
              + self.outfile + " -phred" + self.encoding \
              + ' -trimlog ' + self.logfile \
              + ' SLIDINGWINDOW:5:20' \
              + " ILLUMINACLIP:/home/ubuntu/anaconda3/envs/QC/share/trimmomatic-0.36-6/adapters/TruSeq3-SE.fa:2:30:10"     
        try:
            os.system(cmd)
        except:
            mss = self.infile + " failed to executive"
            logging.warning(mss)

    def TrimmomaticAdaptorTruSeq2(self):         
        cmd = "trimmomatic SE " + self.infile + " " \
              + self.outfile + " -phred" + self.encoding \
              + ' -trimlog ' + self.logfile \
              + " ILLUMINACLIP:/data/anaconda3/envs/QC/share/trimmomatic-0.36-6/adapters/TruSeq2-SE.fa:2:30:10"
        try:
            os.system(cmd)
        except:
            mss = self.infile + " failed to executive"
            logging.warning(mss)
        

class PairEndQualityController(QualityController):
    def Trimmomatic(self):
        pass
    
    


def fastqcReport(root_path, exp, samples): 
    report_name = exp + "_fastqc_report.xlsx"
    save_path = os.path.join(root_path, exp, report_name)
    report = {}
       
    for s in samples:
        s2 = s.split(".")[0]
        s2 += "_fastqc"
        path = os.path.join(root_path, exp, s2, "summary.txt")

        with open(path) as file:
            summary = file.readlines()
        summary = [x.strip("\n").split("\t") for x in summary]
        report[s] = []
        for item in summary:
            if item[0] == "PASS":
                report[s].append("PASS")
            elif item[0] == "WARN":
                report[s].append("WARN")
            elif item[0] == "FAIL":
                report[s].append("FAIL")
            else:
                report[s].append(item[0])
    index = ['Basic Statistics', 
             'Per base sequence quality', 
             'Per tile sequence quality', 
             'Per sequence quality scores',
             'Per base sequence content', 
             'Per sequence GC content',
             'Per base N content',
             'Sequence Length Distribution',
             'Sequence Duplication Levels',
             'Overrepresented sequences', 
             'Adapter Content', 
             'Kmer Content']            
    df = pd.DataFrame(report, index = index)
    
    #write to excel 
    writer = pd.ExcelWriter(save_path)
    df.to_excel(writer,"Sheet1")
    writer.save()        