'''
Aligner module - executes aligner on linux ec2 - assumes sudo and correct path
'''
from rna_seq_tools.FileDealer import *
import os
import logging


class Aligners(object):
    def __init__(self):
        pass
    
    def STAR(self, index_dir, read, gtf, out_dir, read_length):
        if not os.path.exists(out_dir.dir):
            os.system("sudo mkdir -p " + out_dir.dir )
             
        cmd = "sudo $aligners/STAR" + \
             " --runThreadN 12"  + \
             " --genomeDir " + index_dir + \
             " --readFilesIn " + read.path + \
             " --readFilesCommand gunzip -c" + \
             " --sjdbGTFfile " + gtf.path + \
             " --sjdbOverhang " + str(read_length) + \
             " --outFileNamePrefix " + out_dir.path + \
             " --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate" + \
             " --quantMode GeneCounts" 

        try:
            os.system(cmd)
            
        except:
            mss = self.file + " failed to executive"
            logging.warning(mss)
          
    def STARPE(self, index_dir, read1, read2, gtf, out_dir, read_length):
        if not os.path.exists(out_dir.dir):
            os.system("sudo mkdir -p " + out_dir.dir )
             
        cmd = "sudo $aligners/STAR" + \
             " --runThreadN 12"  + \
             " --genomeDir " + index_dir + \
             " --readFilesIn " + read1.path + " " + read2.path + \
             " --readFilesCommand gunzip -c" + \
             " --sjdbGTFfile " + gtf.path + \
             " --sjdbOverhang " + str(read_length) + \
             " --outFileNamePrefix " + out_dir.path + \
             " --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted" + \
             " --quantMode GeneCounts " + \
             " --limitBAMsortRAM 35760514860"
        try:
            os.system(cmd)
            
        except:
            mss = self.file + " failed to executive"
            logging.warning(mss)
    
    
      
    def tophat2(self, index_dir, read, gtf, out_dir):
        if not os.path.exists(out_dir.dir):
            os.system("mkdir -p " + out_dir.dir )
            
        cmd = "$aligners/tophat2" + \
            " -p 8" + \
            " -G " + gtf.path + \
            " -o " + out_dir.path + \
            ""
           
        pass
     
    def novoalign(self):
        cmd = "$aligners -d /path/to/novaindex -c multithreads(no working)" + \
            "--HLimit 5 " + \
            "-f seqfiles"

    
    def getAlignReport(dirs):
        pass
     
