from rna_seq_tools.FileDealer import *
import os
import logging

class PairEndQualityController(object):
    def __init__(self):
        pass
    def TrimTruSeq2PE(self, read1, read2, bname, out_dir):
        if os.path.exists(out_dir.dir):
            os.rmdir(out_dir.dir)
            os.makedirs(out_dir.dir, mode = 0o777)
        elif not os.path.exists(out_dir.dir):
            os.makedirs(out_dir.dir, mode = 0o777)

        cmd = "trimmomatic PE" +\
            " -basein " + read1.path +\
            " -trimlog " + os.path.join(out_dir.dir, bname) + ".log" +\
            " -baseout " + os.path.join(out_dir.dir, bname) + ".fastq.gz" +\
            " ILLUMINACLIP:/data/anaconda3/envs/QC/share/trimmomatic-0.36-6/adapters/TruSeq3-PE.fa:2:30:10" +\
            " LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:40"
        try:
            os.system(cmd)
        except:
            mss = self.file + " failed to execute"
            logging.warning(mss)