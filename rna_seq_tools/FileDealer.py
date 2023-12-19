import os
import logging

class FileObj(object):
    def __init__(self, directory=None, name=None):
        self.dir = directory
        self.fileName = name
        self.path = None
        self.__getPath()
        
    def __getPath(self):
        if self.dir and self.fileName:
            self.path = os.path.join(self.dir, self.fileName)
        else:
            self.path = None

def getSeqFiles(root_path, exp):
    md5Loc = os.path.join(root_path, exp, "md5sums.txt")
    try:
        md5sums = open(md5Loc, "r")
    except:
        logging.warning("no seq files in this directory")
    
    contents = md5sums.readlines()
    md5sums.close()
    
    contents = [x.strip('\n').split('  ') for x in contents]
    seqFiles = []
    for x in contents:
        if len(x) > 1:
            seqFiles.append(x[1])  
    return(seqFiles)
    
def makeFastDirs(root_path, exp, sample):
    raw_path = [root_path, exp, sample, "raw"]
    trim_path = [root_path, exp, sample, "trimmed"]
    
    raw_path = os.path.join(*raw_path)
    trim_path = os.path.join(*trim_path)
    
    if not os.path.exists(raw_path):
        os.makedirs(raw_path)
    
    if not os.path.exists(trim_path):
        os.makedirs(trim_path)
            
def sam2bam(samFile):
    if ".sam.gz" in samFile.fileName:
    # decompress fastq.gz file and save it into tempDir 
        base_name = samFile.fileName.strip(".sam.gz")
        unzip_cmd = "sudo gzip -d " + samFile.path
        print(unzip_cmd)
        try:
            os.system(unzip_cmd)
        except:
            mss = samFile.fileName + " failed to executive"
            logging.warning(mss)
    else:
        base_name = samFile.fileName.strip(".sam")
    #convert sam to sorted bam file
    tempFile = FileObj(samFile.dir, base_name+".sam")
    bamFile = FileObj(samFile.dir, base_name+ ".bam")
    sam2bam = "sudo sh -c \"$tools/samtools view -S -b " + tempFile.path + " > " + bamFile.path + "\"" 
    rm_cmd = "sudo rm " + tempFile.path
    print(rm_cmd, sam2bam)
    try:
        os.system(sam2bam)
        os.system(rm_cmd)
    except:
        mss = samFile.fileName + " failed to executive"
        logging.warning(mss)
         
def sortbam(bamFile):
    base_name = bamFile.fileName.strip(".bam")
    sortbamFile = FileObj(bamFile.dir, base_name+ ".sorted.bam")
    sortbam = "sudo $tools/samtools sort " + bamFile.path + " -o " + sortbamFile.path
    try:
        os.system(sortbam)
    except:
        mss = bamFile.fileName + " failed to executive"
        logging.warning(mss)

def getColumns(cols, common, output):
    cmd = "for f in *.common; do cut -f " + ",".join(cols) + "$f > " + output + "; done"



    
    
