import re
import requests
import os.path
import logging
import os


#Create a file object
class File(object):
    def __init__(self):
        self.expNum = None
        self.fileName = None
        self.fileLink = None


#Create a file parser, this is to get links in .txt files 
class LinkParser(object):
    #fileLoc is the mother directory of target link file
    #expNum is the name of the link file
    def __init__(self, fileLoc, expNum):
        self.files = []
        self.__getLinks(fileLoc, expNum)
        self.files = self.files[1:]
    def __getLinks(self, fileLoc, expNum):
        fLoc = os.path.join(fileLoc, expNum)
        with open(fLoc, "r") as linkFile:
            content = linkFile.readlines()
        content = [x.strip('\n') for x in content]
        link = ""
        file = File()
        for l in content:
            if self.__isFileName(l):
                file.fileLink = link
                self.files.append(file)
                file = File()
                link = ""
                file.fileName = l
                file.expNum = expNum        
            else:
                link += l
        
        file.fileLink = link
        self.files.append(file)
                
    def __isFileName(self, line):
        if "md5sums.txt"  in line and len(line) < 40:
            return True
        
        elif ".txt.gz" in line and len(line) < 40:
            return True
        
        elif ".fastq.gz" in line and len(line) < 40:
            return True
        
        else:
            return False
        
        
    def writeCSV(self, location, fileName):
        if len(self.files) == 0:
            print("No links to write")
            return
        
        file = os.path.join(location, fileName)
        with open(file, "a") as outCSV:
            for f in self.files:
                row = ",".join([f.expNum, f.fileName, f.fileLink])
                outCSV.write(row + "\n")
    
class DownLoadFileTool(object):
    #function to download file
    def downloadFile(self, File, dLocation):
        link = File.fileLink
        fileName= os.path.join(dLocation, File.fileName)
        try:
            r = requests.get(link)
        
        except:
            logging.warning(File.fileName+" fails")
                
        if r.status_code == requests.codes.ok:
            with open(fileName, "wb" ) as seqFile:
                seqFile.write(r.content)
        else:
            logging.warning(File.fileName + " " + str(r.status_code) + " fails")
    
    def downloadFileCurl(self, File, dLocation):
        cmd = 'sudo curl -o ' + os.path.join(dLocation, File.fileName) + " '" + File.fileLink + "' "
        try:
            #print(cmd)
            os.system(cmd)
        except:
            logging.warning(File.fileName+"fails")
    def downloadFileWget(self, File, dLocation):
        cmd = 'sudo wget -a logfile -O ' + os.path.join(dLocation, File.fileName) + " '" + File.fileLink + "' "
        try:
            os.system(cmd)
        except:
            logging.warning(File.fileName+"fails")
    
    def createFileObj(self, string):
        dirName, fileName, fileLink = string.split(',')
        dirName = dirName.strip("\.txt")
        res = File()
        res.expNum = dirName
        res.fileName = fileName
        res.fileLink = fileLink
        return(res)

    
def createLinksDocument(file, location):
    path = os.path.join(location, file) 
    
    with open(path, "rb") as r:
        content = r.readlines()
    
    content = [x.strip("\n") for x in content]
    
    for line in content:
        parser = LinkParser("RNA_seq", line) 
        parser.writeCSV("RNA_seq", "documents.csv")
