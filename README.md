# Pipeline for aws ec2 instance
Customized package for bulk RNA-seq processing pipeline to be used on outsourced sequencing projects.  

### Description  
This is a python package consisting of wrappers for command line functions that perform the following:  

* File download from vendor
* Fastq QC metric generation and filtration
* Alignment
* Read Count (quantification)

## Usage
The order of operation for the main scripts are:  

* downloadMainFun.py
* alignMainFun.py
* qcMainFun.py
* quantifierMainFun.py

### downloadMainFun.py  
```shell
python downloadMainFun.py --experiment IWEXP001 --filename download_links.txt --directory /data --process 2
```

### alignMainFun.py  
```shell
python alignMainFun.py --experiment IWEXP001 --species human --function star --pairend yes --readlength 100 --ercc no --process 4
```

### qcMainFun.py
```shell
python qcMainFun.py --experiment IWEXP001 --directory /data --function fastqc trimadaptor --process 2
```

### quantifierMainFun.py
```shell
python quantifierMainFun.py --experiment IWEXP001 --species human --process 2
```