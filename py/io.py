import re

def rep_copy(replist, fn0, fn1):
    """
    copy file 'fn0' to 'fn1' with replace indicated by 'replist'
    """
    with open(fn0, "r") as fi:
        with open(fn1, "w") as fo:
            for line in fi.readlines():
                for (k0, k1) in replist:
                    line = re.sub(k0, str(k1), line)
                fo.write(line)

def extract_lines(key0, key1, fn0, fn1, mode, n=1):
    """
    Extract lines from key0 to key1 in file 'fn0' and
    paste them to file 'fn1'.
    If mode=="a", lines are added to file 'fn1'.

    key0 : string
    .      keyword which means start of extraction
    key1 : string
    .      keyword which means end of extraction
    fn0  : string
    .      input file name
    fn1  : string 
    .      output file name
    mode : string
    .     "w" for rewrite or "a" for append
    n : integer
    .     n th lines are extracted
    """

    with open(fn1, mode) as f1:
        with open(fn0, "r") as f0:
            nfind = 0
            status = 0
            for line in f0.readlines():
                if(status==0 and line.find(key0)>-1):
                    nfind += 1
                if(status==0 and nfind==n):
                    status=1
                if(status==1):
                    f1.write(line)
                if(status==1 and line.find(key1)>-1):
                    status=2


                    
