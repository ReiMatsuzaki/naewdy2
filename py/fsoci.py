import re
import numpy as np

class CIState:
    def __init__(self):
        self.istate = 0
        self.al = []
        self.be = []
        self.ci = []
        self.ene = None
        self.s = None
    def num(self):
        return len(self.ci)
    
    def __str__(self):
        res = "==== CI State ====\n"
        res += "istate = {0}\n".format(self.istate)
        res += "E      = {0}\n".format(self.ene)
        res += "S      = {0}\n".format(self.s)
        res += " alpha | beta | coef \n"
        for i in range(self.num()):
            res += str(self.al[i]) + " | "
            res += str(self.be[i]) + " | "
            res += str(self.ci[i]) + "\n"
        
        return res
        
def get_string(line, ncore):
    res = []
    ss = line.split()
    for i in (range(len(ss[0]))):
        s = ss[0][i]
        if(s=='1'):
            res.append(int(1+i+ncore))
    if(len(ss)==2):
        res.append(int(ss[1]))
    return res

def read_cicoef(f, ncore, cistate):
    print "start reading CI coefficient"
    line = f.readline()
    line = f.readline()
    while line!="\n":
        [al,be,ci] = line.split("|")
        cistate.al.append(get_string(al, ncore))
        cistate.be.append(get_string(be, ncore))
        cistate.ci.append(float(ci))
        line = f.readline()

def get_ncore(line):
    ss = [s for s in line.split()
          if s.find("CORE=")>-1]
    return int(ss[0].split("=")[-1])
    
def CIState_read(fn, sval=None):

    re_ncore = re.compile("NCORE=")
    
    re_begin = re.compile(r" ALL STATES CONVERGED.")
    re_state = re.compile(r" STATE   ([0-9]+)  ENERGY=      ([-0-9\.]+)  S=  ([-0-9\.]+)")
    re_cicoef = re.compile(r"     ALPHA      |     BETA       | CI COEFFICIENT")

    ncore = None
    istate = None
    mode = "none"
    cistates = []
    
    with open(fn, "r") as f:
        line = f.readline()    
        while line:
        
            line = f.readline()

            res = re_ncore.search(line)
            if(res):
                print "ncore"
                ncore = get_ncore(line)

            res=re_begin.search(line)
            if(mode=="none" and res):
                print "none->find"
                mode = "head"

            res=re_state.search(line)
            if(mode=="head" and res):
                print "head->value"
                istate = int(res.group(1))
                ene = res.group(2)
                s = res.group(3)
                cistates.append(CIState())
                cistates[istate-1].istate = istate
                cistates[istate-1].ene = float(ene)
                cistates[istate-1].s = float(s)
                cistates[istate-1].ene = ene
                mode = "value"
                
            res=re_cicoef.search(line)                
            if(mode=="value" and res):
                print "value->head"
                print "istate = ", istate
                read_cicoef(f, ncore, cistates[istate-1])
                mode = "head"

    res = []
    if(sval is None):
        res = cistates
    else:
        for ci in cistates:
            if(abs(ci.s-sval)<0.0001):
                res.append(ci)
    return res


def eij(alI, beI, i, j, alJ, beJ):
    return 1/np.sqrt(2.0) * (aiaj(alI, i, j, alJ) + aiaj(beI, i, j, beJ))

def aiaj(stringI, i, j, stringJ):
    if(i not in stringI or j not in stringJ):
        return 0
    
    sign = sign_ai(stringI, i) * sign_ai(stringJ, j)

    ai_stringI = stringI[:]
    ai_stringI.remove(i)

    aj_stringJ = stringJ[:]
    aj_stringJ.remove(j)

    if(ai_stringI==aj_stringJ):
        
        return sign
    else:
        return 0
    
def sign_ai(string, i):
    if(i not in string):
        return 0
    else:
        return (-1)**(len(string)-1-string.index(i))
