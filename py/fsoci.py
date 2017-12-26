import re
import numpy as np

class CIState:
    def __init__(self):
        self.istate = 0
        self.istate_s = 0
        self.al = []
        self.be = []
        self.ci = []
        self.ene = None
        self.s = None
        self.e_dmom = []
    def num(self):
        return len(self.ci)
    
    def __str__(self):
        res = "==== CI State ====\n"
        res += "istate = {0}\n".format(self.istate)
        res += "E      = {0}\n".format(self.ene)
        res += "S      = {0}\n".format(self.s)
        res += "e_dmom = {0}, {1}, {2}\n".format(self.e_dmom[0],
                                                 self.e_dmom[1],
                                                 self.e_dmom[2])
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
    
def CIState_read(fn, sval=0.0):

    re_ncore = re.compile("NCORE=")
    
    re_begin = re.compile(r" ALL STATES CONVERGED.")
    re_state = re.compile(r" STATE   ([0-9]+)  ENERGY=      ([-0-9\.]+)  S=  ([-0-9\.]+)")
    re_cicoef = re.compile(r"     ALPHA      |     BETA       | CI COEFFICIENT")

    re_prop_begin = re.compile(r"FSOCI CI PROPERTIES...FOR THE WAVEFUNCTION OF STATE    ([0-9]+)")
    re_e_dmom = re.compile(r"ELECTROSTATIC MOMENTS")

    ncore = None
    istate_s = 0
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
                ene = float(res.group(2))
                s = float(res.group(3))
                if(abs(s-sval)<0.0001):
                    istate_s += 1
                    ci = CIState()
                    ci.istate = istate
                    ci.istate_s = istate_s
                    ci.ene = float(ene)
                    ci.s = float(s)
                    ci.ene = ene
                    cistates.append(ci)
                    mode = "value"
                
            res=re_cicoef.search(line)                
            if(mode=="value" and res):
                print "value->head"
                print "istate = ", istate
                read_cicoef(f, ncore, cistates[istate_s-1])
                mode = "head"
                
            res=re_prop_begin.search(line)
            if(res):
                print "prop begin"                
                istate_s = int(res.group(1))
                print "istate(s) = ", istate_s
                mode = "prop_begin"

            res=re_e_dmom.search(line)
            if(mode=="prop_begin" and res):
                print "prop_e_dmom begin"                
                for n in range(6):
                    line = f.readline()
                cistates[istate_s-1].e_dmom = map(float, line.split())
                mode = "none"

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
