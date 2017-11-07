import numpy as np
from csvio import *
import os

# ==== function ====
def asgn(c0,c,e):
    s = np.dot(np.transpose(c), c0)
    decouple(e,s)
    orthnorm(s)
    """
    n = len(s)
    asgn_q = False
    for i in range(n):
        if(abs(1-s[i,i]) > 10.0**(-10.0)):
            asgn_q = True
    """
    (c1,e1) = transf(s,c,e)
    return (c1, e1, s)

def is_used(j,jls):
    for i in range(len(jls)):
        if(jls[i] == 0):
            return False
        if(jls[i] == j):
            return True
    return False

def decouple(e,s):
    (n,m) = s.shape
    jls = np.zeros(n, dtype=int)
    jmax = 0
    for i in range(n):
        vmax = 0.0
        for j in range(n):
            if(is_used(j,jls)):
                next
            if(s[i,j]**2>vmax):
                vmax = s[i,j]**2
                jmax = j
        jls[i] = jmax
        for j in range(n):
            if(j!=jmax):
                s[i,j] = 0.0
                            
def orthnorm(s):
    tmp = np.dot(s[:,0], s[:,0])
    s[:,0] = s[:,0] / np.sqrt(tmp)
    for i in range(1,len(s[:,0])):
        for j in range(0,i):
            tmp = np.dot(s[:,i], s[:,j])
            s[:,i] = s[:,i] - tmp * s[:,j]
        tmp = np.dot(s[:,i], s[:,i])
        s[:,i] = s[:,i] / np.sqrt(tmp)
            
def transf(s,c,e):
    nb = len(s[:,0])
    ns = len(s[0,:])
    c1 = np.dot(c, s)
    return (c1, None)

if __name__ == '__main__':
    c = np.array([[0.98, 0.1],
                  [0.1,  0.92],
                  [0.0,  0.2]])
    c0= np.array([[0.2,  -0.98],
                  [0.92, 0.1],
                  [0.2,  0.0]])
    print c
    (c1,e1) = asgn(c0,c,None)
    print c1



# ==== IO ====
def ci_asgn(dirout0, dirout1):
    print "CI vector phase assignment"
    c0 = ijv2mat(pd.read_csv(dirout0 + "/cci.csv"))
    c1 = ijv2mat(pd.read_csv(dirout1 + "/cci.csv"))

    (c1p,ep,s) = asgn(c0,c1,None)

    print "save old ci vector as {0}/cci.old.csv".format(dirout1)
    os.system("mv {0}/cci.csv {0}/cci.old.csv".format(dirout1))
        
    print "save assigned ci vector as {0}/cci.csv".format(dirout1)
    df = mat2ijv(c1p)
    df.to_csv(dirout1+"/cci.csv", index=False)
    
    print "save asign transformation matrix as {0}/asgnci.csv".format(dirout1)
    df = mat2ijv(s, 0.001)
    df.to_csv(dirout1+"/asgnci.csv",index=False)

def mo_asgn(dirout0, dirout1):
    print "MO vector phase assignment"
    c0 = ijv2mat(pd.read_csv(dirout0 + "/cmo.csv"))    
    c1 = ijv2mat(pd.read_csv(dirout1 + "/cmo.csv"))
    smat0 = ijv2mat(pd.read_csv(dirout0 + "/s.csv"))

    sc0 = np.dot(smat0, c0)
    (c1p,ep,s) = asgn(sc0,c1,None)

    print "save old MO vector as {0}/cmo.old.csv".format(dirout1)
    os.system("mv {0}/cmo.csv {0}/cmo.old.csv".format(dirout1))
        
    print "save assigned ci vector as {0}/cmo.csv".format(dirout1)
    df = mat2ijv(c1p)
    df.to_csv(dirout1+"/cmo.csv", index=False)
    
    print "save asign transformation matrix as {0}/asgnmo.csv".format(dirout1)
    df = mat2ijv(s, 0.001)
    df.to_csv(dirout1+"/asgnmo.csv",index=False)

    
    





    
