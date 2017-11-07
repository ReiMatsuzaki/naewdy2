import numpy as np
import pandas as pd
import os
from naewdy2.nsh import *
from naewdy2.math import *

class MatNearError(Exception):
    def __init__(self, a, b, i, j, label):
        self.i = i
        self.j = j
        self.msg = """
        MatNearError
        difference of two matrix is too big.
        label = {0}
        i = {1}
        j = {2}
        a[i,j] = {3}
        b[i,j] = {4}
        """.format(label,i,j,a[i,j],b[i,j])

class MatShapeError(Exception):
    def __init__(self, a, b, label):
        self.msg = """
        MatShapeError
        shape of two matrixes are different
        label = {0}
        a.shape = {1}
        b.shape = {2}
        """.format(label, a.shape, b.shape)

def check_near_mat(a, b, label, eps):
    if(a.shape!=b.shape):
        raise MatShapeError(a,b,label)
    (n,m) = a.shape
    for i in range(n):
        for j in range(m):
            diff = abs(a[i,j]-b[i,j])
            if(diff>eps):
                raise MatNearError(a,b,i,j,label)

def calc_x(dir_dict, katom, dx, dir_out, m_xci="1", m_xmo="1"):
    """
    dir_dict : output directory list
    .          dir_dict[n] : calculation results for x0+n.dx
    katom : index of moving atom
    dx : float
    dir_out : result directory name
    """
    d0 = dir_dict[0]    
    tmat = ijv2mat(pd.read_csv(d0+"/t.csv"))
    smat = ijv2mat(pd.read_csv(d0+"/s.csv"))
    df_aij = pd.read_csv(d0+"/aij.csv")

    if(not os.path.exists(dir_out)):
        os.makedirs(dir_out)
    
    with open(d0+"/nshel.json") as f:
        j_nshel = json.load(f)
    gs = nshel2gto(j_nshel, normalize=True)

    with open(d0+"/common.json") as f:
        j = json.load(f)
    nao = j["nao"]
    nmo = j["nmo"]
    nwks = j["nwks"]
    nstate=j["nstate"]

    tr = np.transpose
    dot= np.dot

    # == check ==
    n = len(smat)
    eps = 10**(-7)
    
    smat_py = gtomat(gs, op_s())
    tmat_py = gtomat(gs, op_t())
    try:
        check_near_mat(smat_py, smat, "smat(py vs gms)", eps)
        check_near_mat(tmat_py, tmat, "tmat(py vs gms)", eps)
    except MatNearError as e:
        print e.msg
        print "gs[i]:"
        print gs[e.i]
        print "gs[j]"
        print gs[e.j]
    except MatShapeError as e:
        print e.msg
        raise RuntimeError("Error")

    # == X_MO  ==
    # - MO -
    cmo0 = ijv2mat(pd.read_csv(dir_dict[0 ]+"/cmo.csv"))
    csc0 = dot(tr(cmo0), dot(smat, cmo0))
    csc0_ref = np.identity(nmo)
    try:
        check_near_mat(csc0_ref, csc0, "cSc in MO.", 10**(-8))
    except MatNearError as e:
        print "Warning"
        print e.msg
        
    dzmat = gtomat(gs, op_dw(2))    
    for mu in range(nao):
        for nu in range(nao):
            if(gs[nu].katom != katom ):
                dzmat[mu,nu] = 0
    if(m_xmo=="analytic"):
        xmo = 0.5*dot(tr(cmo0), dot(dzmat-tr(dzmat), cmo0))
    elif(m_xmo=="1"):
        p1 = ijv2mat(pd.read_csv(dir_dict[1 ]+"/cmo.csv"))
        m1 = ijv2mat(pd.read_csv(dir_dict[-1]+"/cmo.csv"))
        xmo =  (dot(tr(cmo0), dot(smat, p1-m1))/(2*dx)
                + dot(tr(cmo0), dot(dzmat, cmo0)))
    elif(m_xmo=="1c"):
        p1 = ijv2mat(pd.read_csv(dir_dict[1 ]+"/cmo.csv"))
        m1 = ijv2mat(pd.read_csv(dir_dict[-1]+"/cmo.csv"))
        xmo =  dot(tr(cmo0), dot(smat, p1-m1))/(2*dx)
    elif(m_xmo=="2"):
        p1 = ijv2mat(pd.read_csv(dir_dict[1 ]+"/cmo.csv"))
        p2 = ijv2mat(pd.read_csv(dir_dict[2 ]+"/cmo.csv"))
        m1 = ijv2mat(pd.read_csv(dir_dict[-1]+"/cmo.csv"))
        m2 = ijv2mat(pd.read_csv(dir_dict[-2]+"/cmo.csv"))
        xmo =  (dot(tr(cmo0), dot(smat, -p2+8*p1-8*m1+m2))/(12*dx)
                +  dot(tr(cmo0), dot(dzmat, cmo0)))

    # - MO -> SCF -
    xcsf=np.zeros((nwks,nwks))
    num = len(df_aij)
    for idx in range(num):
        i = df_aij["i"][idx]
        j = df_aij["j"][idx]
        ii = df_aij["I"][idx]
        jj = df_aij["J"][idx]
        v = df_aij["val"][idx]
        if(ii!=jj):
            xcsf[ii-1,jj-1] += v*xmo[i-1,j-1]
            xcsf[jj-1,ii-1] += v*xmo[j-1,i-1]

    with open(dir_out+"/xcsf.csv", "w") as f:
        f.write("i,j,val\n")
        for i in range(nwks):
            for j in range(nwks):
                f.write("{0},{1},{2}\n".format(i+1,j+1,xcsf[i,j]))
        
    with open(dir_out+"/xmo.csv", "w") as f:
        f.write("i,j,val\n")
        n = len(xmo[0,:])
        for i in range(n):
            for j in range(n):
                f.write("{0},{1},{2}\n".format(i+1,j+1,xmo[i,j]))
            
    if(m_xci==None):
        return
    
    # - SCF -> CI -    
    cci0   = ijv2mat(pd.read_csv(d0+"/cci.csv"))
    xci_mo = dot(tr(cci0), dot(xcsf, cci0))

    # == X_CI (CI coef) ==
    if(m_xci=="1"):
        ccip = ijv2mat(pd.read_csv(dir_dict[1]+"/cci.csv"))
        ccim = ijv2mat(pd.read_csv(dir_dict[-1]+"/cci.csv"))
        xci_ci = (dot(tr(cci0), ccip) -dot(tr(cci0),ccim))/(2*dx)
                  
    elif(m_xci=="2"):
        [p1, p2, m1, m2] = [ ijv2mat(pd.read_csv(dir_dict[i]+"/cci.csv"))
                             for i in [1,2,-1,-2] ]
        xci_ci = ( -1*dot(tr(cci0),p2)
                   +8*dot(tr(cci0),p1)
                   -8*dot(tr(cci0),m1)
                   +1*dot(tr(cci0),m2) )/(12*dx)
    elif(m_xci=="3"):
        [p1, p2, p3, m1, m2, m3] = [ ijv2mat(pd.read_csv(dir_dict[i]+"/cci.csv"))
                                     for i in [1,2,3,-1,-2,-3] ]
        xci_ci = dot(tr(cci0),45*(p1-m1) -9*(p2-m2) + (p3-m3))/(60*dx)

        
    with open(dir_out+"/xci.csv","w") as f:
        f.write("i,j,xmo,xci\n")
        n = len(xci_mo[0,:])
        for i in range(n):
            for j in range(n):
                f.write("{0},{1},{2},{3}\n".format(i+1,j+1,xci_mo[i,j],xci_ci[i,j]))

    

                
