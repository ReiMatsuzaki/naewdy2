from datetime import datetime
import os
join = os.path.join
exists = os.path.exists
import sys
import re
import subprocess
import json

import pandas as pd
import numpy as np
dot = np.dot
tr = np.transpose

from enint.nsh import Nshel, nshel_load
from enint.ciwfn import ao2mo

from naewdy2.math import ijklv2ten4, ijv2mat, mat2ijv, mat2csv, csv2mat

gmsout = "."
run_naewdy2 = os.path.expanduser("~/src/naewdy2/gamess/gamess_naewdy/rungms")
run_naewdy3 = os.path.expanduser("~/src/naewdy3/gms/gms/rungms")
run_orig = os.path.expanduser("~/src/naewdy/gamess/gamess_orig/rungms")

class Gamess(object):
    def __init__(self, dir0=".", gmsinp="gms.inp", gmslog="gms.log"):
        self.gmsinp = join(dir0, gmsinp)
        self.gmslog = join(dir0, gmslog)
        self.gmsout = join(dir0, "gmsout")


def run_gms(base, dirname, jobname="gms", naewdy_json=None, oldout="skip",
            label=""):

    print ""
    print "run_gms begin", datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    print label
    print "dirname = ", dirname
    print "jobname = ", jobname

    if(exists(dirname)):
        print "calculation is already done."
        if(oldout=="remove"):
            print "remove old results."
            os.system("rm -rf {0}".format(dirname))
        elif(oldout=="skip"):
            print "Pass this calculation."
            return 1
    
    if(not exists(dirname)):
        os.makedirs(dirname)

    filename = join(dirname, jobname+".inp")
    with open(filename, "w") as f:
        f.write(base)

    if(naewdy_json is not None):
        fn = join(dirname, "naewdy.in.json")        
        if(isinstance(naewdy_json, str)):
            with open(fn, "w") as f:
                f.write(naewdy_json)
        else:
            with open(fn, "w") as f:
                naewdy_json.dump(f)

    if("NAEWDY3" in base):
        run = run_naewdy3
    if("NAEWDY2" in base):
        run = run_naewdy2
    else:
        run = run_orig

    logfile = jobname + ".log"
    errfile = jobname + ".err"
    
    curr_dirname = os.getcwd()

    os.chdir(dirname)
    fo = open(logfile,"w")
    fe = open(errfile,"w")
    subprocess.check_call([run, jobname], stdout=fo, stderr=fe)
    fo.close()
    fe.close()
    fe = open(errfile,"r")
    lines = fe.read()
    fe.close()
    if("ddikick.x: Fatal error detected." in lines):
        raise RuntimeError("error on calculation.")
    os.chdir(curr_dirname)
    
    print "calculation success"
    return 0
    
def get_vecdat(dirname, jobname="gms"):

    rhf_key = "--- CLOSED SHELL ORBITALS"
    mcscf_key = "--- OPTIMIZED MCSCF MO-S"
    
    dat_lines = ""
    fn = join(dirname, gmsout, "uscr/{0}.dat".format(jobname))
    finddat = False    
    with open(fn) as f:
        line = f.readline()
        while line:
            if(rhf_key in line or mcscf_key in line):
                finddat = True
                line = f.readline()
                line = f.readline()
                line = f.readline()
#                line = f.readline()
                while line:
                    if("END" in line):
                        dat_lines += line
                        break
                    dat_lines += line
                    line = f.readline()
                    
            if(finddat):
                break
            
            line = f.readline()
    return dat_lines
def vecdat_norb(vecdat):
    for line in vecdat.split("\n")[1:]:
        if("$END" in line):
            break
        norb = line.split()[0]
        
    return int(norb)


def need_make(target, reference_list):

    if(not isinstance(reference_list, list)):
        return need_make(target, [reference_list])

    for ref in reference_list:
        if(not exists(ref)):
            raise RuntimeError("""reference file is not exist.
file = {0}
            """.format(ref))
        
    if(not exists(target)):
        return True

    for ref in reference_list:
        if(os.stat(target).st_mtime < os.stat(ref).st_mtime):
            return True
    return False

def calc_xci(dir_dict, katom, dx, dir_res, m1=None, m2=None,
             label="", force=False):
    """
    dir_dict : output directory list
    .          dir_dict[n] : calculation results for x0+n.dx
    katom : index of moving atom
    dx : float
    dir_out : result directory name
    """

    if(dir_res is None):
        dir_res = dir_dict[0]
    
    print ""
    print "calc_xci begin", datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    print label
    print "dir_dict[0] = ", dir_dict[0]
    print "dir_res = ", dir_res
    print "katom:", katom

    fn_xci = join(dir_res, "xci.csv")
    if(not force and
       not need_make(fn_xci, join(dir_dict[0], "dm1.csv"))):
        print "pass calculation"
        return

    # -- basic information --
    d0 = dir_dict[0]
    with open(join(d0,"common.json")) as f:
        j = json.load(f)
    nao = j["nao"]
    nmo = j["nmo"]
    nwks = j["nwks"]
    nstate = j["nstate"]
    noccmo = j["nfrozen"] + j["nact"]

    with open(join(d0, "nshel.json")) as f:
        j = json.load(f)
    nsh = nshel_load(j)
    nsh.setup()
    
    d0 = dir_dict[0]
    dm1 = ijklv2ten4(join(d0, "dm1.csv"))
    cci0 = ijv2mat(join(d0, "cci.csv"))
    cmo0 = ijv2mat(join(d0, "cmo.csv"))[:,0:noccmo]
    smat = nsh.smat()
    dzmat = nsh.dwmat(2)
    ia_vec = nsh.ia_vec()

    # -- AO --
    print "AO begin"
    if(katom is None):
        dzmat[:,:] = 0.0
    elif(katom == "CNM"):
        if(m1 is None or m2 is None):
            raise RuntimeError("m1 and m2 is necessary for katom==CNM")
        for mu in range(nao):
            for nu in range(nao):
                if(ia_vec[nu] == 1):
                    dzmat[mu,nu] *= m2/(m1+m2)
                if(ia_vec[nu] == 2):
                    dzmat[mu,nu] *= -m1/(m1+m2)
    else:
        for mu in range(nao):
            for nu in range(nao):
                if(ia_vec[nu] != katom):
                    dzmat[mu,nu] = 0.0

    # -- MO --
    print "MO begin"
    dp = dir_dict[1]
    ccip = ijv2mat(join(dp, "cci.csv"))
    cmop = ijv2mat(join(dp, "cmo.csv"))[:,0:noccmo]

    dm = dir_dict[-1]
    ccim = ijv2mat(join(dm, "cci.csv"))
    cmom = ijv2mat(join(dm, "cmo.csv"))[:,0:noccmo]
    
    xmo =  (dot(tr(cmo0), dot(smat, cmop-cmom))/(2*dx)
            + dot(tr(cmo0), dot(dzmat, cmo0)))

    print "CI begin"
    xci_mo = []
    for n in range(nstate):
        for m in range(nstate):
            xci_mo.append(np.sum(dm1[:,:,n,m]*xmo))
    xci_mo = np.reshape(xci_mo, (nstate, nstate))
    xci_ci = (dot(tr(cci0), ccip) -dot(tr(cci0),ccim))/(2*dx)

    mat2csv(xci_mo, join(dir_res, "xmo.csv"))
    mat2csv(xci_ci, join(dir_res, "xci.csv"))

    # -- CSF --
    print "CSF begin"
    aij = pd.read_csv(join(d0, "aij.csv"))
    xcsf=np.zeros((nwks,nwks))
    num = len(aij)
    for idx in range(num):
        i = aij["i"][idx]
        j = aij["j"][idx]
        ii= aij["I"][idx]
        jj =aij["J"][idx]
        v = aij["val"][idx]
        if(ii!=jj):
            xcsf[ii-1,jj-1] += v*xmo[i-1,j-1]
            xcsf[jj-1,ii-1] += v*xmo[j-1,i-1]

    with open(join(dir_res, "xcsf.csv"), "w") as f:
        f.write("i,j,val\n")
        for i in range(nwks):
            for j in range(nwks):
                f.write("{0},{1},{2}\n".format(i+1,j+1,xcsf[i,j]))    

    print "calc_xci end"                

def calc_no(dir_out, dir_res=None, force=False):

    print ""
    print "calc_no begin"
    if(dir_res is None):
        dir_res = dir_out
        
    fn_res_val = join(dir_res, "no_val.csv")
    fn_res_vec = join(dir_res, "no_vec.csv")
    if(not force and
       not need_make(fn_res_vec, join(dir_out, "dm1.csv"))):
        print "pass calculation"
        return    
    
    with open(join(dir_out,"common.json")) as f:
        j = json.load(f)
    nstate = j["nstate"]
    noccmo = j["nfrozen"] + j["nact"]
    

    f_val = open(fn_res_val, "w")
    f_vec = open(fn_res_vec, "w")
    f_val.write("i,j,val\n")
    f_vec.write("i,j,k,val\n")
    dm1 = ijklv2ten4(join(dir_out, "dm1.csv"))
    for n in range(nstate):
        (ni, nvec) = np.linalg.eigh(dm1[:,:,n,n])
        for i in range(noccmo):
            f_val.write("{0},{1},{2}\n".format(n,i,ni[i]))
            for j in range(noccmo):
                f_vec.write("{0},{1},{2},{3}\n".format(n,i,j,nvec[i,j]))

    print "calc_no end"
                

    
def calc_dip(dir_out, dir_res=None, force=False, label="",
             fn_dm1="dm1.csv", fn_res="dip.csv"):
    
    if(dir_res is None):
        dir_res = dir_out
    
    print ""
    print "calc_dip begin", datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    print label
    print "dir_out = ", dir_out
    print "dir_res = ", dir_res

    fn_res = join(dir_res, fn_res)
    if(not force and
       not need_make(fn_res, join(dir_out, fn_dm1))):
        print "pass calculation"
        return

    with open(join(dir_out,"common.json")) as f:
        j = json.load(f)
    nao = j["nao"]
    nmo = j["nmo"]
    nwks = j["nwks"]
    nstate = j["nstate"]
    noccmo = j["nfrozen"] + j["nact"]

    with open(join(dir_out, "nshel.json")) as f:
        j = json.load(f)
    nsh = nshel_load(j)
    nsh.setup()

    dm1 = ijklv2ten4(join(dir_out, fn_dm1))
    cmo = ijv2mat(join(dir_out, "cmo.csv"))[:,0:noccmo]
    xmat = nsh.rmat(0)
    ymat = nsh.rmat(1)
    zmat = nsh.rmat(2)

    mat2ijv(xmat).to_csv(join(dir_out, "xao.csv"), index=None)
    mat2ijv(ymat).to_csv(join(dir_out, "yao.csv"), index=None)
    mat2ijv(zmat).to_csv(join(dir_out, "zao.csv"), index=None)

    xmat = ao2mo(xmat, cmo)
    ymat = ao2mo(ymat, cmo)
    zmat = ao2mo(zmat, cmo)

    mat = []
    for n in range(nstate):
        x = np.sum(xmat*dm1[:,:,n,n])
        y = np.sum(ymat*dm1[:,:,n,n])
        z = np.sum(zmat*dm1[:,:,n,n])
        mat.append([x,y,z])
    mat = np.array(mat)
    df = mat2ijv(mat)
    df.to_csv(fn_res, index=None)

    print "calc_dip end"


