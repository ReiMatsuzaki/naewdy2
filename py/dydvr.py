import numpy as np
import os
import datetime

from naewdy2.dvr import *
from naewdy2.math import *

dot = np.dot
tr = np.transpose
join = os.path.join
exists = os.path.exists

def write_csvrow(f, xs):
    nx = len(xs)
    ixs = range(nx)
    for ix in ixs:
        f.write(str(xs[ix]))
        if(ix==nx-1):
            f.write("\n")
        else:
            f.write(",")

class DyDvrVal(object):
    def __init__(self, fn=None):
        self.fn = fn
        
    def init(self, name, dtype, out="out"):
        self.name = name
        self.dtype = dtype
        self.out = out
        
    def calc(self, c,**kwargs):
        pass

class DyDvrValOne(DyDvrVal):
    def init(self, name, dtype, out="out"):
        DyDvrVal.init(self, name, dtype, out)
        
        self.fp = open(os.path.join(out, name+".csv"), "w")
        if(type(dtype)!=list):
            raise RuntimeError("dtype must be type of list")
        headers = []
        nd = len(dtype)
        for (d,i) in zip(dtype, range(nd)):
            if(nd==1):
                pre = ""
            else:
                pre = str(i+1)
            if(d==float):                
                headers.append("val"+pre)
            elif(d==complex):
                headers.append("re"+pre)
                headers.append("im"+pre)
            else:
                raise RuntimeError("invalid type")
        write_csvrow(self.fp, headers)

    def calc(self,c,**kwargs):
        """
        c : matrix[nnuc,nel]
        """
        ys = self._calc(c,**kwargs)
        zs = []
        for (y,dt) in zip(ys,self.dtype):
            if(dt==float):
                zs.append(y)
            elif(dt==complex):
                zs.append(y.real)
                zs.append(y.imag)
            else:
                raise RuntimeError("invalid dtype")
                
        write_csvrow(self.fp, zs)

    def _calc(self,c,**kwargs):
        pass
        
class Norm(DyDvrValOne):
    def init(self, **kwargs):
        DyDvrValOne.init(self, "norm", [float])

    def _calc(self,c,**kwargs):
        return [(np.sum(c.conj()*c)).real]
        
class Ac(DyDvrValOne):
    def init(self, **kwargs):
        DyDvrValOne.init(self, "ac", [complex])
        self.c0 = kwargs["c0"]

    def _calc(self,c,**kwargs):        
        return [(np.sum(c.conj()*self.c0))]

class R(DyDvrValOne):    
    def init(self, **kwargs):
        """
        <R> = sum_{ai} <c_ai g_ai | R | c_ai g_ai>
        .   = sum_{ai} |c_ai|^2 R_a
        """
        DyDvrValOne.init(self, "r", [float])

    def _calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        rs = np.array([dvr.xs]).transpose()
        rc = rs*c
        y = np.sum(c.conj()*rc)
        return [y.real]


class P(DyDvrValOne):
    def init(self, **kwargs):
        """
        <P> = -i <psi| D |psi>
        .   = -i sum_{abI} c_{aI}^* D_{ab} c_{bI}
        """        
        DyDvrValOne.init(self, "p", [float])

    def _calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        d = dvr.dmat(1)
        cd = np.dot(d, c)
        y = -1.0j * np.sum(c.conj()*cd)
        return [y.real]

class Prob(DyDvrValOne):
    def init(self, **kwargs):
        self.nel = kwargs["nel"]
        DyDvrValOne.init(self, "prob", [float]*self.nel)

    def _calc(self,c,**kwargs):
        ys = []
        for i in range(self.nel):
            y = np.dot(c[:,i].conj(), c[:,i])
            ys.append(y.real)
        return ys
    
class Psi(DyDvrVal):
    def __init__(self, xs, fn="psi.csv"):
        self.xs = xs
        self.fn = fn
        
    def init(self, **kwargs):
        self.nel = kwargs["nel"]
        if(self.nel > 20):
            raise RuntimeError("nel exceed 20. You should not use this as outpu")

    def calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        dir0 = kwargs["dir0"]        
        ilist = range(self.nel)
        columns = ["x"]
        dat = { "x":self.xs }
        for i in ilist:
            lre = "re"+str(i+1)
            lim = "im"+str(i+1)
            columns.append(lre)
            columns.append(lim)
            y = dvr.at(c[:,i], self.xs)
            dat[lre] = y.real
            dat[lim] = y.imag

        df = pd.DataFrame(dat, columns=columns)
        df.to_csv(join(dir0, self.fn), index=None)
        
class ProbCI(DyDvrValOne):
    def __init__(self, nstate):
        self.nstate = nstate
        
    def init(self, **kwargs):
        """
        P_n = |<Psi_n|Psi>|^2
        .   = sum_a C_a_in^* C_ia
        """
        DyDvrValOne.init(self, "probci", [float for n in range(self.nstate)])        

    def _calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        out_qc = kwargs["out_qc"]

        ys = np.zeros(self.nstate)
        nnuc = dvr.num
        for a  in range(nnuc):
            fn = os.path.join(out_qc, str(a), "0/out/cci.csv")
            if(not exists(fn)):
                continue
            
            cci = ijv2mat(pd.read_csv(fn))

            for n in range(self.nstate):
                y = np.dot(cci[:,n].conj(), c[a,:])
                ys[n] += (abs(y))**2
                
        return ys

class PsiCI(DyDvrVal):
    def init(self,**kwargs):
        DyDvrVal.init("psici", [complex])

    def _calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        out_qc = kwargs["out_qc"]
        out_dy = kwargs["out_dy"]
        nnuc = dvr.num

        fs = None
        for  a in range(nnuc):
            fn = os.path.join(out_qc, str(a), "0/out/cci.csv")
            if(not exists(fn)):
                continue
            
            cci = ijv2mat(pd.read_csv(fn))

            if(len(fs)==None):
                nstate = len(cci[0,:])
                for n in range(nstate):                    
                    fs.append(open(os.path.join(out_dy, "")))

class Rho(DyDvrVal):
    def __init__(self, xs):
        """
        sum_I |<Phi_I|Psi(t)>|^2 = sum_I |c_aI.u_a(x)|^2
        .                        =
        """        
        self.xs = xs

    def init(self,**kwargs):
        DyDvrVal.init(self, "rho", [float])

    def calc(self,c,**kwargs):
        dvr = kwargs["dvr"]
        dir0 = kwargs["dir0"]
        nel = len(c[0,:])

        nx = len(self.xs)
        
        rho = np.zeros(nx)
        for i in range(nel):
            psii = dvr.at(c[:,i],self.xs)
            rho[:] += (psii.conj()*psii).real
        
        with open(join(dir0, "rho.csv"), "w") as f:
            f.write("x,val\n")
            for (x,y) in zip(self.xs, rho):
                f.write("{0},{1}\n".format(x, y))
        
def analyze(vals, dvr, nel, out_dy="out/dy", out="out", out_qc="out/qc"):
    nnuc = dvr.num
    df = pd.read_csv(os.path.join(out_dy, "0/c.csv"))
    c0 = ijv2mat(df)

    for v in vals:
        v.init(nel=nel,c0=c0)
    
    dirs = os.listdir(out_dy)
    dirs.sort(key=int)
    for dir0 in dirs:
        print "calculation for ", dir0
        df = pd.read_csv(os.path.join(out_dy, dir0, "c.csv"))
        c1 = ijv2mat(df)
        for v in vals:
            v.calc(c1,dvr=dvr,out_qc=out_qc,out_dy=out_dy,dir0=join(out_dy, dir0))
    
def analyze_old(cmd_list, nel, dvr, out="out"):
    """
    read calculation result of DyDVR, calculate physical values.
    Assume that common.py is located 
    
    Inputs
    ------
    cmd_list : [String or dict]
    nel: integer
    dvr: DVR object
    out : str
    .    calculation results
    """
    
    cmd_dict = {}
    for cmd in cmd_list:
        if(isinstance(cmd, str)):
            cmd_dict[cmd] = {}
        elif(isinstance(cmd, dict)):
            cmd_dict[cmd["type"]] = cmd

    if(not os.path.exists(out)):
        print "Error. {0} not found".format(out)
        
    nnuc = dvr.num
    df = pd.read_csv(out + "/dy/0/c.csv")
    c0 = np.reshape(iv2vec(df), (nnuc, nel))

    for cmd in cmd_dict:
        if(cmd in ["norm", "ac", "prob"]):
            cmd_dict[cmd]["file_pointer"] = open("{0}/{1}.csv".format(out,cmd), "w")
            f = cmd_dict[cmd]["file_pointer"]
        if(cmd=="norm"):
            f.write("val\n")
        elif(cmd=="ac"):
            f.write("re,im\n")
        elif(cmd=="prob"):
            for i in range(nel):
                if(i==nel-1):
                    f.write("val{0}\n".format(i+1))
                else:
                    f.write("val{0},".format(i+1))
    
    dirs = os.listdir(out+"/dy")
    dirs.sort(key=int)
    for dir0 in dirs:
        print "calculation for ", dir0
        df = pd.read_csv("{0}/dy/{1}/c.csv".format(out, dir0))
        c1 = np.reshape(iv2vec(df), (nnuc, nel))
        if("norm" in cmd_dict):
            norm2 = np.sum(c1.conj()*c1)
            f = cmd_dict["norm"]["file_pointer"]
            f.write("{0}\n".format(norm2.real))

        if("ac" in cmd_dict):
            ac = np.sum(c1.conj()*c0)
            f = cmd_dict["ac"]["file_pointer"]
            f.write("{0},{1}\n".format(ac.real,ac.imag))

        if("prob" in cmd_dict):
            for i in range(nel):
                f = cmd_dict["prob"]["file_pointer"]
                y = np.dot(c1[:,i].conj(), c1[:,i]).real
                if(i==nel-1):
                    f.write("{0}\n".format(y))
                else:
                    f.write("{0},".format(y))
                    
        if("psi" in cmd_dict):
            dict_psi = cmd_dict["psi"]
            xs = dict_psi["xs"]
            for i in range(nel):
                ys = dvr.at(c1[:,i], xs, nd=0)
                df = pd.DataFrame({"x":xs, "re":ys.real, "im":ys.imag},
                                  columns=["x","re","im"])
                df.to_csv("{0}/dy/{1}/psi{2}.csv".format(out,dir0,i+1),index=None)

def dvrc2molc(nucc, nel, iel):
    """
    gives nuclear-electronic coefficient from nuclear coefficient and electron state index
    """
    nnuc = len(nucc)
    c = np.zeros(nnuc*nel, dtype=complex)
    for a in range(nnuc):
        c[idx(nel,a,iel-1)] = nucc[a]
    return c
    
def hdot(x, y):
    return dot(x.conj(), y)

def norm(x):
    return np.sqrt(hdot(x,x).real)

def get_now():    
    return datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    
    
class DyDVR:
    def __init__(self, mass, dvr, c0, nel, out="out",
                 inte="diag", krylov_num=None):
        print "DyDVR.init begin"

        self.inte = inte
        print "inte:", self.inte
        if(inte not in ["diag", "krylov", "eig"]):
            raise RuntimeError("inte<-[diag, krylov]")

        self.krylov_num = krylov_num
        if(inte == "krylov" and krylov_num == None):
            raise RuntimeError("using inte==krylov, krylov_num is necessary ")
        
        self.mass = mass
        self.dvr = dvr

        self.D1 = self.dvr.dmat(1)
        self.D2 = self.dvr.dmat(2)

        if((dvr.num, nel) != c0.shape):
            raise RuntimeError("""(dvr.num,nel)!=c0.shape
dvr.num = {0}
nel = {1}            
c0.shape = {2}
""".format(dvr.num, nel, c0.shape))
            

        self.c = np.reshape(c0, (nel*dvr.num))
        norm2 = dot(self.c.conj(), self.c)
        self.c = self.c/np.sqrt(norm2)

        self.out = out
        self.nel = nel

        if(not os.path.exists(out)):
            os.makedirs(out)
        
        self.ft = open("{0}/t.csv".format(out),"w")
        self.ft.write("val\n")

        print "DyDVR.init end"

    def precalc_data(self, hel):
        """
        hel : hel[a,i,j] gives electronic Hamiltonian matrix element (i,j) at dvr.xs[a]
        """
        self.nel = hel.shape[1]
        D1 = self.dvr.dmat(1)
        D2 = self.dvr.dmat(2)

        nnuc = self.dvr.num
        nel = self.nel
        h = np.zeros((nel*nnuc, nel*nnuc))
        for a in range(nnuc):
            for i in range(nel):
                for b in range(nnuc):
                    for j in range(nel):
                        ai = self.idx(a,i)
                        bj = self.idx(b,j)
                        if(a==b):
                            tmp = hel[a,i,j]
                            h[ai,bj] = h[ai,bj] + tmp
                        if(i==j):
                            h[ai,bj] = h[ai,bj] -1/(2*self.mass)*D2[a,b]
                            
        (self.e, self.u) = eigh_sort(h)
        self.uH = tr(self.u)
        #(self.e, self.u) = np.linalg.eigh(h)
        #(self.e, self.u) = np.linalg.eig(h)
        #df = vec2iv(self.e)
        #df.to_csv(self.out+"/e.csv",index=None)
#
#        df = mat2ijv(h)
#        df.to_csv(self.out+"/h.csv",index=None)
#        
#        self.uH = np.conjugate(tr(self.u))

        
    def precalc_data_nac1_slow(self, hel, xij):
        """
        hel : hel[a,i,j] gives electronic Hamiltonian matrix element (i,j) at dvr.xs[a]
        xij : xij[a,i,j] gives derivative coupling matrix element (i,j) at dvr.xs[a]
        """
        print "DyDVR.precalc_data_nac1_slow begin"
        self.nel = hel.shape[1]
        D1 = self.dvr.dmat(1)
        D2 = self.dvr.dmat(2)

        nnuc = self.dvr.num
        nel = self.nel
        h = np.zeros((nel*nnuc, nel*nnuc))
        for a in range(nnuc):
            for i in range(nel):
                for b in range(nnuc):
                    for j in range(nel):
                        ai = self.idx(a,i)
                        bj = self.idx(b,j)
                        if(a==b):
                            tmp = hel[a,i,j]
                            h[ai,bj] = h[ai,bj] + tmp
                        if(i==j):
                            h[ai,bj] = h[ai,bj] -1/(2*self.mass)*D2[a,b]
                        h[ai,bj] = h[ai,bj] -1/self.mass * xij[a,i,j]*D1[a,b]
                            
        (self.e, self.u) = eigh_sort(h)
        self.uH = tr(self.u)
        print "DyDVR.precalc_data_nac1_slow end"

    def precalc_data_nac1(self, hel, xij):
        """
        hel : hel[a,i,j] gives el Hamiltonian matrix element (i,j) at dvr.xs[a]
        xij : xij[a,i,j] gives derivative coupling matrix element (i,j) at dvr.xs[a]
        """
        print "DyDVR.precalc_data_nac1 begin", get_now()
        print "inte:", self.inte
        self.nel = hel.shape[1]

        if(self.inte == "krylov"):
            self.hel = hel
            self.xij = xij

        if(self.inte == "eig" or self.inte == "diag"):
            print "DyDVR.precalc_data_nac1 build H mat begin", get_now()
            self.h = build_H(hel, xij, self.D1, self.D2, self.mass)
            print "DyDVR.precalc_data_nac1 build H mat end", get_now()

        if(self.inte=="eig"):
            print "DyDVR.precalc_data_nac1 diag H begin", get_now()
            (self.e, self.u) = np.linalg.eig(h)
            if(dtype==complex):
                self.uH = tr(self.u.conjugate())
            elif(dtype==float):
                self.uH = tr(self.u)
            else:
                raise RuntimeError("invalid dtype. dtype={0}".format(dtype))
            print "DyDVR.precalc_data_nac1 diag H end"
        
        print "DyDVR.precalc_data_nac1 end", get_now()
        
    def update(self, it, dt, ntskip):

        t = it*dt*ntskip
        self.ft.write("{0}\n".format(t))

        dn = "{0}/dy/{1}".format(self.out, it)
        if(os.path.exists(dn)):
            print "directory exists. pass calculation and read coefficient"
            cmat = ijv2mat(pd.read_csv(dn+"/c.csv"))
            self.c = np.reshape(cmat, (self.dvr.num*self.nel))
            return
        else:
            print "make directory and calculation begin"
            print "dir:", dn
            os.makedirs(dn)

            print "dump current data"
            cmat = np.reshape(self.c, (self.dvr.num,self.nel))
            df = mat2ijv(cmat)
            df.to_csv(dn+"/c.csv",index=None)

            for j in range(ntskip):
                if(self.inte=="eig"):
                    opts = {"inte":"eig"}
                    self.c = uni_inte((self.e,self.u,self.uH), dt, self.c,
                                      opts)
                elif(self.inte=="diag" ):
                    opts = {"inte":"diag"}
                    self.c = uni_inte(self.h, dt, self.c, opts)
                elif(self.inte=="krylov"):
                    Hc = sigma_H(self.hel, self.xij, self.D1, self.D2, self.mass)
                    self.c = uni_inte_krylov(
                        Hc, dt, self.krylov_num, self.c)
                else:
                    raise RuntimeError("invalid self.inte")
        
    # == Utils ==
    def idx(self,a,i):
        return idx(self.nel,a,i)


    
