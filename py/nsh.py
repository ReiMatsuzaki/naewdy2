import json
import numpy as np
from numpy import sqrt, exp, pi
from scipy.special import gamma, erf, gammainc

class GTO:
    def __init__(self, ex, c, R, n, katom=0, do_normalize=False):
        self.ncont = len(ex)
        self.ex = np.array(ex)
        self.c  = np.array(c)
        self.w = np.array(R)
        self.n = np.array(n)
        self.katom = 0
        if(do_normalize):
            smat = gtomat([self], op_s())
            self.c = self.c / np.sqrt(smat[0,0])
            
    def __str__(self):
        return """
        ==== GTO ====
        ex: {0}
        c:  {1}
        w:  {2}
        n:  {3}
        """.format(self.ex, self.c, self.w, self.n)
                
def gtoint(n, z):
    if(n%2==1):
        return 0
    if(n==0):
        return np.sqrt(np.pi/z)
    return (n-1)/(2*z) * gtoint(n-2,z)

def mole_gammainc(m, z):
    """
    Compute incomplete gamma function defined below
    .    F_m(z) = Int_0^1 t^{2m} Exp[-zt^2]dt
    Impled with scipy.special.gammainc(a,x) defined below
    .    1/Gamma[a] Int_0^x t^{a-1} Exp[-t] dt
    Variable transformation 
    .    t->(t^2)/x
    .    dt -> (2t/x) dt
    .    gammainc(a,x) = x^{1-a} /Gamma[a]  Int_0^1 t^{2a-2} (2t/x) Exp[-(t^2)/x] dt
    .                  = 2 (x^(-a))/gamma[a] Int_0^1 t^{2a-1} Exp[-(t^2)/x] dt
    .   F_m(z) = (1/2) gamma[a] x^a 

    Special case
    .  z = 0
    .  F_m(0) = Int_0^1 t^{2m} dt 
    .         = (2m+1)^{-1} [t^{2m+1}]_0^1 
    .         = (2m+1)^{-1}
    """

    eps = 10.0**(-14.0)
    
    if(abs(z)<eps):
        return 1/(2*m+1.0)
    
    if(m==0):
        a = m+0.5
        res = gamma(a)/(2*z**a)*gammainc(a, z)
#        res = sqrt(pi/z)/2 * erf(sqrt(z))
    else:
        a = m+0.5
        res = gamma(a)/(2*z**a)*gammainc(a, z)
        
    if((not res<1) and (not res>-1)):
        raise RuntimeError("""gammain failed
        res: {2}
        m: {0}
        z: {1}
        """.format(m,z,res))
    return res
    
def coef_d(zp,wpk,wak,wbk,nak,nbk,nk):
    if(nak==0 and nbk ==0 and nk ==0):
        return 1.0
    if(nk<0 or nk>nak+nbk):
        return 0.0
    if(nak>0):
        return (1/(2*zp) * coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk-1) +
                (wpk-wak)* coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk)   +
                (nk+1.0) * coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk+1)   )
    else:
        return (1/(2*zp) * coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk-1) +
                (wpk-wbk)* coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk)   +
                (nk+1.0) * coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk+1)   )

def coef_R(zp,wp,c,m,j):
    """
    Compute function R_{m,j}(zp,w,c). See T.Kuchitsu, J.Okuda and M.Tachikawa, Int.J.Q.Chem. 109, 540 (2008)

    Inputs
    -------
    zp: float
    wp: [float]*3
    c : [float]*3
    m : [int]*3
    j : int

    Returns
    --------
    R_{m,j}(zp,w,c) : float
    """

    wc = wp-c

    if(m[0]==0 and m[1]==0 and m[2]==0):
        try:
            return (-2*zp)**j * mole_gammainc(j,zp*dist2(wc))
        except RuntimeError as e:
            raise RuntimeError("""
failed at mole_gammainc. 
zp = {0}
wp = {1}
c  = {2}
m  = {3}
j  = {4}
Below are errror message from mole_grammainc
{5}
            """.format(zp,wp,c,m,j,e.message))

    for i in range(3):
        im = np.zeros(3); im[i] = 1
        if(m[i]>0):
            res = wc[i] * coef_R(zp,wp,c,m-im,j+1)
            if(m[i]>1):
                res += (m[i]-1) * coef_R(zp,wp,c,m-2*im,j+1)
            return res

    raise RuntimeError("one of m is negative integer")
        

def coef_d_list(zp,wp,wa,wb,mna,mnb):
    """
    ds = np.zeros((3,ma+1,mb+1,ma+mb+2))
    for i in range(3):
        for na in range(ma+1):
            for nb in range(mb+1):
                for n in range(ma+mb+2):
                    ds[i,na,nb,n] = coef_d(zp,wp[ir],wa[ir],wb[ir],nj,nk,n)
    """
    ds = [[[[coef_d(zp,wp[ir],wa[ir],wb[ir],nj,nk,n)
             for n in range(mna+mnb+2)]
            for nk in range(mnb+1)]
           for nj in range(mna+1)]
          for ir in range(3)]
    ds = np.array(ds)
    return ds

def coef_R_list(zp,wp,wc,maxn,n):
    rs = [[[coef_R(zp, wp, wc, [nx,ny,nz], n)
             for nz in range(maxn+1)]
            for ny in range(maxn+1)]
           for nx in range(maxn+1)]
    rs = np.array(rs)
    return rs
    
def dist2(xs):
    acc = 0
    for x in xs:
        acc += x*x
    return acc
    

def prim_s(na, wa, za, nb, wb, zb):
    zp = za+zb
    wp = (za*wa + zb*wb)/zp
    d2 = sum([x*x for x in wa-wb])
    ep = np.exp(-za*zb/zp*d2)
    c = 1    
    for i in range(3):
        cd = coef_d(zp,wp[i],wa[i],wb[i],na[i],nb[i],0)
        c = c * cd
    return ep * (np.sqrt(np.pi/zp))**3 * c

def prim_r(i):
    def __func__(na, wa, za, nb, wb, zb):
        naa = [n for n in na]
        naa[i] = naa[i] + 1
        return (prim_s(naa,wa,za, nb,wb,zb)
            + wa[i] * prim_s(na,wa,za, nb,wb,zb))
    return __func__

def prim_dw(i):
    def __func__(na, wa, za, nb, wb, zb):
        #   Dw.(x-w)^n Exp[-z(x-w)^2]
        # = {-n(x-w)^{n-1} +2z(x-w)^{n+1}} Exp[-z(x-w)^2]
        nbp = [n for n in nb]
        nbp[i] = nbp[i]+1
        acc = 2*zb*prim_s(na,wa,za, nbp,wb,zb)
        if(nb[i]!=0):
            nbm = [n for n in nb]
            nbm[i] = nbm[i]-1
            acc = acc - nb[i]*prim_s(na,wa,za, nbm,wb,zb)
        return acc
    return __func__

def prim_na(wc):
    def __func__(na, wa, za, nb, wb, zb):
        zp = za+zb
        wp = (za*wa + zb*wb)/zp
        d = wp-wc
        d2p = dist2(d)
        d2  = dist2(wa-wb)
        ep  = np.exp(-za*zb*d2/zp)

        res = 0
        ns = np.zeros(3, dtype=int)
        for nx in range(na[0]+nb[0]+1):
            cx = coef_d(zp,wp[0],wa[0],wb[0],na[0],nb[0],nx)
            for ny in range(na[1]+nb[1]+1):
                cy = coef_d(zp,wp[1],wa[1],wb[1],na[1],nb[1],ny)
                for nz in range(na[2]+nb[2]+1):
                    ns[0]=nx; ns[1]=ny; ns[2]=nz; 
                    cz = coef_d(zp,wp[2],wa[2],wb[2],na[2],nb[2],nz)
                    cr = coef_R(zp,wp,wc,ns,0)
                    res += cx*cy*cz*cr
        return -2*np.pi*ep*res/zp
    return __func__

def prim_t(na, wa, za, nb, wb, zb):    
    
    s000 = prim_s(na, wa, za, nb, wb, zb)

    sp = [0,0,0]
    for i in range(3):
        npb = [n for n in nb]
        npb[i] = nb[i]+2
        sp[i] = prim_s(na, wa, za, npb, wb, zb)

    res = -2*zb*(2*nb[0] + 2*nb[1] + 2*nb[2] + 3) * s000
    res+=  4*zb*zb*(sp[0] + sp[1] + sp[2])

    for i in range(3):
        if(nb[i]>1):
            npb = [n for n in nb]
            npb[i] = nb[i]-2
            s = prim_s(na, wa, za, npb, wb, zb)
            res += nb[i]*(nb[i]-1) * s
    return -0.5*res

class OneOp():
    def __init__(self,prim,op_type="H"):
        self.op_type = op_type
        self.prim = prim

def op_s():
    return OneOp(prim_s)

def op_t():
    return OneOp(prim_t)

def op_r(i):
    return OneOp(prim_r)

def op_dw(i):
    return OneOp(prim_dw(i), "A")

def op_na(wc):
    return OneOp(prim_na(wc))

def gtoele(ga, op, gb):
    acc = 0.0
    for i in range(ga.ncont):
        for j in range(gb.ncont):
            acc += (ga.c[i]*gb.c[j]*
                    op.prim(ga.n, ga.w, ga.ex[i],
                            gb.n, gb.w, gb.ex[j]))
    return acc
    
def gtomat(gs, op):
    num = len(gs)
    mat = np.zeros((num,num))
    for a in range(num):
        for b in range(num):
            mat[a,b] = gtoele(gs[a], op, gs[b])
    return mat
            
def nshel2gto(j, normalize=False):
    ex = np.array(j["ex"])
    nshell = j["nshell"]
    ng = j["ng"]
    ian = j["ian"]
    c = np.array(j["c"])
    gs = []
    l_kmin_dict = {1:0, 2:1, 5:2, 11:3} # see inputa.src
    l_kmax_dict = {1:0, 4:1,10:2, 20:3} # see inputa.src
    for k in range(nshell):
        kstart = j["kstart"][k]
        katom =  j["katom"][k]
        ktype =  j["ktype"][k]
        kng   =  j["kng"][k]
        kex = [ex[kk-1] for kk in range(kstart,kstart+kng)]
        kks = range(kstart,kstart+kng)
        ls = range(l_kmin_dict[j["kmin"][k]],
                   l_kmax_dict[j["kmax"][k]]+1)
        kc = c[:,katom-1]
        if(0 in ls):
            kcs = [j["cs"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcs, kc, [0,0,0], katom, normalize))
        if(1 in ls):
            kcp = [j["cp"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcp, kc, [1,0,0], katom, normalize) )
            gs.append(GTO(kex, kcp, kc, [0,1,0], katom, normalize) )
            gs.append(GTO(kex, kcp, kc, [0,0,1], katom, normalize) )
        if(2 in ls):
            kcd = [j["cd"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcd, kc, [2,0,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,2,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,0,2], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [1,1,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [1,0,1], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,1,1], katom, normalize))
                      
    return gs


class Nshel_v1:
    def __init__(self):
        self.nshel = 0
        self.ng = 0
        
        self.ex = []        
        
        self.c  = [[], [], []]
        self.ian = []
        self.zan = []

        self.kstart = []
        self.katom  = []
        self.ktype  = []
        self.kng    = []
        self.kmin   = []
        self.kmax   = []
        self.kloc   = []

        self.ntable = np.array([
            [0,0,0],
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
            [2,0,0],
            [0,2,0],
            [0,0,2],
            [1,1,0],
            [1,0,1],
            [0,1,1]])
        
        self.ltable = np.sum(self.ntable, axis=1)

        num_n = len(self.ntable[:,0])
        self.coef = [[] for n in range(num_n)]
            
    def load(self, j):
        self.ex = np.array(j["ex"])
        self.nshel = j["nshell"]
        self.ng = j["ng"]
        self.ian = j["ian"]
        self.c  = np.array(j["c"])

        self.kstart = j["kstart"]
        self.katom  = j["katom"]
        self.ktype  = j["ktype"]
        self.kng    = j["kng"]
        self.kmin   = j["kmin"]
        self.kmax   = j["kmax"]
        self.kloc   = j["kloc"]

        cs = j["cs"]
        cp = j["cp"]
        cd = j["cd"]

        for k in range(len(self.ntable)):
            l = self.ltable[k]
            if(l==0):
                self.coef[k] = cs
            elif(l==1):
                self.coef[k] = cp
            elif(l==2):
                self.coef[k] = cd

    def add_atom(self, ci, anum, z):
        ia = len(self.c[0])
        for ir in range(3):
            self.c[ir].append(ci[ir])
        self.ian.append(anum)
        self.zan.append(z)
        return ia
    
    def add_shel(self, ltype, cs, ex, ia):
        self.ex += ex
        self.nshel += 1
        self.ng += len(ex)
            
        self.katom.append(ia)        
        self.ktype.append(0)
        self.kng.append(len(ex))

        
        if(isinstance(ltype, str)):
            if(ltype=="s"):
                self.kmin.append(1)
                self.kmax.append(1)
                l = 0                
            elif(ltype=="p"):
                self.kmin.append(2)
                self.kmax.append(4)
                l = 1
            elif(ltype=="d"):
                self.kmin.append(5)
                self.kmax.append(10)
                l = 2
            else:
                raise RuntimeError("not impl")
        elif(isinstance(ltype, list)):
            ins = self.get_ins(ltype)
            self.kmin.append(ins)
            self.kmax.append(ins)
            l = np.sum(ltype)
        else:
            raise RuntimeError("not impl")

        inlist = range(self.kmin[-1], self.kmax[-1]+1)
        num_n = len(self.ntable[:,0])
        c0 = [0 for n in range(len(cs))]
        for i in range(num_n):
            if(i in inlist):
                self.coef[i] += cs
            else:
                self.coef[i] += c0[:]

    def setup(self):

        self.kstart.append(1)
        self.kloc.append(1)
        for k in range(1,self.nshel):
            self.kstart.append(self.kstart[k-1] + self.kng[k-1])
            nloc = self.kmax[k-1] - self.kmin[k-1] + 1
            self.kloc.append(self.kloc[k-1] + nloc)

        self.c = np.array(self.c)

        self.coef = np.array(self.coef)
        coef_old = np.copy(self.coef)
        
        s = self.smat()
        idx = -1
        for ishel in range(self.nshel):
            ips = self.get_iprims(ishel)
            for k in range(self.kmin[ishel], self.kmax[ishel]+1):
                idx += 1
                for ip in ips:
                    self.coef[k,ip] = coef_old[k,ip]/np.sqrt(s[idx,idx])
                            
    def smat(self):
        mat = np.zeros((21,21))

        for ishel in range(self.nshel):
            ips = self.get_iprims(ishel)
            ns_i = self.ntable[self.kmin[ishel]:self.kmax[ishel]+1,:]
            max_ni = sum(ns_i[-1])
            iatom = self.katom[ishel]-1
            ci = self.c[:,iatom]
            
            
            for jshel in range(self.nshel):
                jps = self.get_iprims(jshel)
                ns_j = self.ntable[self.kmin[jshel]:self.kmax[jshel]+1,:]
                max_nj = sum(ns_j[-1])
                jatom = self.katom[jshel]-1
                cj = self.c[:,jatom]
                locj = self.kloc[jshel]-1
                
                d2 = sum([x*x for x in ci-cj])
                
                for ip in ips:
                    for jp in jps:
                        zi = self.ex[ip]
                        zj = self.ex[jp]                        
                        zp = zi+zj
                        wp = (zi*ci+zj*cj)/zp
                        ep = np.exp(-zi*zj/zp*d2)
                        ds = [ [ [ coef_d(zp,wp[ir],ci[ir],cj[ir],ni,nj,0)
                                   for nj in range(max_nj+1)]
                                 for ni in range(max_ni+1)]
                               for ir in range(3)]

                        i = self.kloc[ishel]-1-1
                        for ik in range(self.kmin[ishel],
                                       self.kmax[ishel]+1):
                            i += 1
                            j = self.kloc[jshel]-1-1
                            for jk in range(self.kmin[jshel],
                                       self.kmax[jshel]+1):
                                j += 1
                                acc = 1.0                                
                                for ir in range(3):
                                    ni = self.ntable[ik,ir]
                                    nj = self.ntable[jk,ir]
                                    acc *= ds[ir][ni][nj]
                                acc *= ep*(np.pi/zp)**(1.5)
                                li = self.ntable[ik]
                                lj = self.ntable[jk]
                                coef = self.coef[ik,ip]*self.coef[jk,jp]
                                mat[i,j]+=acc*coef
        return mat
                
    def get_iprims(self, ishel):
        start_i = self.kstart[ishel]
        ng_i = self.kng[ishel]
        ips = range(start_i-1,start_i+ng_i-1)
        return ips
        
    def get_ns(self, ishel):
        l_kmin_dict = {1:0, 2:1, 5:2, 11:3} # see inputa.src
        l_kmax_dict = {1:0, 4:1,10:2, 20:3} # see inputa.src
        ls = range(l_kmin_dict[self.kmin[ishel]],
                   l_kmax_dict[self.kmax[ishel]]+1)
        ns = []
        if(0 in ls):
            ns.append([0,0,0])
        if(1 in ls):
            ns.append([1,0,0])
            ns.append([0,1,0])
            ns.append([0,0,1])
        if(2 in ls):
            ns.append([2,0,0])
            ns.append([0,2,0])
            ns.append([0,0,2])
            ns.append([1,1,0])
            ns.append([1,0,1])
            ns.append([0,1,1])
        return (ns, np.max(ls))

    def get_ins(self, n):
        for ins in range(1, len(self.ntable[:,0])):
            if(n[0]==self.ntable[ins,0] and
               n[1]==self.ntable[ins,1] and
               n[2]==self.ntable[ins,2]):
                return ins

def prod(xs):
    acc = 1
    for x in xs:
       acc *= x
    return acc

class Shel:
    def __init__(self, ntypes, ex, coef_l, w):

        if(isinstance(ntypes, str)):
            ntypes = [ntypes]

        self.ns = []
        if("s" in ntypes):
            self.ns.append([0,0,0])
        if("px" in ntypes):
            self.ns.append([1,0,0])            
        if("p" in ntypes):
            self.ns.append([1,0,0])
            self.ns.append([0,1,0])
            self.ns.append([0,0,1])
        if("dxx" in ntypes):
            self.ns.append([2,0,0])
        if("d" in ntypes):
            self.ns.append([2,0,0])
            self.ns.append([0,2,0])
            self.ns.append([0,0,2])
            self.ns.append([1,1,0])
            self.ns.append([1,0,1])
            self.ns.append([0,1,1])
        if("f" in ntypes):
            self.ns.append([3,0,0])
            self.ns.append([0,3,0])
            self.ns.append([0,0,3])
            self.ns.append([2,1,0])
            self.ns.append([2,0,1])
            self.ns.append([1,2,0])
            self.ns.append([0,2,1])
            self.ns.append([1,0,2])
            self.ns.append([0,1,2])
            self.ns.append([1,1,1])

        self.num = len(self.ns)

        self.ns = np.array(self.ns)
        self.max_n = max([sum(n) for n in self.ns])
        
        self.ex = ex
        self.ng = len(ex)

        self.coef = np.zeros((self.num, self.ng))
        for jn in range(self.num):
            nj = self.ns[jn]
            l = sum(nj)
            self.coef[jn,:] = coef_l[l]
        
        if((self.num, self.ng) != self.coef.shape):
            raise RuntimeError("size mismatch")

        self.w = np.array(w)

        self.j0 = 0
    
class Nucs:
    def __init__(self):
        self.ws = []
        self.anum = []
        self.zs = []
        self.num = 0
    def add_atom(self,w,anum,z):
        self.ws.append(w)
        self.anum.append(anum)
        self.zs.append(z)
        self.num += 1
        return self.num-1


class Nshel:
    def __init__(self, nucs):
        self.shels = []
        self.nucs = nucs

    def __str__(self):
        line = "jshel  w   n  ex  coef   \n"
        for jshel in range(len(self.shels)):
            shel = self.shels[jshel]
            
            line += " {0}  {1}\n".format(jshel, shel.w)

            for jn in range(shel.num):
                line += "    {0}\n".format(shel.ns[jn])
                for ig in range(shel.ng):
                    line += "     {0} {1} {2}\n".format("A", shel.ex[ig],
                                                        shel.coef[jn,ig])
        return line
            
    def add_shel(self,ntypes,ex,coef,ia):
        self.shels.append(Shel(ntypes,ex,coef,self.nucs.ws[ia]))
        
    def setup(self, normalize=True):
        jn = 0
        for shel in self.shels:
            shel.j0 = jn
            jn += shel.num

        if(normalize):
            s = self.smat()
            idx = -1
            for shel in self.shels:
                coef_old = np.copy(shel.coef)
                for jj in range(shel.num):
                    idx += 1
                    shel.coef[jj,:] = coef_old[jj,:]/np.sqrt(s[idx,idx])
    
    def smat(self):
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])                

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = ep*(np.pi/zp)**(1.5)

                        ds = [ [ [ coef_d(zp,wp[ir],wj[ir],wk[ir],nj,nk,0)
                                   for nk in range(sk.max_n+1)]
                                 for nj in range(sj.max_n+1)]
                               for ir in range(3)]

                        for jj in range(sj.num):
                            for kk in range(sk.num):
                                nj = sj.ns[jj]
                                nk = sk.ns[kk]
                                acc = prod([ds[ir][nj[ir]][nk[ir]]
                                            for ir in range(3)])
                                coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                j = sj.j0 + jj
                                k = sk.j0 + kk
                                mat[j,k] += acc*coef
        return mat
    
    def tmat(self):
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = ep*(np.pi/zp)**(1.5)

                        ds = [ [ [ coef_d(zp,wp[ir],wj[ir],wk[ir],nj,nk,0)
                                   for nk in range(sk.max_n+2+1)]
                                 for nj in range(sj.max_n+1)]
                               for ir in range(3)]

                        for jj in range(sj.num):
                            for kk in range(sk.num):
                                nj = sj.ns[jj]
                                nk = sk.ns[kk]
                                i3 = range(3)

                                acc = 0
                                s000 = prod([ds[i][nj[i]][nk[i]] for i in i3])
                                acc += -2*zk*(2*nk[0]+2*nk[1]+2*nk[2]+3)*s000
                                
                                for jr in i3:
                                    nkp = np.copy(nk[:])
                                    nkp[jr] += 2
                                    s = prod([ds[i][nj[i]][nkp[i]] for i in i3])
                                    acc += 4*zk*zk*s
                                    if(nk[jr]>1):
                                        nkm = np.copy(nk[:])
                                        nkm[jr] -= 2
                                        s = prod([ds[i][nj[i]][nkm[i]] for i in i3])
                                        acc += nk[jr]*(nk[jr]-1)*s
                                        
                                coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                j = sj.j0 + jj
                                k = sk.j0 + kk
                                mat[j,k] += acc*coef
        mat *= -0.5
        return mat        
        
    def vmat(self, nucs=None):
        if(nucs==None):
            nucs = self.nucs
            
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = -2*np.pi*ep/zp
                        
                        ds = coef_d_list(zp,wp,wj,wk,sj.max_n,sk.max_n)
                        
                        for ic in range(nucs.num):
                            nmax = sj.max_n+sk.max_n
                            rs = coef_R_list(zp,wp,nucs.ws[ic],nmax,0)
                              
                            for jj in range(sj.num):
                                for kk in range(sk.num):
                                    nj = sj.ns[jj]
                                    nk = sk.ns[kk]
                                    i3 = range(3)
                                    
                                    acc = 0
                                    for nx in range(nj[0]+nk[0]+1):
                                        for ny in range(nj[1]+nk[1]+1):
                                            for nz in range(nj[2]+nk[2]+1):
                                                acc += (ds[0,nj[0],nk[0],nx]*
                                                        ds[1,nj[1],nk[1],ny]*
                                                        ds[2,nj[2],nk[2],nz]*
                                                        rs[nx,ny,nz]*
                                                        nucs.zs[ic])
                                    coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                    j = sj.j0 + jj
                                    k = sk.j0 + kk
                                    mat[j,k] += acc*coef
        return mat                
            
    def to_gtos(self):
        gtos = []
        for shel in self.shels:            
            for jn in range(shel.num):
                gtos.append(GTO(shel.ex,
                                shel.coef[jn,:],
                                shel.w,
                                shel.ns[jn,:]))
        return gtos
        
def nshel_load(j):
    
    zan = j["ian"]
    ian = j["ian"]
    w  = np.array(j["c"])
    natom = len(zan)
    nucs = Nucs()
    for ia in range(natom):
        nucs.add_atom(w[:,ia], ian[ia], zan[ia])

    nshel = Nshel(nucs)
    
    ex = np.array(j["ex"])
    num_shel = j["nshell"]
    ng = j["ng"]
    kstart = j["kstart"]
    katom  = j["katom"]
    ktype  = j["ktype"]
    kng    = j["kng"]
    kmin   = j["kmin"]
    kmax   = j["kmax"]
    kloc   = j["kloc"]

    cs = j["cs"]
    cp = j["cp"]
    cd = j["cd"]
    cf = j["cf"]

    l_kmin = {1:0, 2:1, 5:2, 11:3} # see inputa.src
    l_kmax = {1:0, 4:1,10:2, 20:3} # see inputa.src

    for ishel in range(num_shel):
        ls = range(l_kmin[kmin[ishel]], l_kmax[kmax[ishel]]+1)
        spd = ["s", "p", "d", "f"]
        ntypes = [spd[l] for l in ls]
        ig0 = kstart[ishel]-1
        ig1 = kstart[ishel] + kng[ishel]-1
        
        ex_i = ex[ig0:ig1]
        cs_i = cs[ig0:ig1]
        cp_i = cp[ig0:ig1]
        cd_i = cd[ig0:ig1]
        cf_i = cf[ig0:ig1]

        coef_l = [None,None,None,None]
        coef_l[0] = cs_i
        coef_l[1] = cp_i
        coef_l[2] = cd_i
        coef_l[3] = cf_i
        
        nshel.add_shel(ntypes, ex_i, coef_l, katom[ishel]-1)
        
    return nshel

        
