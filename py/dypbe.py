import csv
import json
import os

from naewdy2.dvr import DVR

join = os.path.join
exists = os.path.exists
tr = np.transpose

def idx(nel, a, i):
    return a*nel+i
    
class DyPBE(object):
    def __init__(self, dvr, ehren, c0, dt, nt
                 out_dy = "out/dy", inte="diag", krylov_num=None):
        print "DyPBE.init begin"

        self.dvr = dvr
        self.mass = ehren.mass

        if((dvr.num,nel) != c0.shape):
            raise RuntimeError("(dvr.num,nel)!=c0.shape")

        self.c = np.reshape(c0, nel*dvr.num)
        self.nel = len(ehren.c)
        self.dt = dt
        ehren.dt = dt
        self.nt = nt
        self.it = 0

        self.out_dy = out_dy
        self.out = out
        self.fp_pq = open(join(out, "pq.csv"), "w")

        if(not os.path.exists(out)):
            os.makedirs(out)
        
        self.inte = inte
        if(inte not in ["diag", "krylov"]):
            raise RuntimeError("inte <- [diag, krylov] ")

        self.tmat = -1/(2*mass) * dvr.dmat(2)
        
        print "DyPBE.init end"

    def update(self, hIJq, dhIJq, xIJq):
        """
        hIJq : matrix(float)
        .      Electronic Hamiltonian evaluated at q
        dhIJq : matrix(float)
        .      Derivative of Electronic Hamiltonian evaluated at q
        xIJq : matrix(float)
        .      Derivative coupling evaluated at q
        """

        print "DyPBE.update begin. ({0}/{1})".format(self.it, self.nt)

        dn_out = join(self.out_dy, str(self.it))
        fn_c   = join(dn_out, "pbe_c.csv")
        fn_ec  = join(dn_out, "ehren_c.csv")
        fn_dat = join(dn_out, "data.json")

        if(exists(fn_c)):
            print "calculation result file exists. pass calculation and read this"
            j = json.load(open(fn_dat, "w"))
            ehren.q = j["q"]
            ehren.p = j["p"]

            df = pd.read_csv(fn_ec)
            ehren.c = iv2vec(df)

            df = pd.read_csv(fn_c)
            cmat = ijv2mat(df)
            self.c = np.reshape(cmat, self.dvr.num*self.nuel)
            
        else:
            if(not exists(dn_out)):
                os.makedirs(dn_out)
            
            if(self.it !=0 ):
                if(self.inte == "diag"):
                    update_diag(hIJq, dhIJq, xIJq)
                elif(self.inte == "krylov"):
                    update_krylov(hIJq, dhIJq, xIJq)

            try:
                ehren.update(hIJq, dhIJq, xIJq)
            except:
                print "Error on ehren. Stop calculation ..."
                sys.exit()

            t = self.it*self.dt
            with open(fn_dat, "w") as f:
                json.dump({"t":t, "q":ehren.q, "p":ehren.p}, f)

            df = vec2iv(ehren.c)
            df.to_csv(fn_ec,index=None)
            
            cmat = np.reshape(self.c, (self.dvr.num, self.nel))
            df = mat2ijv(cmat)
            df.to_csv(fn_c,index=None)

        print "DyPBE.update end"

    def update_diag(self, hIJq, dhIJq, xIJq):
        
        print "DyPBE.update_diag begin"
        
        n = self.nel*self.dvr.num
        h = np.zeros((n,n))

        t = -1/(2*self.mass)*self.dvr.dmat(2)
        dq = self.p/self.mass

        elid = np.identity(nel)
        for a in range(nnuc):
            a0 = self.idx(a,0)
            a1 = self.idx(a+1,0)
            h[a0:a1,a0:a1] = -dq*xIJq + hIJq
            for b in range(nnuc):
                b0 = self.idx(b,0)
                b1 = self.idx(b+1,0)
                h[a0:a1,b0:b1] += t[a,b]

        (e,u) = np.linalg.eigh(h)
        uH = tr(u.conj())

        self.c = dot(uH, self.c)
        self.c = np.exp(-1.0j*e*self.dt) * self.c
        self.c = dot(u, self.c)
        
        print "DyPBE.update_diag end"
            
    def idx(self,a,i):
        return idx(self.nel,a,i)

