import numpy as np
import json
from numpy import exp, sin, cos, pi, sqrt, dot
tr = np.transpose

class DVR:
    def __init__(self, num):
        self.num = num
        self.xs = np.zeros(num)
        self.ws = np.zeros(num)
        self.u  = np.zeros((num, num))

    def __str__(self):
        return "<DVR>"

    def phi(self, xs, n=0):
        """
        gives n derivative values of basis function set at x
        xs : [float]
        n : integer 
        """
        raise RuntimeError("abstract subroutine")
        
    def dmat_phi(self, n):
        """
        gives matrix of derivative operator in basis function
        n : integer [1,2]
        """

    def dmat(self, n):
        """
        gives matrix of derivative operator in DVR basis
        """
        

    def at(self, c, x_or_xs, nd=0):
        """
        gives function value with coefficient c.

        Input
        ------
        c : [scalar]
        x_or_xs : scalar or [scalar]
        nd : integer

        Result
        ------
        y_or_ys : 'same for x_or_xs'
        """
        if(len(c)!=self.num):
            raise RuntimeError("""
            len(c)!=self.num
            len(c)={0}
            self.num={1}""".format(len(c), self.num))
        
        if(isinstance(x_or_xs, list)):
            xs = np.array(x_or_xs)
        elif(isinstance(x_or_xs, np.ndarray)):
            xs = x_or_xs
        else:
            xs = np.array([x_or_xs])
            pass

        phiss = self.phi(xs, nd)
        uc = dot(self.u, c)
        ys = np.dot(uc, phiss)
#        for x in xs:
#            phis = self.phi(x, nd)
#            y = dot(uc, phis)
#            ys.append(y)
#            
        if(len(ys)==1):
            return ys[0]
        else:
            return np.array(ys)

    def fit(self, fun):
        return np.array([fun(x)*sqrt(w) for (x,w) in zip(self.xs, self.ws)])

    def to_json(self, fn):
        pass
        
class SinDVR(DVR):
    def __init__(self, num, x0, xNp1):
        DVR.__init__(self, num)
        self.x0 = x0
        self.xNp1 = xNp1
        self.L = (xNp1-x0)
        self.dx = self.L/(num+1)
        
        for al in range(num):
            self.xs[al] = x0 + al*self.dx
            self.ws[al] = self.dx

        for j in range(num):
            for al in range(num):
                self.u[j,al] = sqrt(2.0/(num+1)) * sin(j*al*pi/(num+1))
                
    def phi(self, x, n=0):
        
        ys = np.zeros(self.num)
        if(x<self.x0 or self.xNp1<x):
            return ys
        
        z = pi/self.L
        zz = z*(x-self.x0)
        a = sqrt(2/self.L)
        if(n==0):
            for j in range(self.num):
                ys[j] = a*sin(zz*j)
        elif(n==1):
            for j in range(self.num):
                ys[j] = a*j*cos(zz*j)
        elif(n==2):
            for j in range(self.num):
                ys[j] = a*j*j*cos(zz*j)
        else:
            raise RuntimeError("n<-[0,1,2]")
        return ys

    def dmat_phi(self, n):
        m = np.zeros((self.num, self.num))
        if(n==1):
            for j in range(self.num):
                for k in range(self.num):
                    if((j-k)%2==0):
                        m[j,k] = 0.0
                    else:
                        m[j,k] = 4/self.L * (1.0*j*k)/(j*j-k*k)
        elif(n==2):
            for j in range(self.num):
                m[j,j] = -(j*pi/self.L)**2
        else:
            raise RuntimeError("n<-[1,2]")
        return m

    def dmat(self, n):
        if(n!=2):
            raise RuntimeError("n<-[2]")

        n = self.num
        m = np.zeros((n,n))
        for a in range(n):
            for b in range(n):
                sa = sin(a*pi/(n+1))
                ca = cos(a*pi/(n+1))
                if(a==b):
                    m[a,a] = -1.0/3.0 + 1.0/(6*(n+1)**2) - 1.0/(2*(n+1)**2 *sa**2)
                else:
                    sb = sin(b*pi/(n+1))
                    cb = cos(b*pi/(n+1))
                    m[a,b] = (2*(-1)**(a-b))/((n+1)**2) * (sa*sb)/((ca-cb)**2)
        m = -(pi/self.dx)**2 * m
        return m

    def to_json(self, fn):
        raise RuntimeError("not impl")
    
class ExpDVR(DVR):
    def __init__(self, n, x0, xN):
        DVR.__init__(self, 2*n+1)
        self.u  = np.zeros((self.num, self.num), dtype=complex)
        self.n = n
        self.x0 = x0
        self.xN = xN
        self.L = xN-x0
        self.dx = self.L/self.num
        for al in range(self.num):
            self.xs[al] = x0 + (al+1)*self.dx
            self.ws[al] = self.dx
        for j in range(-n,n+1):
            for al in range(self.num):
                self.u[n+j,al] = sqrt(1.0/self.num) * exp(-2.0j*j*(al+1)*pi/(self.num))

    def __str__(self):
        return "ExpDVR(n={0}, x0={1}, x1={2})".format(self.n, self.x0, self.xN)

    def phi(self, xs, nd=0):
        yss = np.zeros((self.num, len(xs)), dtype=complex)
        r = 2.0j*pi/self.L
        s = sqrt(1/self.L)
        for j in range(-self.n, self.n+1):
            z = r*j
            a = s * (z**nd) * exp(-z*self.x0)
            yss[j+self.n,:] = a*exp(z*xs)
            #for ix in range(len(xs)):                
            #    yss[j+self.n,ix] = a*exp(z*xs[ix])
                
            #z = 2.0j*pi*j/self.L
            #ys[j+self.n] = sqrt(1/self.L) * (z**nd) * exp(z*(x-self.x0))
        return yss
    
    def dmat(self, nd):
        m = np.zeros((self.num, self.num))
        if(nd==1):
            for a in range(self.num):
                for b in range(self.num):
                    if(a==b):
                        m[a,a] = 0.0
                    else:
                        m[a,b] = pi/self.L * ((-1)**(a-b)) / sin(pi*(a-b)/self.num)
        elif(nd==2):
            for a in range(self.num):
                for b in range(self.num):
                    if(a==b):
                        m[a,a] = -pi*pi/(3*self.L**2) * (self.num**2-1)
                    else:
                        s = sin(pi*(a-b)/self.num)
                        c = cos(pi*(a-b)/self.num)
                        m[a,b] = -2*pi*pi/(self.L**2)*(-1)**(a-b)*c/(s*s)
        else:
            raise RuntimeError("nd<-[1,2]")
        return m


    
    def to_json(self, fn):
        dic0 = {}
        dic0["type"] = "exp"
        dic0["n"] = self.n
        dic0["x0"] = self.x0
        dic0["xN"] = self.xN
        json_str = json.dumps(dic0)
        with open(fn, "w") as f:
            f.write(json_str)

def dvr_load(fn):
    with open(fn, "r") as f:
        str_json = f.readlines()
    dict_json = json.loads(str_json)
    if(dict_json["type"]=="exp"):
        dvr = ExpDVR(dict_json["n"],
                     dict_json["x0"],
                     dict_json["xN"] )
        return dvr
    else:
        raise RuntimeError("not supported")
            
