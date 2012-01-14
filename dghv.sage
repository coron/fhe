# 
# dghv.sage: an implementation of DGHV-FHE with compressed public-key
# 
# as described in:
#
# J.S. Coron, D. Naccache and M. Tibouchi, "Public-key Compression and Modulus 
# Switching for Fully Homomorphic Encryption over the Integers"
#
# available at http://eprint.iacr.org/2011/440
#
# Copyright (c) 2012 Jean-Sebastien Coron <jean-sebastien.coron@uni.lu> 
# and Mehdi Tibouchi <mehdi.tibouchi@normalesup.org>
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation.

load "scalprod.spyx"
load "utils.sage"
from itertools import chain,izip,islice,count
from time import time

theta=15
n=4

class Pk(object):
  "The PK of the DGHV somewhat homomorphic scheme, without the xi's"
  def __init__(self,rho,eta,gam,modx0=True,verbose=True,*args,**kwargs):
    print "  Pk: generation of p"
    self.rho,self.eta,self.gam,self.modx0=rho,eta,gam,modx0
    self.p=random_prime(2^self.eta,lbound=2^(self.eta-1),proof=False)
    if self.modx0:
      self.x0=self.p*RandomOdd(self.gam-self.eta)

  def encrypt(self,m):
    return self.p*ZZ.random_element(2^(self.gam-self.eta-1))+2*ZZ.random_element(-2^self.rho+1,2^self.rho)+m

  def noise(self,c):
    return modNear(c,self.p)

  def decrypt(self,c):
    return mod(self.noise(c),2)

  def add(self,c1,c2):
    return self.modx0 and mod(c1+c2,self.x0).lift() or (c1+c2)

  def mult(self,c1,c2):
    return self.modx0 and mod(c1*c2,self.x0).lift() or (c1*c2)

  def __repr__(self):
    return "<Pk with rho=%d, eta=%d, gam=%d>" % (self.rho,self.eta,self.gam) 

class Ciphertext():
  def __init__(self,val_,pk_,degree_=1):
    self.val,self.pk,self.degree=val_,pk_,degree_

  @staticmethod
  def encrypt(pk,m=None):
    if m==None: m=ZZ.random_element(2)
    return Ciphertext(pk.encrypt(m),pk)
  
  def noise(self):
    return self.pk.noise(self.val)

  def decrypt(self,verbose=False):
    t=cputime(subprocesses=True)
    m=self.pk.decrypt(self.val)
    if verbose: print "Decrypt",cputime(subprocesses=True)-t
    return m

  def __repr__(self):
    return "<Ciphertext with m=%d size=%d noise=%d deg=%d>" % (self.decrypt(),self.val.nbits(),self.noise().nbits(),self.degree)

  def __add__(self,x):
    return self.__class__(self.pk.add(self.val,x.val),self.pk,max(self.degree,x.degree))

  def __mul__(self,x):
    return self.__class__(self.pk.mult(self.val,x.val),self.pk,self.degree+x.degree)

  def scalmult(self,x):
    if isinstance(x,array):
      return array([self.__class__(self.val*xi,self.pk,self.degree) for xi in x])
    else:
      return self.__class__(self.val*x,self.pk,self.degree)

  def expand(self):
    return self.pk.expand(self.val)

class PRIntegers: 
  "A list of pseudo-random integers."
  def __init__(self,gam,ell):
    self.gam,self.ell=gam,ell
    self.li=[None for i in range(self.ell)]
    set_random_seed()
    self.se=initial_seed()

  def __getitem__(self,i):
    return self.li[i]
    
  def __setitem__(self,i,val):
    self.li[i]=val

  def __iter__(self):
    set_random_seed(self.se)
    for i in range(self.ell):
      a=ZZ.random_element(2^self.gam)
      if self.li[i]!=None:
        yield self.li[i]
      else:
        yield a
    set_random_seed()


class PRIntegersDelta(PRIntegers):
  """A list of pseudo-random integers, with their delta corrections"""
  def __iter__(self):
    return (c+d for c,d in izip(PRIntegers.__iter__(self),self.delta))

  def ciphertexts(self,pk):
    return (Ciphertext(cd,pk) for cd in self)

  @staticmethod
  def encrypt(pk,v):
    pr=PRIntegersDelta(pk.gam,len(v))
    r=[ZZ.random_element(-2^pk.rho+1,2^pk.rho) for i in range(len(v))]
    pr.delta=[0 for i in range(len(v))]
    
    # We have to do it in two passes, otherwise malloc error.
    # Namely if we write:
    #   temp=[-mod(xi,pk.p).lift() for xi in pr]
    # then for some reason the xi's are not garbage collected.

    temp=[-mod(xi,pk.p) for xi in PRIntegers.__iter__(pr)]
    pr.delta=[te.lift()+2*ri+vi for te,ri,vi in izip(temp,r,v)]
    return pr

class PkRecrypt(Pk):
  "The Pk of the DGHV scheme, with the xi's, the yi's and the encrypted secret-key bits"
  def __init__(self,rho,eta,gam,Theta,pkRecrypt=None,*args,**kwargs):
    t=cputime(subprocesses=True)
    super(PkRecrypt,self).__init__(rho,eta,gam,*args,**kwargs)
    self.kap=64*(self.gam//64+1)-1
    self.Theta=Theta
    self.alpha,self.tau=kwargs['alpha'],kwargs['tau']

    xp=QuotientNear(2^self.kap,self.p)
    B=self.Theta/theta

    assert self.Theta % theta==0,"Theta must be a multiple of theta"
    print "  PkRecrypt: generation of s"
    self.s=[1]+[0 for i in range(B-1)]+sum2([randSparse(B,1) for i in range(theta-1)])

    t2=cputime(subprocesses=True)
    self.y=PRIntegers(self.kap,self.Theta)
    self.y[0]=0
    self.y[0]=mod(xp-prodScal(self.s,self.y),2^(self.kap+1)).lift()
    assert mod(prodScal(self.s,self.y),2^(self.kap+1))==mod(xp,2^(self.kap+1)),"Equality not valid"
    print "  PkRecrypt: generation of yis",cputime(subprocesses=True)-t2

    t2=cputime(subprocesses=True)
    self.x=PRIntegersDelta.encrypt(self,[0 for i in range(self.tau)])
    print "  PkRecrypt: generation of xis",cputime(subprocesses=True)-t2
    
    t2=cputime(subprocesses=True)
    self.pkRecrypt=pkRecrypt
    if not self.pkRecrypt: self.pkRecrypt=self
    self.se=PRIntegersDelta.encrypt(self.pkRecrypt,self.s)
    print "  PkRecrypt: generation of seis",cputime(subprocesses=True)-t2

    print "  PkRecrypt: genkey",cputime(subprocesses=True)-t

  def encrypt_pk(self,m=None,verbose=True):
    t=cputime(subprocesses=True)
    if m==None: m=ZZ.random_element(2)
    rhop=self.eta-self.rho

    f=[Ciphertext(ZZ.random_element(2^self.alpha),self) for i in range(self.tau)]

    c=sum2(ci*fi for ci,fi in izip(self.x.ciphertexts(self),f))+\
           Ciphertext(ZZ.random_element(2^rhop),self) + Ciphertext(m,self)
    if verbose: print "Encrypt",cputime(subprocesses=True)-t
    c.degree=1
    return c

  def expand(self,c,verbose=True):
    m=2^(n+1)-1
    t=cputime(subprocesses=True)

    ce=[(((directProd(c,yi,self.kap) >> (32-(n+2)))+1) >> 1) & m for yi in self.y]

    if verbose: print "expand",cputime(subprocesses=True)-t
    return ce

  def recrypt(self,c,verbose=False):
    t=cputime(subprocesses=True)
    B=self.Theta/theta
    cebin=[array(toBinary(cei,n+1)) for cei in c.expand()]

    sec=(c.__class__(ci,self.pkRecrypt) for ci in self.se)

    li=(ski.scalmult(cei) for ski,cei in izip(sec,cebin))
    ly=[sum2(islice(li,B)) for i in range(theta)]

    res=reduceSteps(sumBinary,ly,verbose)
    v=res[-1]+res[-2]
    v.val+=c.val & 1
    if verbose: print "recrypt",cputime(subprocesses=True)-t
    return v

def testPkRecrypt():
  toy={'ty':"toy",'lam':42,'rho':26,'eta':988,'gam':147456,'Theta':150,'pksize':0.076519,'seclevel':42.0,'alpha':936,'tau':158}
  small={'ty':"small",'lam':52,'rho':41,'eta':1558,'gam':843033,'Theta':555,'pksize':0.437567,'seclevel':52.0,'alpha':1476,'tau':572}
  medium={'ty':"medium",'lam':62,'rho':56,'eta':2128,'gam':4251866,'Theta':2070,'pksize':2.207241,'seclevel':62.0,'alpha':2016,'tau':2110}
  large={'ty':"large",'lam':72,'rho':71,'eta':2698,'gam':19575950,'Theta':7965,'pksize':10.303797,'seclevel':72.0,'alpha':2556,'tau':7659}

  for param in [toy,small,medium,large]:
    ty,rho,eta,gam,Theta=[param[x] for x in ['ty','rho','eta','gam','Theta']]
    print "type=",ty,"lam=",param['lam'],"rho=",rho,"eta=",eta,"gamma=",gam,"Theta=",Theta
    pk=PkRecrypt(rho,eta,gam,Theta,tau=param['tau'],alpha=param['alpha'])

    c1=Ciphertext.encrypt(pk)
    c2=Ciphertext.encrypt(pk)
    print "c1=",c1
    print "c2=",c2
    print "c1+c2=",c1+c2
    t=cputime(subprocesses=True)
    print "c1*c2=",c1*c2
    print "tmult:",cputime(subprocesses=True)-t

    c=pk.encrypt_pk()
    c.decrypt(verbose=True)
    print "c=",c
    newc=pk.recrypt(c,verbose=True)
    print "newc=",newc
    print

