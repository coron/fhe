# 
# utils.sage
# 
# Copyright (c) 2012 Jean-Sebastien Coron <jean-sebastien.coron@uni.lu> 
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation.


from itertools import chain,islice,izip,izip_longest
from sys import stdout

def RandomOdd(n):
  return 2*(2^(n-2)+ZZ.random_element(2^(n-2)))+1

def sum2(li):
  empty=True
  for x in li:
    if empty:
      res=x
      empty=False
    else:
      res=res+x
  if empty: 
    return 0
  else:
    return res

class array(list):
  def __mul__(self,x):
    return array([x*ci for ci in self])

  def __add__(self,xx):
    return array([ai+xi for ai,xi in zip(self,xx)])

def reduceSteps(f,li,verbose=True):
  i=0
  res=li[0]
  if verbose: print i,
  for x in li[1:]:
    stdout.flush()
    i=i+1
    res=f(res,x)
    if verbose: print i,
  if verbose: print
  return res

def sumBinary(a,b):
  "Computes the sum of the binary vectors a and b, modulo 2^n where n is the length of the vectors a and b"
  c=[a[0]+b[0]]
  carry=a[0]*b[0]

  for i in range(1,len(a)-1):
    carry2=(a[i]+b[i])*carry+a[i]*b[i]    
    c.append(a[i]+b[i]+carry)             
    carry=carry2
  
  c.append(a[-1]+b[-1]+carry)
  return c

def QuotientNear(a,b):
  "Gives the nearest integer to a/b"
  return (2*a+b)//(2*b)

def modNear(a,b):
  "Computes a mod b with a \in ]-b/2,b/2]"
  return a-b*QuotientNear(a,b)

def toBinary(x,l):
  "Converts a positive integer x into binary with l digits"
  return (x+2^l).digits(2)[:-1]
  
def prodScal(x,y):
  return sum2((xi*yi for xi,yi in izip(x,y)))

def randSparse(n,h):
  v=[0 for i in range(n)]
  while sum(v)<h:
    i=int(ZZ.random_element(n))
    if v[i]==0: v[i]=1
  return v
