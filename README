	An implementation of the DGHV fully homomorphic scheme

This is an implementation of the DGHV fully homomorphic scheme with
compressed public-key. This implementation is described in the
following article: 

[1] J.S. Coron, D. Naccache and M. Tibouchi, "Public-key Compression and
Modulus Switching for Fully Homomorphic Encryption over the Integers",
Proceedings of Eurocrypt 2012.

available at http://eprint.iacr.org/2011/440

The implementation is done with the SAGE 4.7.2 mathematical library
under Python, available at http://www.sagemath.org/.

WHAT IS FULLY HOMOMORPHIC ENCRYPTION ?
--------------------------------------

An encryption scheme is fully homomorphic when it is possible to
perform implicit addition and multiplication of plaintexts while
manipulating only ciphertexts. The first construction of a fully
homomorphic encryption (FHE) scheme was described by Gentry in 2009.  

WHAT IS THE DGHV SCHEME ?
-------------------------

The DGHV scheme is a FHE scheme published in:

[2] M. van Dijk, C. Gentry, S. Halevi and V. Vaikuntanathan, "Fully
  Homomorphic Encryption over the Integers". Proceedings of Eurocrypt
  2010.

and available at http://eprint.iacr.org/2009/616

In DGHV a ciphertext has the form:

c= q*p + 2*r + m

where p is the secret-key, m is the bit plaintext (m=0 or 1), q is a
large random, and r is a small random.

To decrypt, compute m=(c mod p) mod 2. It is easy to see that
decryption works as long as the noise 2*r is smaller than p.

Given two ciphertexts 
c1= q1*p + 2*r1 + m1
c2= q2*p + 2*r2 + m2

we have:

c1+c2=(q1+q2)*p+2*(r1+r2)+m1+m2

Therefore one can obtain the encryption of m1+m2 (mod 2)=m1 xor m2 simply by
computing c1+c2.  

Similarly we have:

c1*c2=q12*p+2*(2*r1*r2+r1*m2+r2*m1)+m1*m2

Therefore one can obtain the encryption of m1*m2 simply by computing
c1*c2. However one gets a new ciphertext with noise roughly
twice larger than in the original ciphertexts c1 and c2. Since the
noise must remain below p, the number of permitted multiplications on
ciphertexts is therefore limited. This is called a somewhat
homomorphic encryption scheme. 
 
To obtain a fully homomorphic encryption scheme, i.e. unlimited
addition and multiplication on ciphertexts, one must be able to reduce
the amount of noise in a ciphertext; this is called a ciphertext
refresh.

Gentry's key idea to refresh a ciphertext is to homomorphically
evaluate the decryption circuit on the ciphertext bits, using an
encryption of the secret-key bits. This is called bootstrapping. 
Then instead of getting the  plaintext bit (as we would get if we 
would evaluate the decryption circuit with the secret-key bits in
clear), we get an __encryption__ of the plaintext bit, i.e. a new
ciphertext for the same plaintext. Now if the decryption circuit has a
small enough depth, then the amount of noise in the new ciphertext can
be actually smaller than in the original plaintext, hence a ciphertext refresh.

However the previous decryption algorithm m=(c mod p) mod 2 does not
have a small depth, therefore the decryption procedure must be "squashed" so
that it can be expressed as a low depth circuit. How this is done in
explained in [2].

So far we have only described a secret-key encryption scheme, i.e. to
encrypt one must know the secret-key p. However it is easy to obtain a
public-key encryption scheme. For this generate a public set of
ciphertexts xi which are all different encryptions of 0: 

xi=qi*p+2*ri

and to encrypt a bit m compute:

c=m+2*r+random_subset_sum(xi)

It is easy to see that c is indeed an encryption of the bit m.

WHAT IS IMPLEMENTED ?
---------------------

We provide an implementation of the DGHV scheme with the fully
homomorphic capability, i.e. we implement the key generation,
encryption, decryption, add, mult and ciphertext refresh procedures.

The implementation is done using the Sage 4.7.2 mathematical library
under Python, available at http://www.sagemath.org/

First run sage, then type:

sage: load "dghv.sage"
sage: testPkRecrypt()

The testPkRecrypt() function demonstrates key generation, addition and
multiplication of ciphertexts, and ciphertext refresh. It runs for
four sets of parameters (toy, small, medium and large), as described
in [1].


WHAT IS A COMPRESSED PUBLIC-KEY ?
---------------------------------

In DGHV the ciphertext size must be huge to prevent attacks based on
lattice reduction algorithms; at least 10^7 bits. The secret p is
comparatively smaller, roughly 2000 bits. Since roughly 10^4
ciphertexts xi must be included in the public-key, that would give a
public-key size of 10^11 bits, i.e. 12.5 GB.

To reduce the public-key size we implement the following technique
described in [1]. Instead of generating xi as xi=qi*p+2*ri, one first
generates a pseudo-random Xi of the same size, and computes a small
correction di such that xi=Xi-di is small modulo p. Then only these
small corrections need to be stored in the public-key, with the seed
of the PRNG. Using this compression technique the public-key size
becomes 10^4*(2*10^3)=2*10^7 bits, i.e. 2.5 MB, which is more
manageable.
