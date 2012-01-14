/* 
* directscale.c: faster implementation of ciphertext expand
* 
* Copyright (c) 2012 Mehdi Tibouchi <mehdi.tibouchi@normalesup.org>
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License version 2 as published
* by the Free Software Foundation.
*/

#include <gmp.h>

#define w		(GMP_NUMB_BITS/2)
#define _BOTMASK	((1ul << w)-1)
#define _TOPMASK	(~_BOTMASK)
#define BOT(x)		((x) & _BOTMASK)
#define TOP(x)		((x) >> w)
#define LIMB(z,i)	(((i)<((z)->_mp_size))?((z)->_mp_d[i]):(0L))
#define BOTL(z,i)	(BOT(LIMB(z,i)))
#define TOPL(z,i)	(TOP(LIMB(z,i)))
#define HLIMB(z,j)	((j&1)?(TOPL(z,j>>1)):(BOTL(z,j>>1)))

unsigned getGMP_NUMB_BITS()
{
  return GMP_NUMB_BITS;
}

unsigned long directScal(unsigned long kap, mpz_t cz, mpz_t yz)
{
  unsigned long nW=(kap+1)/(2*w), val=0, i;

  if(nW*w*2 != kap+1)
      return 0;

  for(i = 0; i < nW-1; i++) {
    val += BOTL(cz,i) * LIMB(yz,nW-1-i);
    val += (BOTL(cz,i) * TOPL(yz,nW-2-i)) >> w;
    val += TOPL(cz,i) * ((BOTL(yz,nW-1-i) << w) + TOPL(yz,nW-2-i));
    val += (TOPL(cz,i) * BOTL(yz,nW-2-i)) >> w;
  }

  val += BOTL(cz,nW-1) * LIMB(yz,0);
  val += TOPL(cz,nW-1) * (BOTL(yz,0) << w);

  return val;
}

