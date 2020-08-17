// Time-stamp: <2008-08-21 16:25:21 zaykind> (Dmitri Zaykin)
//
// Uniform random numbers generators. Based on Marsaglia's
// generators posted to sci.stat.math (with my fixing of 
// the bug in the cyclic table)
//

#ifndef DZ_MARSA_H
#define DZ_MARSA_H

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#define SEMIOPEN01  2.3283064365386963e-10 /* [0, 1) */
#define CLOSED01    2.3283064370807974e-10 /* [0, 1] */
#define ORIGINAL01  2.328306e-10

class Rndm
{
 public:
    typedef /* DZ: WAS: unsigned long*/ uint32_t Unlo;
    typedef unsigned char Uc;
 private:
    Unlo z, w, jsr, jcong, t[UCHAR_MAX+1], x, y, a, b;
    unsigned char c;
    Unlo znew() { return (z = 36969UL*(z & 65535UL)+(z >> 16)); }
    Unlo wnew() { return (w = 18000UL*(w & 65535UL)+(w >> 16)); }
 public:
    Rndm() {}
    virtual ~Rndm() {}
    Rndm(const Rndm& src) {
        z=src.z, w=src.w, jsr=src.jsr, jcong=src.jcong;
        for(size_t i=0; i<UCHAR_MAX+1; i++) t[i]=src.t[i];
        x=src.x, y=src.y, a=src.a, b=src.b;
        c=src.c;
    }
    Rndm (Unlo i1, Unlo i2, Unlo i3, Unlo i4, Unlo i5, Unlo i6) {
        Init (i1, i2, i3, i4, i5, i6);
    }
    Rndm (Unlo i1) {
        Init (i1);
    }
    virtual void Init(Unlo i1, Unlo i2, Unlo i3, Unlo i4, Unlo i5, Unlo i6) {
        z=i1, w=i2, jsr=i3, jcong=i4, x=0, y=0, a=i5, b=i6, c=0;
        for(size_t i=0; i<UCHAR_MAX+1; i++) t[i] = Fib();
    }
    virtual void Init(Unlo i1) {
        z=i1, w=i1, jsr=i1, jcong=i1, x=0, y=0, a=i1, b=i1, c=0;
        for(size_t i=0; i<UCHAR_MAX+1; i++) t[i] = Fib();
    }
    Unlo Mwc() { return (znew() << 16) + wnew(); }
    Unlo Shr3 () {
        jsr=jsr^(jsr<<17);
        jsr=jsr^(jsr>>13); 
        return (jsr=jsr^(jsr<<5));
    }
    Unlo Cong() { return (jcong = 69069UL*jcong + 1234567UL); }
    Unlo Kiss() { return (Mwc() ^ Cong()) + Shr3(); }
    Unlo Swb () {
        x = t[(Uc)(c+15)];
        t[(Uc)(c+237)] = x - (y = t[(Uc)(c+1)] + (x < y));
        return t[++c];
    }
    Unlo Lfib4() {
        t[c]=t[c]+t[(Uc)(c+58)]+t[(Uc)(c+119)]+t[(Uc)(c+179)];
        return t[++c];
    }
    Unlo Fib() { b=a+b; return (a=b-a); }
    double Closed01 () { return (Kiss()+Swb()) * CLOSED01; }
    double Op () { return (Kiss()+Swb()) * SEMIOPEN01; }
    double Ori () { return (Kiss()+Swb()) * ORIGINAL01; }
    double Vni() { return long(Kiss()+Swb()) * 4.656613e-10; }
    double operator () () { return Op(); }
    Unlo operator () (Unlo n) {
        return int(Op() * n);
    }
    double operator () (double Min, double Max) { return Min+Op()*(Max-Min); }
};

#endif
