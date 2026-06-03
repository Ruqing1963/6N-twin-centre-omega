#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omega_{>3} statistics of twin centres vs ordinary integers (S10).
Confirms three things:
 (1) per-prime: P(q|N | N twin) = 1/(q-2) exactly for small q (ratio -> 1),
     versus 1/q for ordinary integers -- the microscopic root of the Part I
     enrichment factor prod (q-1)/(q-3).
 (2) mean shift <omega>_twin - <omega>_ord, measured unbiasedly (full
     factorisation), compared to the convergent closed form sum_{q>3}(1/(q-2)-1/q)
     ~ 0.2604. The measured global shift is slightly smaller; the deficit comes
     from large-q terms truncated by finite N (each contributes ~2/q^2). With S10's
     larger N range the measured shift should move toward 0.2604.
 (3) normality diagnostics (mean, var, skew, excess kurtosis). NOTE: in this range
     lnln(6N) ~ 3 is small, so omega_{>3} has NOT converged to the Erdos-Kac normal
     (mean != var); the limiting distribution is left as an open problem. We report
     the diagnostics honestly rather than claim normality.
Default S10. Requires: numpy.
"""
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
QLADDER=[q for q in [5,7,11,13,17,19,23,29,31,37,41,43,47,97,199,401,1009,2003,
                     5003,10007,100003,1000003,10000019] if q<=HI]
div_tw={q:0 for q in QLADDER}; n_tw=0
# moments
sa=ssa=m3a=m4a=ca=0.0
st=sst=m3t=m4t=ct=0.0
import time; t0=time.time()
n=LO
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int32)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3: ob[idx]+=1
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    obf=ob.astype(np.float64)
    sa+=obf.sum(); ssa+=(obf**2).sum(); m3a+=(obf**3).sum(); m4a+=(obf**4).sum(); ca+=sz
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st2=max(p*p,((vlo+p-1)//p)*p)
        if st2>vhi: continue
        comp[st2-vlo:span:p]=True
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    Ntw=Narr[np.nonzero(tw)[0]]; obt=obf[np.nonzero(tw)[0]]
    st+=obt.sum(); sst+=(obt**2).sum(); m3t+=(obt**3).sum(); m4t+=(obt**4).sum(); ct+=len(obt)
    for q in QLADDER: div_tw[q]+=int((Ntw%q==0).sum())
    n_tw+=len(Ntw)
    n=nh
print(f"S{MAXK}: ordinary {int(ca):,}, twins {n_tw:,}; scan {time.time()-t0:.0f}s")
def moments(s,ss,m3,m4,c):
    mean=s/c; var=ss/c-mean**2
    mu3=m3/c-3*mean*ss/c+2*mean**3; mu4=m4/c-4*mean*m3/c+6*mean**2*ss/c-3*mean**4
    return mean,var,mu3/var**1.5,mu4/var**2-3
ma,va,ska,ka=moments(sa,ssa,m3a,m4a,ca)
mt,vt,skt,kt=moments(st,sst,m3t,m4t,ct)
print(f"\n(1) per-prime P(q|N|twin) vs 1/(q-2)  [ordinary 1/q]:")
print(f"{'q':>9}{'meas':>12}{'1/(q-2)':>12}{'ratio':>9}")
for q in QLADDER:
    pm=div_tw[q]/n_tw
    print(f"{q:>9}{pm:>12.7f}{1/(q-2):>12.7f}{pm*(q-2):>9.4f}")
qs=primes_upto(10**7); qs=qs[qs>3].astype(float)
closed=np.sum(1/(qs-2)-1/qs)
print(f"\n(2) mean shift: measured {mt-ma:+.4f}  vs closed-form sum(1/(q-2)-1/q)={closed:.4f}")
print(f"    (S9 measured was +0.2220; if deficit is finite-N tail, S10 should be larger)")
print(f"    ordinary mean {ma:.4f}, twin mean {mt:.4f}")
print(f"\n(3) normality diagnostics (NOT converged in this range; reported honestly):")
print(f"    {'':>9}{'mean':>9}{'var':>9}{'skew':>8}{'exkurt':>8}")
print(f"    {'ordinary':>9}{ma:>9.4f}{va:>9.4f}{ska:>8.3f}{ka:>8.3f}")
print(f"    {'twin':>9}{mt:>9.4f}{vt:>9.4f}{skt:>8.3f}{kt:>8.3f}")
print(f"    mean != var and skew,kurt != 0 => Erdos-Kac limit not yet reached (open).")

# ---- emit CSV (omega_twin_S{K}.csv: per-q probabilities + summary) ----
import csv as _csv
qs=primes_upto(10**7); qs=qs[qs>3].astype(float); closed=float(np.sum(1/(qs-2)-1/qs))
with open(f'omega_perq_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['q','P_ord','P_tw','inv_qm2','ratio'])
    for q in QLADDER:
        pm=div_tw[q]/n_tw
        _w.writerow([q,f'{1/q:.7f}',f'{pm:.7f}',f'{1/(q-2):.7f}',f'{pm*(q-2):.4f}'])
with open(f'omega_summary_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f)
    _w.writerow(['shell','ord_mean','twin_mean','mean_shift','closed_form_shift','ord_var','twin_var','twin_skew','twin_exkurt'])
    _w.writerow([f'S{MAXK}',f'{ma:.4f}',f'{mt:.4f}',f'{mt-ma:.4f}',f'{closed:.4f}',f'{va:.4f}',f'{vt:.4f}',f'{skt:.3f}',f'{kt:.3f}'])
print(f"\n[ok] wrote omega_perq_S{MAXK}.csv and omega_summary_S{MAXK}.csv")
