#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 2-panel omega figure from ../data. Left: P(q|twin)/(1/(q-2)) vs q for
S9 and S10 (flat at 1 for small q). Right: omega_{>3} distribution, twin vs
ordinary (S9 histogram). Reads omega_perq_data.csv and omega_dist_data.csv.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
perq=list(csv.DictReader(open('../data/omega_perq_data.csv')))
dist=list(csv.DictReader(open('../data/omega_dist_data.csv')))
fig,axes=plt.subplots(1,2,figsize=(14,5.4))
q=np.array([int(r['q']) for r in perq])
rS10=np.array([float(r['P_tw_S10'])/float(r['inv_qm2']) for r in perq])
rS9=np.array([float(r['P_tw_S9'])/float(r['inv_qm2']) for r in perq])
axes[0].axhspan(0.99,1.01,color='#2ca25f',alpha=.12,label='$\\pm1\\%$')
axes[0].plot(q,rS9,'o-',color='#185FA5',lw=1.4,ms=6,label='S9',alpha=.7)
axes[0].plot(q,rS10,'s-',color='#c0392b',lw=1.6,ms=7,label='S10')
axes[0].axhline(1,color='gray',ls=':',lw=1.2); axes[0].set_xscale('log')
axes[0].set_xlabel('prime $q$ (log scale)',fontsize=11)
axes[0].set_ylabel(r'$P(q\mid N,\ \mathrm{twin})\ /\ (1/(q-2))$',fontsize=11)
axes[0].set_title('$P(q\\mid N,\\,\\mathrm{twin})=1/(q-2)$ exactly for small $q$',fontsize=12)
axes[0].set_ylim(0.96,1.04); axes[0].legend(fontsize=10); axes[0].grid(alpha=.25,which='both')
k=np.array([int(r['k']) for r in dist])
po=np.array([float(r['ord_S9']) for r in dist]); pt=np.array([float(r['twin_S9']) for r in dist])
w=0.38
axes[1].bar(k-w/2,po,w,color='#185FA5',label='ordinary integers',alpha=.85)
axes[1].bar(k+w/2,pt,w,color='#c0392b',label='twin centres',alpha=.85)
axes[1].set_xlabel(r'$\omega_{>3}(N)$',fontsize=11); axes[1].set_ylabel('probability',fontsize=11)
axes[1].set_title('$\\omega_{>3}$ distribution: twins shifted right\n(mean $+0.22$, from $1/q\\to1/(q-2)$)',fontsize=12)
axes[1].legend(fontsize=10); axes[1].grid(alpha=.25,axis='y'); axes[1].set_xticks(k)
plt.suptitle('Twin-centre $\\omega_{>3}$ statistics: each prime divides a twin centre with probability $1/(q-2)$, not $1/q$',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper12_omega.pdf',bbox_inches='tight')
plt.savefig('fig_paper12_omega.png',dpi=160,bbox_inches='tight')
print("figure saved")
