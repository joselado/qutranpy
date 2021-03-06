import numpy as np
import klist
import scipy.linalg as lg


def ipr2d(h,nk=10,random=True,window=[-0.1,0.1]):
  """Calculate the inverse participation ratio"""
  hkgen = h.get_hk_gen() # get generator
  ks = range(nk*nk)
  out = 0.0 # start
  for k in ks:
    k = np.random.random(2) # random k-point
    hk = hkgen(k) # get hamiltonian
    es,ws = lg.eigh(hk) # diagonalize
    ws = ws.transpose() # transpose
    nume = 0
    tmp = 0.0
    for (e,w) in zip(es,ws):
      if window[0]<e<window[1]:
        nume += 1
        d = w*np.conjugate(w) # density
        d2 = (d*d).real # ipr
        d2 = np.sum(d2) # ipr
        tmp += d2 # store
    if nume>0: tmp/= nume # normalize
    out += tmp # add
  return out/len(ks)
        
