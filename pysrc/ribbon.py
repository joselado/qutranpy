from __future__ import print_function,division
import sculpt
import numpy as np
import supercell

def bulk2ribbon(g,boundary=[1,0],n=10,clean=True):
  """Return the geometry of a ribbon"""
  go = g.copy() # copy
  m = [[boundary[0],boundary[1],0],[0,1,0],[0,0,1]] # supercell
  if boundary[0]!=1 or boundary[1]!=0:
    go = supercell.non_orthogonal_supercell(go,m,mode="brute") 
  go = go.supercell((1,n)) # create a supercell
  go = sculpt.rotate_a2b(go,go.a1,np.array([1.,0.,0.]))
  go.dimensionality = 1 # zero dimensional
  go.a2 = np.array([0.,1.,0.])
  if clean: 
    go = go.clean(iterative=True)
    if len(go.r)==0:
      print("Ribbon is not wide enough")
      raise
  go.real2fractional()
  go.fractional2real()
  go.celldis = go.a1[0]
  go.center()
  return go


