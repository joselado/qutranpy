import multicell
import sculpt

def build(h,nz=1):
  """Create Hamiltonian of a film from a 3d geometry"""
#  if not h.dimensionality==3: raise
  ho = multicell.supercell(h,nsuper=[1,1,nz],sparse=False,ncut=3)
  ho.dimensionality = 2 # reduce dimensionality
  ho.geometry.dimensionality = 2 # reduce dimensionality
  ho.geometry = sculpt.set_xy_plane(ho.geometry) # put in the xy plane
  hoppings = [] # empty list
  for t in ho.hopping: # loop over hoppings
    if t.dir[2]== 0: hoppings.append(t) # remove in the z direction
  ho.hopping = hoppings # overwrite hoppings
  return ho



def geometry_film(g,nz=1):
  """Create the geometry of a film"""
  go = g.supercell([1,1,nz]) # create the supercell
  go.dimensionality = 2 # reduce dimensionality
  go = sculpt.set_xy_plane(go) # put in the xy plane
  return go


