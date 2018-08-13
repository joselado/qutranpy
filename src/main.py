#!/usr/bin/python

from __future__ import print_function
import matplotlib.pyplot as plt


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1])
  xmlpath = ""
else: # default path
  raise


def plotbands():
  """Function to plot the band structure"""
  os.system("plotbands  &")

def plotlandauer(name):
  """Function to plot the band structure"""
  os.system("plotlandauer  "+name+"  &")

def plotadiabatic(name):
  """Function to plot the band structure"""
  os.system("plotadiabatic  "+name+"  &")




from gi.repository import Gtk as gtk
builder = gtk.Builder()

from qh_interface import * # import all the libraries needed


def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())


def modify(name,value,active=False):
  """Modify a certain text entry"""
  entry = builder.get_object(name)
  entry.set_text(str(value))
  if active: entry.activate() # set active









def get_geometry():
  """Get geometry used"""
  lattice_name = builder.get_object("geometry").get_active_text()
  if lattice_name=="Honeycomb zigzag":
    geometry_builder = geometry.honeycomb_zigzag_ribbon
  elif lattice_name=="Honeycomb armchair":
    geometry_builder = geometry.honeycomb_armchair_ribbon
  elif lattice_name=="Square":
    geometry_builder = geometry.square_tetramer_ribbon
  elif lattice_name=="Chain":
    geometry_builder = geometry.chain
  g = geometry_builder(int(get("width"))) # get the unit cell
  return g # return geometry



def update_mxyz2m(self):
  """Update different components taking as inputs mxyz"""
  for name in ["L","C","R"]: # loop over leads
    jx = get("Jx_"+name) # get component 
    jy = get("Jy_"+name) # get component
    jz = get("Jz_"+name) # get component
    jj = np.sqrt(jx*jx + jy*jy + jz*jz) # absolute value
    theta = np.arctan2(np.sqrt(jx*jx+jy*jy),jz)/np.pi # theta angle
    phi = np.arctan2(jy,jx)/np.pi # phi angle
    jj = round(jj,6)
    theta = round(theta,6)
    phi = round(phi,6)
    modify("JJ_"+name,jj)
    modify("theta_"+name,theta)
    modify("phi_"+name,phi)


def update_m2mxyz(self):
  """Update different components taking as inputs angles"""
  for name in ["L","C","R"]: # loop over leads
    jj = get("JJ_"+name) # get component 
    theta = get("theta_"+name) # get component
    phi = get("phi_"+name) # get component
    jz = jj*np.cos(theta*np.pi)
    jx = jj*np.sin(theta*np.pi)*np.cos(phi*np.pi)
    jy = jj*np.sin(theta*np.pi)*np.sin(phi*np.pi)
    jx = round(jx,6)
    jy = round(jy,6)
    jz = round(jz,6)
    modify("Jx_"+name,jx)
    modify("Jy_"+name,jy)
    modify("Jz_"+name,jz)






def get_hamiltonian(name="C"):
  """Get a certain hamiltonian""" 
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.shift_fermi(get("fermi"+"_"+name)) # get fermi energy
  zeeman = [get("Jx"+"_"+name),get("Jy"+"_"+name),get("Jz"+"_"+name)]
  h.add_zeeman(zeeman) # add zeeman field
  h.add_rashba(get("rashba"+"_"+name)) # get fermi energy
  h.add_sublattice_imbalance(get("mass"+"_"+name)) # get fermi energy
  h.add_antiferromagnetism(get("maf"+"_"+name)) # get fermi energy
  h.add_peierls(get("peierls"+"_"+name)) # get fermi energy
  if builder.get_object("has_eh").get_active():
    h.add_swave(get("super"+"_"+name)) # get fermi energy
  if not builder.get_object("has_spin").get_active():
    h.remove_spin()
  # scale all the energy scales
  h.intra *= get("hopping_"+name)
  h.inter *= get("hopping_"+name)
  return h



def get_heterostructure():
  """Get a heterostructure"""
  # get the three hamiltonians
  while gtk.events_pending(): # don't freeze the interface
    gtk.main_iteration_do(False) # update the interface
  hl = get_hamiltonian("L")
  hc = get_hamiltonian("C")
  hr = get_hamiltonian("R")
  nc = int(get("length"))
  print("Central part has",nc,"cells")
  hlist = [hc for i in range(nc)] # central part
  ht = heterostructures.create_leads_and_central_list(hr,hl,hlist)
  # modify couplings to the leads
  ht.left_coupling *= get("coupling_L")
  ht.right_coupling *= get("coupling_R")
  # finish
  return ht


def get_operator_bands(h,name):
  """Get the operator for the bands"""
  part = builder.get_object(name).get_active_text()
  if part=="Sx": return operators.get_sx(h)
  if part=="Sy": return operators.get_sy(h)
  if part=="Sz": return operators.get_sz(h)
  return None

def get_bands_right(self):
  h = get_hamiltonian(name="R")
  op = get_operator_bands(h,"operator_R")
  h.get_bands(operator=op)
  os.system("tb90-bands &")



def get_bands_central(self):
  h = get_hamiltonian(name="C")
  op = get_operator_bands(h,"operator_C")
  h.get_bands(operator=op)
  os.system("tb90-bands &")



def get_bands_left(self):
  h = get_hamiltonian(name="L")
  op = get_operator_bands(h,"operator_L")
  h.get_bands(operator=op)
  os.system("tb90-bands &")


##################################
# different types of sweeps
################################
def get_name_parameter(nameparam,lead):
  namep = builder.get_object(nameparam).get_active_text() # active name
  if namep=="Magnetic field": name = "peierls" 
  elif namep=="Exchange x": name = "Jx" 
  elif namep=="Exchange y": name = "Jy" 
  elif namep=="Exchange z": name = "Jz" 
  elif namep=="Total exchange": name = "JJ" 
  elif namep=="Theta exchange": name = "theta" 
  elif namep=="Phi exchange": name = "phi" 
  elif namep=="Fermi": name = "fermi" 
  elif namep=="Rashba": name = "rashba" 
  elif namep=="Sublattice imbalance": name = "mass" 
  elif namep=="SC pairing": name = "super" 
  elif namep=="Hopping": name = "hopping" 
  elif namep=="Lead coupling": name = "coupling" 
  else: raise
  name += "_" # add to the name
  namep = builder.get_object(lead).get_active_text() # active lead
  if namep=="Left lead": name += "L" 
  elif namep=="Right lead": name += "R" 
  elif namep=="Central": name += "C" 
  else: raise
  return name # return name



def run_hamiltonian(self):
  ht = get_heterostructure() # get the current heterostructure
  # name of the parameter in the interface
  name = get_name_parameter("parameter_ham","which_lead") # get label of param
  # now name for the axis
  nameaxis = builder.get_object("parameter_ham").get_active_text() 
  nameaxis +=" "+ builder.get_object("which_lead").get_active_text() 
  # start the loop over values
  param_ham = np.linspace(get("min_ham"),get("max_ham"),
                            get("step_ham")) # get the energies
  Gs = [] # empty list
  # now modify that parameter accordingly
  for p in param_ham:
    modify(name,p,active=True) # modify the parameter
    ht = get_heterostructure() # get the current heterostructure
    Gs.append(heterostructures.didv(ht,energy = get("energy")))
  namef = "TRANSPORT_PARAMETER.OUT" # name of the output file
  inout.write(param_ham,Gs,output_file=namef,comment="xaxis = "+nameaxis)
  plotlandauer(namef)
#  os.system("qh-landauer  "+namef+ "  &")
  # nor the adiabatic part
  if builder.get_object("do_adiabatic").get_active():
    At = [] # initialize list
    for p in param_ham:
      modify(name,p,active=True) # modify energy parameter
      At.append(get_adiabatic()) # add the adibatic part
    namef = "ADIABATIC_PARAMETER.OUT" # name of the output file
    inout.write(param_ham,At,output_file=namef,comment="xaxis = "+nameaxis)
    plotadiabatic(namef)
#    os.system("qh-adiabatic  "+namef+ "  &")




def run_energy(self):
  ht = get_heterostructure() # get the current heterostructure
  energies = np.linspace(get("min_energy"),get("max_energy"),
                            get("step_energy")) # get the energies
  Gs = [heterostructures.didv(ht,energy = e) for e in energies] # get transport
  name = "LANDAUER_ENERGY.OUT" # name of the output file
  inout.write(energies,Gs,output_file=name,comment="xaxis = Energy")
  plotlandauer(name)
#  os.system("qh-landauer  "+name+ "  &")
  # now adiabatic part
  if builder.get_object("do_adiabatic").get_active():
    At = [] # initialize list
    for e in energies:
      modify("energy",e) # modify energy parameter
      At.append(get_adiabatic()) # add the adibatic part
    name = "ADIABATIC_ENERGY.OUT" # name of the output file
    inout.write(energies,At,output_file=name,comment="xaxis = Energy")
    plotadiabatic(name)
#    os.system("qh-adiabatic  "+name+ "  &")




def run_length(self):
  lengths = range(int(get("min_length")),int(get("max_length")),
                     int(get("step_length"))) # get the energies
  Gs = [] # empty list
  for l in lengths: # loop over lengths
    modify("length",l) # modify the length
    ht = get_heterostructure() # get the current heterostructure
    Gs.append(heterostructures.didv(ht,energy = get("energy")))
  name = "TRANSPORT_LENGTHS.OUT" # name of the output file
  inout.write(lengths,Gs,output_file=name,comment="xaxis = Central region length")
  plotlandauer(name)
#  os.system("qh-landauer  "+name+ "  &")



def run_disorder(self):
  ht = get_heterostructure() # get the current heterostructure



###########################
# now adiabatic transport #
###########################


def get_adiabatic():
  """Calcualte derivative of smatrix"""
  # names of the modified parameters
  name1 = get_name_parameter("parameter_ham_adiabatic1","which_lead_adiabatic1")
  name2 = get_name_parameter("parameter_ham_adiabatic2","which_lead_adiabatic2")
  # now perform the derivative
  eps = 0.001 # value for the derivative
  energy = get("energy") # get energy
  def get_s(ht):
    return heterostructures.get_smatrix(ht,energy=energy,
                                       as_matrix=True,check=True)
  def get_ds(name):
    modify(name,get(name)-eps,active=True) 
    hm = get_heterostructure()
    modify(name,get(name)+2*eps,active=True) 
    hp = get_heterostructure()
    modify(name,get(name)-eps,active=True) # original value 
    sm = get_s(hm) # get smatrix
    sp = get_s(hp) # get smatrix
    return (sp-sm)/(2.*eps) # return derivative
  s1 = get_ds(name1)
  s2 = get_ds(name2)
  c12 = np.conjugate(s1)*s2.T
  value = np.sum([c12[i,i] for i in range(c12.shape[0]/2)]).imag
  return value # return adiabatic pumping
  # different smatrices


###########################
###########################


# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["bands_right"] = get_bands_right  # initialize and run
signals["bands_left"] = get_bands_left  # initialize and run
signals["bands_central"] = get_bands_central  # initialize and run
signals["run_hamiltonian"] = run_hamiltonian  # initialize and run
signals["run_energy"] = run_energy  # initialize and run
signals["run_length"] = run_length  # initialize and run
signals["run_disorder"] = run_disorder  # initialize and run

# update signals
signals["update_mxyz2m"] = update_mxyz2m  # initialize and run
signals["update_m2mxyz"] = update_m2mxyz  # initialize and run

class RibbonApp(object):       
	def __init__(self):
	    builder.add_from_file(xmlpath+"main.xml")
	    builder.connect_signals(signals)
	    self.window = builder.get_object("main")
# creates a temporal folder in /tmp
#            folder = create_folder()  
	    self.window.show()



#if __name__ == "__main__":
app = RibbonApp()  # launch the app
gtk.main()  # infinite loop

