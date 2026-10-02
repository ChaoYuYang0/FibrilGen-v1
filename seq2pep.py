# run seq2pep.py

import numpy as np
import pymol
import math

def get_b_strand(seq):
	pymol.cmd.fab(seq,'pep')
	for i in range(len(seq)-1):
		# Set psi angle
		pymol.cmd.select('po1','resi %s and name n'%(i+1))
		pymol.cmd.select('po2','resi %s and name ca'%(i+1))
		pymol.cmd.select('po3','resi %s and name c'%(i+1))
		pymol.cmd.select('po4','resi %s and name n'%(i+2))
		pymol.cmd.set_dihedral('po1','po2','po3','po4',135)
		# Set phi angle
		pymol.cmd.select('po1','resi %s and name c'%(i+1))
		pymol.cmd.select('po2','resi %s and name n'%(i+2))
		pymol.cmd.select('po3','resi %s and name ca'%(i+2))
		pymol.cmd.select('po4','resi %s and name c'%(i+2))
		pymol.cmd.set_dihedral('po1','po2','po3','po4',221)
	pymol.cmd.delete('po* or pk*')
	pymol.cmd.hide('all')
	pymol.cmd.show('sticks')
	return 

def set_chi_angle(resi,a1,a2,a3,a4,angle):
	# Set chi angle
	pymol.cmd.select('po1','resi %s and name %s'%(resi,a1))
	pymol.cmd.select('po2','resi %s and name %s'%(resi,a2))
	pymol.cmd.select('po3','resi %s and name %s'%(resi,a3))
	pymol.cmd.select('po4','resi %s and name %s'%(resi,a4))
	pymol.cmd.set_dihedral('po1','po2','po3','po4',angle)
	pymol.cmd.delete('po* or pk*')
	return