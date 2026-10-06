# run pep2unit.py

import numpy as np
import pymol
import math

def get_ca(name):
	pymol.cmd.select('lo_ca','name ca and '+name)
	pos_ca = pymol.cmd.get_coords('lo_ca',1)
	pymol.cmd.delete('*ca*')
	return np.array(pos_ca)

def get_boundary(lo_name):
	lop = []
	for name in lo_name:
		lop += pymol.cmd.get_coords(name,1).tolist()
	lop = np.array(lop)
	return [max(lop[:,0]),min(lop[:,0]),max(lop[:,1]),min(lop[:,1]),max(lop[:,2]),min(lop[:,2])]

def rotate_coordinate(name,coord):
	mf = coord[0].tolist()+[0]+coord[1].tolist()+[0]+coord[2].tolist()+[0]+[0,0,0,1]
	pymol.cmd.transform_selection(name,mf)

def example():
	## Reset
	pymol.cmd.delete('all')
	pymol.cmd.reset()
	## Load PDB
	pymol.cmd.load('structures/input/F8.pdb')

	## Create a periodic unit
	# INPUT: (peptide, start residue index, end residue index, does the start residue locate at the face)
	unit = create_pep_unit('F8',2,8,0)

	## --- Examples of sheet structures ---
	# INPUT: (peptide unit)
	sheet = create_sheet(unit)
	# INPUT: (peptide alignment aaa/apa/aap/app/paa/ppa/pap/ppp, does two beta-sheet aligned face-to-face, num of units per sheet)
	sheet.build_a_plain_sheet('ppp',1,5)
	## --- End ---
	
	sheet.get_dimension()
	pymol.cmd.zoom()


class create_pep_unit():
	# INPUT (molecule name, start residue index on the linear segment, end residue index on the linear segment)
	def __init__(self,name,resi_start,resi_end,is_face):
		# Create a cooridinate
		coord = self.get_coordinate(resi_start,resi_end,is_face)
		# Align the unit to the coordinate
		rotate_coordinate(name,coord)
		self.name,self.resi_start,self.resi_end = name,resi_start,resi_end
		self.set_boundary()

	def set_x_by_ca(self,resi_start,resi_end):
		resi_end_even_space = (resi_end-resi_start)//2*2+resi_start
		pymol.cmd.select('x_po1',f'resi {resi_start} and name ca')
		pymol.cmd.select('x_po2',f'resi {resi_end_even_space} and name ca')
		pos_x_po1 = np.array(pymol.cmd.get_coords('x_po1',1)).reshape(3)
		pos_x_po2 = np.array(pymol.cmd.get_coords('x_po2',1)).reshape(3)
		pymol.cmd.delete('*po*')
		return pos_x_po2 - pos_x_po1

	def set_y_by_nh(self,resi_start,resi_end,is_face):
		acc,vec_acc = 1,0
		if is_face:
			resi_start_from_face = resi_start+1
		else:
			resi_start_from_face = resi_start
		for idx_i in range(resi_start_from_face,resi_end+1,2):
			pymol.cmd.select(f'y{acc}_po1',f'resi {idx_i} and name n')
			pymol.cmd.select(f'y{acc}_po2',f'resi {idx_i} and name h')
			pos_yi_po1 = np.array(pymol.cmd.get_coords(f'y{acc}_po1',1)).reshape(3)
			pos_yi_po2 = np.array(pymol.cmd.get_coords(f'y{acc}_po2',1)).reshape(3)
			y_i = pos_yi_po2 - pos_yi_po1
			vec_acc += y_i
			acc += 1
			pymol.cmd.delete('*po*')
		return vec_acc/acc

	def get_coordinate(self,resi_start,resi_end,is_face):
		x = self.set_x_by_ca(resi_start,resi_end)
		y = self.set_y_by_nh(resi_start,resi_end,is_face)
		z = np.cross(x,y)
		y = np.cross(z,x)
		unit_x,unit_y,unit_z = x/np.linalg.norm(x),y/np.linalg.norm(y),z/np.linalg.norm(z)
		return [unit_x,unit_y,unit_z]

	def set_boundary(self):
		# Set dimensions
		boundary = get_boundary([self.name])
		self.y = 4.8 # Set peptide to peptide dist as 0.48 nm 
		self.z = boundary[4]-boundary[5]


class create_sheet():
	def __init__(self,unit): 
		self.unit = unit
		self.ca_start = self.map_idx_res2ca(self.unit.resi_start)
		self.ca_end = self.map_idx_res2ca((self.unit.resi_end-self.unit.resi_start)//2*2+self.unit.resi_start)

	def map_idx_res2ca(self,idx):
		pymol.cmd.select('all_ca',f'name ca')
		idx_all_ca = pymol.cmd.index('all_ca')
		pymol.cmd.select('this_ca',f'resi {idx} and name ca')
		idx_this_ca = pymol.cmd.index('this_ca')
		idx_new = idx_all_ca.index(idx_this_ca[0])
		pymol.cmd.delete('*ca*')
		return idx_new

	def build_a_plain_sheet(self,b_alignment,s_alignment,num_half):
		b_flip_angle = self.get_b_flip_angle(b_alignment)
		s_flip_angle = self.get_s_flip_angle(s_alignment)
		z_offset = self.unit.z/2.0
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s1_pep1_'+str(i),'p_s1_pep2_'+str(i)
			pymol.cmd.create(name_pep1,self.unit.name)
			self.affine_transformation_backbone(name_pep1,[0,0])
			pymol.cmd.translate([0,self.unit.y*2*i,-z_offset],name_pep1)
			pymol.cmd.create(name_pep2,self.unit.name)
			self.affine_transformation_backbone(name_pep2,b_flip_angle[0])
			pymol.cmd.translate([0,self.unit.y*(2*i+1),-z_offset],name_pep2)
		for i in range(num_half):
			name_pep1,name_pep2 = 'p_s2_pep1_'+str(i),'p_s2_pep2_'+str(i)
			pymol.cmd.create(name_pep1,self.unit.name)
			self.affine_transformation_sidechain(name_pep1,s_flip_angle)
			self.affine_transformation_backbone(name_pep1,b_flip_angle[1])
			pymol.cmd.translate([0,self.unit.y*2*i,z_offset],name_pep1)
			pymol.cmd.create(name_pep2,self.unit.name)
			self.affine_transformation_sidechain(name_pep2,s_flip_angle)
			self.affine_transformation_backbone(name_pep2,b_flip_angle[2])
			pymol.cmd.translate([0,self.unit.y*(2*i+1),z_offset],name_pep2)
		pymol.cmd.color('green','p_s1_*')
		pymol.cmd.color('orange','p_s2_*')
		pymol.cmd.group('plain_sheet','p_*')
		lon = ['p_s1_pep1_'+str(i) for i in range(num_half)]+['p_s1_pep2_'+str(i) for i in range(num_half)]+\
				['p_s2_pep1_'+str(i) for i in range(num_half)]+['p_s2_pep2_'+str(i) for i in range(num_half)]
		self.set_dimension(lon)
		return

	def set_dimension(self,lon):
		boundary = get_boundary(lon)
		self.x = boundary[0]-boundary[1]
		self.y = boundary[2]-boundary[3]
		self.z = boundary[4]-boundary[5]
		return

	def get_dimension(self):
		print ('Dimension x: '+str(self.x)+' nm')
		print ('Dimension y: '+str(self.y)+' nm')
		print ('Dimension z: '+str(self.z)+' nm')
		return

	def get_b_flip_angle(self,b_alignment):
		# assign flip angle for [s1_pep2,s2_pep1,s2_pep2]
		if b_alignment=='aaa':
			return [[0,180],[180,0],[180,180]]
		elif b_alignment=='apa':
			return [[0,180],[180,180],[180,0]]
		elif b_alignment=='aap':
			return [[0,180],[180,0],[180,0]]
		elif b_alignment=='app':
			return [[0,180],[180,180],[180,180]]
		elif b_alignment=='paa':
			return [[0,0],[180,0],[180,180]]
		elif b_alignment=='ppa':
			return [[0,0],[180,180],[180,0]]
		elif b_alignment=='pap':
			return [[0,0],[180,0],[180,0]]
		elif b_alignment=='ppp':
			return [[0,0],[180,180],[180,180]]
		else:
			return None

	def get_s_flip_angle(self,s_alignment):
		if s_alignment:
			return [0,0]
		else:
			return [180,180]

	def affine_transformation_sidechain(self,name,angle): # rotate z-axis than x-axis
		angle_y,angle_z = angle
		pymol.cmd.rotate('z',angle_z,name)
		pymol.cmd.rotate('y',angle_y,name)
		# Correct pos
		pos_com = np.mean(get_ca(name)[self.ca_start:self.ca_end+1],0)
		pymol.cmd.translate([0,-pos_com[1],0],name)	

	def affine_transformation_backbone(self,name,angle): # rotate y-axis than z-axis
		angle_y,angle_z = angle
		pymol.cmd.rotate('y',angle_y,name)
		pymol.cmd.rotate('z',angle_z,name)
		# Correct pos
		pos_com = np.mean(get_ca(name)[self.ca_start:self.ca_end+1],0)
		pymol.cmd.translate([0,-pos_com[1],0],name)
		if  (angle_y+angle_z)%360 != 0:
			pos_x_init = get_ca(name)[self.ca_end][0]
		else:
			pos_x_init = get_ca(name)[self.ca_start][0]
		pymol.cmd.translate([-pos_x_init,0,0],name)		

	def show_unit(self):
		pymol.cmd.color('cyan','p_s1_pep2_*')
		pymol.cmd.color('purple','p_s2_pep2_*')
		pymol.cmd.hide('all')
		pymol.cmd.show('sticks','backbone or resn ace+nhe')
		return





