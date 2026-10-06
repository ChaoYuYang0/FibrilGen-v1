import numpy as np
import pymol
import math

def padding_white_spaces(length_str,input_str):
	out = ''
	for _ in range(length_str-len(input_str)):
		out += ' '
	return out+input_str

def ordering_indices(root,file_in,file_out):
	# Open files
	f_in = open(root+file_in,'r')
	f_out = open(root+file_out,'w')
	# Change atom and residue indices
	lines = f_in.readlines()
	current_res,acc_atom,acc_res = '',0,0
	for line in lines:
		if line[:4] == 'ATOM':
			idx_atom,idx_res,name_atom,name_res = line[7:11],line[23:26],line[13:16],line[17:20]
			# print (name_atom)
			if ((name_res != current_res) or (name_atom == 'N  ')):
				current_res = name_res
				acc_res += 1
			acc_atom += 1
			out = line[:7]+padding_white_spaces(4,str(acc_atom))+line[11:23]+padding_white_spaces(3,str(acc_res))+line[26:]
		else:
			out = line
		f_out.write(out)
	f_out.write('TER \n')
	f_out.write('END \n')

	f_in.close()
	f_out.close()
	return

def reordering_indices(root,filename_ref,filename_in,filename_out,num_pep,num_res):
	# Open files
	pymol.cmd.load(root+filename_ref+'.pdb')
	pymol.cmd.load(root+filename_in+'.pdb')
	# Get coordinates
	pos_ref = pymol.cmd.get_coords(filename_ref,1).reshape((num_pep,-1,3))
	pos_in = pymol.cmd.get_coords(filename_in,1).reshape((num_pep,-1,3))
	# Compare indices
	shift = num_pep*num_res
	for (idx_in_i,pos_in_i) in enumerate(pos_in):
		sub = pos_ref-pos_in_i.reshape((1,-1,3))
		dist = np.sum(np.linalg.norm(sub,axis=2),axis=1)
		idx_out_i = np.argmin(dist)
		shift_i = (idx_out_i-idx_in_i)*num_res+shift
		pymol.cmd.select(f'pep_in_{idx_in_i}',f'resi {int(idx_in_i*num_res+1)}-{int((idx_in_i+1)*num_res)}')
		pymol.cmd.alter(f'pep_in_{idx_in_i}',f'resi=str(int(resi)+{int(shift_i)})')
	pymol.cmd.alter(f'pep_in_*',f'resi=str(int(resi)-{int(shift)})')
	pymol.cmd.delete('*pep_in_*')
	pymol.cmd.save(f'{root}{filename_out}.pdb',f'{filename_in}')
	return

