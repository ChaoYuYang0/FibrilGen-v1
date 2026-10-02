import os
abspath = os.path.abspath('')
import sys
sys.path.append(f'{abspath}/../../')

from seq2pep import *

## Create a b-strand structure
# Initialize a b-strand
get_b_strand('FEFKFEFK')
# Align bulky residues on the cross-section
set_chi_angle(1,'n','ca','cb','cg',-20)
set_chi_angle(3,'n','ca','cb','cg',-20)
set_chi_angle(5,'n','ca','cb','cg',-20)
set_chi_angle(7,'n','ca','cb','cg',-20)
# Save structure
pymol.cmd.save('F8.pdb','pep')