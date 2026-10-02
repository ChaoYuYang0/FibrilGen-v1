pymol.cmd.run('seq2pep.py')

## Reset
pymol.cmd.delete('all')
pymol.cmd.reset()

## Create a b-strand structure
# Initialize a b-strand
get_b_strand('FEFKFEFK')
# Align bulky residues on the cross-section
set_chi_angle(1,'n','ca','cb','cg',-20)
set_chi_angle(3,'n','ca','cb','cg',-20)
set_chi_angle(5,'n','ca','cb','cg',-20)
set_chi_angle(7,'n','ca','cb','cg',-20)