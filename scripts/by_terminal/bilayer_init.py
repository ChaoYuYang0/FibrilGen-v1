import os
abspath = os.path.abspath('')
import sys
sys.path.append(f'{abspath}/../../')

from pep2unit import *

## Load PDB
pymol.cmd.load(f'F8.pdb')
## Select a reference coordinate
pymol.cmd.select('po1','resi 1 and name ca')
pymol.cmd.select('po2','resi 7 and name ca')
pymol.cmd.select('po3','resi 2 and name c')
pymol.cmd.select('po4','resi 2 and name o')

## Create a periodic unit
unit = create_pep_unit('F8','po1','po2','po3','po4')
# unit.rotate_sidechain() ## list of angles

## --- Examples of sheet structures ---
sheet = create_sheet(unit,[0,7],[1,7])
# INPUT: ([aaa/apa/aap/app/paa/ppa/pap/ppp,sidechain flip],num of units per sheet)
sheet.build_a_plain_sheet(['aaa','s'],5)
## --- End ---

## Save the structure
pymol.cmd.save('bilayer.pdb','plain_sheet')

# Import openbabel and optimization library
from optimize import *
# Remove steric clashes
pymol.cmd.create('plain_sheet_em','plain_sheet')
minimize('plain_sheet_em',nsteps0=100)
pymol.cmd.save('bilayer_em.mol','plain_sheet_em')
