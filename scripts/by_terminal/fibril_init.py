import os
abspath = os.path.abspath('')
import sys
sys.path.append(f'{abspath}/../../')

from builder import *

## Load PDB
pymol.cmd.load('bilayer_em_indexed2.pdb')

## Create a periodic unit
# INPUT: (peptide length, the number of peptides per sheet, start residue index, end residue index)
unit = create_sheet_unit(8,10,2,8)

## Create a fibril object
fibril = create_fibril(unit)

## --- Examples of structures ---
## Build a plain sheet (num of units per sheet)
# fibril.build_a_flat_sheet(10)
## Build a stacked sheet (stacking pattern, num of units per sheet)
# fibril.build_a_stacked_sheet([[0,1],[1,1]],10)								
## Build a rod (tilt angle, num of units per sheet, twist sign)
# fibril.build_a_rod(20,50,1) 	
## Build a stacked rod (tilt angle, stacking pattern, num of units per sheet, twist sign)								
# fibril.build_a_stacked_rod(10,[[0,1],[1,1]],10,1)	
## Build a ribbon (tilt angle, radius, num of units per sheet, twist sign)
# fibril.build_a_ribbon(30,60,10,-1)
## Build a stacked ribbon (tilt angle, radius, stacking angle, stacking number, num of units per sheet, twist sign)
fibril.build_a_stacked_ribbon(10,30,70,3,20,1)
## --- End ---

## Save structures
if fibril.write_output == 1:
	# pymol.cmd.save('a_sheet.pdb','plain_sheet')
	# pymol.cmd.save('s_sheet.pdb','plain_sheet')
	# pymol.cmd.save('a_rod.pdb','a_rod')
	# pymol.cmd.save('s_rod.pdb','s_rod')
	# pymol.cmd.save('a_ribbon.pdb','a_ribbon')
	pymol.cmd.save('s_ribbon.pdb','s_ribbon')
