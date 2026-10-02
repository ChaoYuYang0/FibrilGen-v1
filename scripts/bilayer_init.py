pymol.cmd.run('pep2unit.py')

## Reset
pymol.cmd.reset()

## Select a reference coordinate
pymol.cmd.select('po1','resi 1 and name ca')
pymol.cmd.select('po2','resi 7 and name ca')
pymol.cmd.select('po3','resi 2 and name c')
pymol.cmd.select('po4','resi 2 and name o')

## Create a periodic unit
unit = create_pep_unit('pep','po1','po2','po3','po4')

## --- Examples of sheet structures ---
sheet = create_sheet(unit,[0,7],[1,7])
# INPUT: ([aaa/apa/aap/app/paa/ppa/pap/ppp,sidechain flip],num of units per sheet)
sheet.build_a_plain_sheet(['aaa','s'],5)
## --- End ---

# Import openbabel and optimization library
pymol.cmd.run('optimize.py')
# Remove steric clashes
pymol.cmd.create('plain_sheet_em','plain_sheet')
minimize('plain_sheet_em',nsteps0=100)

sheet.get_dimension() # Print the refined fibril dimension
pymol.cmd.zoom()