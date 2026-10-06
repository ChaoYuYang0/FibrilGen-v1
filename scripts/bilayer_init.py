pymol.cmd.run('pep2unit.py')

## Reset
pymol.cmd.reset()

## Create a periodic unit
# INPUT: (peptide, start residue index, end residue index, does the start residue locate at the face)
unit = create_pep_unit('pep',2,8,0)

## --- Examples of sheet structures ---
# INPUT: (peptide unit)
sheet = create_sheet(unit)
# INPUT: (peptide alignment aaa/apa/aap/app/paa/ppa/pap/ppp, does two beta-sheet aligned face-to-face, num of units per sheet)
sheet.build_a_plain_sheet('ppa',1,5)
## --- End ---

# Import openbabel and optimization library
pymol.cmd.run('optimize.py')
# Remove steric clashes
pymol.cmd.create('plain_sheet_em','plain_sheet')
minimize('plain_sheet_em',nsteps0=100)

sheet.get_dimension() # Print the refined fibril dimension
pymol.cmd.zoom()