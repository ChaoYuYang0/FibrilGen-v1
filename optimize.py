'''
Molcular optimization Python library

PyMOL wiki: http://www.pymolwiki.org/index.php/optimize
Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: august 2014
License: MIT License
Version 0.9
'''

import sys
import pymol

try:
    from openbabel import openbabel as ob
except ImportError:
    print('<' * 80 + '''

Optimize plug-in needs openbabel to be installed in your system, please follow the instructions at
http://openbabel.org/wiki/Get_Open_Babel

''' + '>' * 80)

def minimize(selection: str ='all', forcefield: str ='MMFF94s',
             method: str ='Conjugate Gradients', nsteps0: int = 500,
             conv: float = 0.0001, cutoff: bool = False,
             cut_vdw: float = 6.0, cut_elec: float = 8.0) -> None:
    """
DESCRIPTION

    Minimize the energy of a molecule using OpenBabel's force fields.

ARGUMENTS

    selection = string: The selection string for the molecule to minimize.

    forcefield = string: The force field to use (e.g., 'GAFF', 'MMFF94s', 'UFF', 'Ghemical').

    method = string: The optimization method ('Conjugate Gradients' or 'Steepest Descent').

    nsteps0 = int: The number of optimization steps.

    conv = float: The convergence criterion.

    cutoff = bool: Whether to use cutoff for van der Waals and electrostatic interactions.

    cut_vdw = float: The cutoff distance for van der Waals interactions.

    cut_elec = float: The cutoff distance for electrostatic interactions.

SEE ALSO

    conf_search

    """
    mol_string = pymol.cmd.get_str('mol',selection)
    name = pymol.cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('mol', 'mol')
    mol = ob.OBMol()
    obconversion.ReadString(mol, mol_string)
    ff = ob.OBForceField.FindForceField(forcefield) ## GAFF, MMFF94s, MMFF94, UFF, Ghemical
    ff.Setup(mol)
    if cutoff == True:
        ff.EnableCutOff(True)
        ff.SetVDWCutOff(cut_vdw)
        ff.SetElectrostaticCutOff(cut_elec)
    if method == 'Conjugate Gradients':
        ff.ConjugateGradients(nsteps0, conv)
    else:
        ff.SteepestDescent(nsteps0, conv)
    ff.GetCoordinates(mol)
    nrg = ff.Energy()
    mol_string = obconversion.WriteString(mol)
    pymol.cmd.delete(name)
    if name == 'all':
        name = 'all_'
    pymol.cmd.read_molstr(mol_string, name,state=0,finish=1,discrete=1)

pymol.cmd.extend('minimize', minimize)
