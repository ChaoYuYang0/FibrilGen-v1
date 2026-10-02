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
    print('#########################################')
    print('The Energy of %s is %8.2f %s       '  % (name, nrg, ff.GetUnit()))
    print('#########################################')

# def convert_mol2pdb(selection,file):
#     mol_string = pymol.cmd.get_str('mol',selection)
#     obconversion = ob.OBConversion()
#     mol = ob.OBMol()
#     obconversion.ReadString(mol, mol_string)
#     obconversion.WriteFile(mol, file)

def conf_search(selection: str = 'all', forcefield: str = 'MMFF94s',
                method: str = 'Weighted', nsteps1: int = 500,
                conformers: int = 25, lowest_conf: int = 5):
    """
DESCRIPTION

    Perform a conformational search on a molecule using OpenBabel's force fields.

ARGUMENTS

    selection = string: The selection string for the molecule to search.

    forcefield = string: The force field to use (e.g., 'GAFF', 'MMFF94s', 'UFF', 'Ghemical').

    method = string: The search method for global-minimum ('Weighted', 'Random', or 'Systematic').

    nsteps1 = int: The number of optimization steps for each conformer.

    conformers = int: The number of conformers to be analyzed.

    lowest_conf = int: The number of lowest energy conformers to keep.

SEE ALSO

    minimize
    """
    mol_string = pymol.cmd.get_str('mol', selection)
    name = pymol.cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('mol', 'mol')
    mol = ob.OBMol()
    obconversion.ReadString(mol, mol_string)
    ff = ob.OBForceField.FindForceField(forcefield) ## GAFF, MMFF94s, MMFF94, UFF, Ghemical
    ff.Setup(mol)
    if method == 'Weighted':
        ff.WeightedRotorSearch(conformers, nsteps1)
    elif method == 'Random':
        ff.RandomRotorSearch(conformers, nsteps1)
    else:
        ff.SystematicRotorSearch(nsteps1)
    if name == 'all':
        name = 'all_'
    if method in ['Weighted', 'Random']:
        ff.GetConformers(mol)
        print('##############################################')
        print('   Conformer    |         Energy      |  RMSD')
        nrg_unit = ff.GetUnit()
        rmsd = 0
        ff.GetCoordinates(mol)
        nrg = ff.Energy()
        conf_list = []
        for i in range(conformers):
            mol.SetConformer(i) 
            ff.Setup(mol)
            nrg = ff.Energy()
            conf_list.append((nrg, i))
        conf_list.sort()
        lenght_conf_list = len(conf_list)
        if lowest_conf > lenght_conf_list:
            lowest_conf = lenght_conf_list
        for i in range(lowest_conf):
            nrg, orden = conf_list[i]
            name_n = '%s%02d' % (name, i)
            pymol.cmd.delete(name_n)
            mol.SetConformer(orden) 
            mol_string = obconversion.WriteString(mol)
            pymol.cmd.read_molstr(mol_string, name_n,state=0,finish=1,discrete=1)
            if i != 0:
                rmsd = pymol.cmd.fit(name_n, '%s00' % name, quiet=1)
            print('%15s | %10.2f%9s |%6.1f'    % (name_n, nrg, nrg_unit, rmsd))
        print('##############################################')
    else:
        ff.GetCoordinates(mol)
        nrg = ff.Energy()
        mol_string = obconversion.WriteString(mol)
        pymol.cmd.delete(name)
        pymol.cmd.read_molstr(mol_string, name,state=0,finish=1,discrete=1)
        print('#########################################')
        print('The Energy of %s is %8.2f %s       '  % (name, nrg, ff.GetUnit()))
        print('#########################################')

pymol.cmd.extend('minimize', minimize)
pymol.cmd.extend('conf_search', conf_search)
