# vim: set expandtab shiftwidth=4 softtabstop=4:

# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/atomic/__init__.py
# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/atomic/struct_edit.py

import sqlite3
import numpy as np
import os
from   chimerax.atomic import AtomicStructure
from   chimerax.atomic.struct_edit import add_atom, add_bond
from   .sqlcmd   import *

from chimerax.atomic import all_atoms
# atoms = all_atoms(session)
# atoms = session.models.list(type=AtomicStructure)[0].atoms
# bonds = session.models.list(type=AtomicStructure)[0].bonds

def open_dms(session, file_name):
    """Read DMS
    Returns the 2-tuple return value expected by the
    "open command" manager's :py:meth:`run_provider` method.
    """
    struct = AtomicStructure(session)

    conn   = sqlite3.connect(file_name)
    atoms  = conn.execute(("SELECT name, anum, resname, resid, chain, "
                           "x, y, z, id " 
                           "FROM particle ORDER BY CHAIN, RESID;")).fetchall()
    #atoms  = conn.execute('SELECT name, anum, resname, resid, chain, x, y, z FROM particle;').fetchall()
    bonds  = conn.execute('SELECT * FROM bond;').fetchall()
    

    ### Add atoms
    n_chains   = 0
    n_residues = 0
    chain_prev = None
    resi_prev  = None
    old2new    = {}

    for i, atom in enumerate(atoms):
        name  = atom[0]
        anum  = atom[1]
        resn  = atom[2]
        resi  = atom[3]
        chain = atom[4]
        posx  = atom[5]
        posy  = atom[6]
        posz  = atom[7]
        index = atom[8]
        old2new[index] = i
        
        if chain_prev == chain and resi_prev == resi:
            # same residue
            pass
        else:
            if chain_prev != chain: n_chains += 1
            residue     = struct.new_residue(resn, chain, resi)
            chain_prev  = chain
            resi_prev   = resi
            n_residues += 1
        
        a = add_atom(name, anum, residue, np.array([posx, posy, posz]))


    ### Add bonds (if written in DMS)
    new_atoms = struct.atoms
    for bond in bonds:
        b0 = old2new[bond[0]]
        b1 = old2new[bond[1]]
        add_bond(new_atoms[b0], new_atoms[b1])
    

    ### Make bonds based on distance if bonds are not in DMS
    if len(bonds) == 0:
        struct.connect_structure()
    else:
        pass
    
    status = (f"Opened {file_name} containing "
              f"{n_chains} chains, "
              f"{n_residues} residues, "
              f"{len(atoms)} atoms, "
              f"{len(bonds)} bonds")

    return [struct], (status)


def save_dms(session, path, models=None):

    if os.path.exists(path): os.remove(path)

    conn = sqlite3.connect(path)
    conn.executescript(sql_create)

    msys_ct = 0
    vx=0.0; vy=0.0; vz=0.0; insertion=''; atype='TBD'; segname=''
    charge = 0.0; formal_charge = 0; nbtype=0; mass=0.0

    if models is None:
        models = session.models.list(type=AtomicStructure)

    num_atoms = 0
    num_bonds = 0
    
    assert len(models) == 1, ('Currently only supports saving one model. '
                              'Specify model (e.g. save xxx.dms models #3)')

    for s in models:
        # We get the list of atoms and transformed atomic coordinates
        # as arrays so that we can limit the number of accesses to
        # molecular data, which is slower than accessing arrays directly
        atoms  = s.atoms
        bonds  = s.bonds
        coords = atoms.scene_coords
        index_match = {}

        for i in range(len(atoms)):
            atom = atoms[i]
            pos  = coords[i]

            # if you read pdb -> index starts from 1
            # atom.serial_number becomes discontinuous at TER
            index_match[atom.serial_number] = i
            elem  = atom.element.number
            name  = atom.name
            #mass  = atom.element.mass
            chain = atom.residue.chain_id
            resn  = atom.residue.name
            resi  = atom.residue.number

            conn.execute(sql_insert_particle, (
                i, elem, name, resn,
                chain, resi, mass, charge,
                pos[0], pos[1], pos[2], vx, vy, vz, segname,
                insertion, msys_ct, nbtype, atype, atom.bfactor))
            
            num_atoms += 1


        for bond in bonds:
            i0 = index_match[bond.atoms[0].serial_number]
            i1 = index_match[bond.atoms[1].serial_number]
            conn.execute(sql_insert_bond.format(i0, i1, 1))
            num_bonds += 1

    conn.commit()
    conn.close()

    status = (f"Saved {path} containing "
              f"{num_atoms} atoms, "
              f"{num_bonds} bonds")

    session.logger.status(status)


