# vim: set expandtab shiftwidth=4 softtabstop=4:
# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/pdb/__init__.py
# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/pdb/pdb.py
# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/atomic/__init__.py
# /Applications/ChimeraX-1.5.app/Contents/Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/chimerax/atomic/struct_edit.py

import sqlite3
import numpy as np
import os
#from   chimerax.atomic import AtomicStructure
from   chimerax.atomic.structure import AtomicStructure
from   chimerax.atomic.structure import Structure
from   chimerax.atomic.struct_edit import add_atom, add_bond
from   .sqlcmd   import *

from chimerax.atomic import all_atoms
# atoms = all_atoms(session)
# atoms = session.models.list(type=AtomicStructure)[0].atoms
# bonds = session.models.list(type=AtomicStructure)[0].bonds

def open_dms(session, path, file_name, *, atomic=True, sort=False, connect=False):
    """Read DMS
    Returns the 2-tuple return value expected by the
    "open command" manager's :py:meth:`run_provider` method.
    """
    name = '.'.join(os.path.basename(path).split('.')[:-1])

    if atomic:
        struct = AtomicStructure(session, name=name)
    else:
        struct = Structure(session, name=name)

    conn  = sqlite3.connect(path)
    cols  = [row[1] for row in conn.execute("PRAGMA table_info(particle);").fetchall()]  # Extract column names
    if 'id' not in cols: raise IOError("id is not availabe in DMS")

    data  = {'id': [], 'name': 'TBD', 'anum': 6, 'resname': 'TBD', 'resid': 1, 'chain': 'X', 'x': 0.0, 'y': 0.0, 'z': 0.0, 'charge': 0.0}
    for key in cols: data[key] = []
    
    avail = [key for key in data.keys() if isinstance(data[key], list)]
    if sort:
        query   = "SELECT " + " , ".join(avail) + " FROM particle ORDER BY CHAIN, RESID;"
    else:
        query   = "SELECT " + " , ".join(avail) + " FROM particle;"
    
    print(query)
    atoms = conn.execute(query).fetchall()
    bonds = conn.execute('SELECT * FROM bond;').fetchall()
    
    ### Add atoms
    for n_atoms, atom in enumerate(atoms, 1):
        for j in range(len(avail)):
            data[avail[j]].append(atom[j])
    
    ### Check columns that did not exist in structure
    for key, value in data.items():
        if isinstance(value, list):
            data[key] = np.array(value)
        else:
            data[key] = np.array([value] * n_atoms)
    data['newid'] = np.arange(0, n_atoms)

    ### Count residues
    bA = (data['resid'][1:] == data['resid'][:-1]) & (data['chain'][1:] == data['chain'][:-1])
    data['resn'] = np.insert(np.cumsum(~bA), 0, 0)

    ### Make structure
    for i in range(n_atoms):
        if i == 0:
            residue = struct.new_residue(data['resname'][i], data['chain'][i], data['resid'][i])
        elif data['resn'][i] != data['resn'][i-1]:
            residue = struct.new_residue(data['resname'][i], data['chain'][i], data['resid'][i])
        a = add_atom(data['name'][i], int(data['anum'][i]), residue, np.array([data['x'][i], data['y'][i], data['z'][i]]), bfactor=data['charge'][i])

    ### Add bonds (if written in DMS)
    new_atoms = struct.atoms
    for bond in bonds:
        b0 = data['newid'][np.where(data['id'] == bond[0])[0][0]]
        b1 = data['newid'][np.where(data['id'] == bond[1])[0][0]]
        if b0 < b1:
            add_bond(new_atoms[b0], new_atoms[b1])
        else:
            add_bond(new_atoms[b1], new_atoms[b0])

    ### Make bonds based on distance if bonds are not in DMS
    #if len(bonds) == 0:
    if connect:
        struct.connect_structure()
    
    status = (f"Opened {path} containing "
              f"{len(set(data['chain']))} chains, "
              f"{len(set(data['resn']))} residues, "
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


