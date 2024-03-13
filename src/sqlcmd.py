# in nonbonded, param is the atom type in integer

sql_create = """

CREATE TABLE global_cell (id integer primary key, x float, y float, z float);

CREATE TABLE particle (
  id integer primary key,
  anum integer,
  name text not null,
  x float,
  y float,
  z float,
  vx float,
  vy float,
  vz float,
  resname text not null,
  resid integer,
  chain text not null,
  segname text not null,
  mass float,
  charge float,
  formal_charge integer,
  insertion text not null,
  msys_ct integer not null,
  'm_grow_name' text,
  'm_mmod_type' integer,
  nbtype integer not null,
  type text not null,
  bfactor float
);

CREATE TABLE bond (p0 integer, p1 integer, 'order' integer);
"""

sql_insert_particle = """
INSERT INTO particle
    (id, anum, name, resname, chain, resid, mass, charge, x, y, z, vx, vy, vz, segname, insertion, msys_ct, nbtype, type, bfactor)
    VALUES
    (?,  ?,    ?,    ?,       ?,     ?,     ?,    ?,      ?, ?, ?, ?,  ?,  ?,  ?,       ?,         ?,       ?,      ?,    ?);
"""

sql_insert_bond = "INSERT INTO bond      VALUES('{:d}', '{:d}', '{:d}');"

