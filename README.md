# ChimeraX_DMS
Reading and writing Desmond Molecular Structure (DMS) format files in ChimeraX

# Install
In ChimeraX,
```
devel build   ~/ChimeraX_DMS
devel install ~/ChimeraX_DMS
```

# Remove
```
toolshed uninstall ChimeraX_DMS
```

# Notes
1. Each residue should have a unique (chain, resid) in DMS.
2. When opening a DMS file, atoms are reordered (!)
