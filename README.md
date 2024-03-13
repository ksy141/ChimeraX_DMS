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

# Open
```
open a.dms b.dms
```

# Save
```
save a.dms models #1
```

# Notes
Use sort true if you want your atoms reordered. If you use sorting, each residue should have a unique (chain, resid) in DMS.

