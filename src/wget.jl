"""

```
wget(PDBid; selection)
```
Retrieves a PDB file from the protein data bank. Selections may be applied.

### Example

```julia-repl
julia> protein = wget("1LBD","chain A")
   Array{Atoms,1} with 1870 atoms with fields:
   index name resname chain   resnum  residue        x        y        z  beta occup model segname index_pdb
       1    N     SER     A      225        1   45.228   84.358   70.638 67.05  1.00     1       -         1
       2   CA     SER     A      225        1   46.080   83.165   70.327 68.73  1.00     1       -         2
       3    C     SER     A      225        1   45.257   81.872   70.236 67.90  1.00     1       -         3
                                                       ⋮ 
    1868  OG1     THR     A      462      238  -27.462   74.325   48.885 79.98  1.00     1       -      1868
    1869  CG2     THR     A      462      238  -27.063   71.965   49.222 78.62  1.00     1       -      1869
    1870  OXT     THR     A      462      238  -25.379   71.816   51.613 84.35  1.00     1       -      1870

```
"""
function wget(pdb_id::String, selection::String )  
  query = parse_query(selection)
  return wget(pdb_id, only = atom -> apply_query(query,atom) )
end

function wget(pdb_id::String; only = all)
  file = download("https://files.rcsb.org/download/$(pdb_id).pdb")
  return readPDB(file,only=only)
end

