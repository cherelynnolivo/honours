from Bio.PDB import PDBParser
import sys
AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
pdb_file = sys.argv[1]
s = PDBParser(QUIET=True).get_structure('p', pdb_file)
for m in s:
    for c in m:
        seq = ''.join([AA.get(r.resname,'X') for r in c if r.id[0]==' '])
        if seq: print(f">Chain_{c.id}\n{seq}")
