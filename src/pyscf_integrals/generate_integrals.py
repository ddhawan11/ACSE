import numpy as np
import pyscf
import os
from pyscf import gto, ao2mo

molecule = "H6_"
mol = gto.Mole()
R = 0.75
label = molecule + str(R)
atom = [
        ["h",   (0,       0.0,        0.0)],   
        ["h",   ( R,       0.0,        0.0)],
        ["h",   ( 2*R,      0.0,        0.0)],   
        ["h",   ( 3*R,       0.0,        0.0)], 
        ["h",   ( 4*R,        0.0,      0.0)],   
        ["h",   ( 5*R,        0.0,       0.0)],  
        ]

basis = "sto-3g"
mol.atom = atom
mol.build()
mf = pyscf.scf.RHF(mol).run(conv_tol=1e-8)
h0 = mf.energy_nuc()
print(mf.mo_energy)
C = mf.mo_coeff    
norb = C.shape[1]
h1 = C.T @ mf.get_hcore(mol) @ C
h2 = pyscf.ao2mo.kernel(mol, C, aosym="s4",compact=False)
h2.shape = (norb, norb, norb, norb)

os.mkdir(label)
np.save(label+'/h0.npy', h0)
np.save(label+'/h1.npy', h1)
np.save(label+'/h2.npy', h2)
