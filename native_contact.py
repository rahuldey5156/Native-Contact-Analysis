import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF,DCD


u = mda.Universe('protein_0ns.pdb','try1.xtc')

sel_basic = "(segid D) and not (name H* or name 1H* or name 2H* or name 3H*)"
sel_acidic = "(segid E) and not (name H* or name 1H* or name 2H* or name 3H*)"

acidic = u.select_atoms(sel_acidic)
basic = u.select_atoms(sel_basic)


ca3 = contacts.Contacts(u, select=(sel_acidic, sel_basic),
                        refgroup=(acidic, basic),
                        radius=5.0,
                        method='soft_cut',
                        kwargs={'beta': 5.0,
                                'lambda_constant': 1.8}).run()

ca3_df = pd.DataFrame(ca3.timeseries,
                      columns=['Frame',
                                'Contacts from first frame'])


ca3_df.to_csv('D_E_v1.csv', sep='\t', index=False)

