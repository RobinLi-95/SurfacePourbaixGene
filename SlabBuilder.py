
import numpy as np
from ase.build import bulk, add_adsorbate, surface, make_supercell


class SlabBuilder():
    """
    Class constructing the bulk crystal structure and correlated slab model

    Attributes:
            @BulkStructureBuilder(): Construct the bulk structure of the nominated element
            @SlabBuilder(list,int,float,list): Construct the slab
            @StructureWriter(str,str): Output the constructed structure with a specific format
    """    

    def __init__(self,atom,lattice,lattice_parameter):
        self.atomtype = atom # atom type

        self.lattice = lattice # lattice type: fcc, hcp, bcc, ...
        self.latticeA = lattice_parameter[0]
        self.latticeB = lattice_parameter[1]
        self.latticeC = lattice_parameter[2]

    def BulkStructureBuilder(self):
        self.bulkstructure =  bulk(self.atomtype, self.lattice, self.latticeA, self.latticeB, self.latticeC)       
        return
    
    def SlabBuilder(self,surfaceplane,Nlayer,Vacuum,supercell):
        self.surfaceplane = surfaceplane
        self.Nlayer = int(Nlayer)
        self.Vacuum = float(Vacuum)
        self.supercell = np.array([[supercell[0],0,0],[0,supercell[1],0],[0,0,supercell[2]]])

        self.slab = surface(self.bulkstructure, (surfaceplane[0],surfaceplane[1],surfaceplane[2]), layers = self.Nlayer, vacuum = self.Vacuum)
        self.supercellslab = make_supercell(self.slab,self.supercell)

        return