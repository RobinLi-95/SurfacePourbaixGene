import numpy as np
from ase.build import bulk, add_adsorbate, surface, make_supercell

class Mg_10m10_AdsorptionFinder():
    """
    Class finding the typical adsorption site on the Mg (10m10) surface

    Attributes:
            @LayerFinder(): Find the atomic position of atoms in the top 4 layers
            @OnSurfAdsorptionFinder(): Find typical adsorption site (SB,LB,OT,HW) on the top surface
            @SubSurfAdsorptionFinder(bulklatticeA,bulklatticeC): Find typical subsurface adsorption site 
                                                                 (OCTA/TETRA) between the N(1st) and N (2nd)

    """ 
    def __init__(self,slab):
        self.slab = slab

        return
    
    def LayerFinder(self):
        self.toplayer = self.slab.get_positions()[self.slab.get_positions()[:,2].argsort()][-4:]
        self.sortedtl = self.toplayer[self.toplayer[:,0].argsort()]

        self.secondlayer = self.slab.get_positions()[self.slab.get_positions()[:,2].argsort()][-8:-4]
        self.sortedsl = self.secondlayer[self.secondlayer[:,0].argsort()]

        self.thirdlayer = self.slab.get_positions()[self.slab.get_positions()[:,2].argsort()][-12:-8]
        self.sortedthl = self.thirdlayer[self.thirdlayer[:,0].argsort()]

        self.fourthlayer = self.slab.get_positions()[self.slab.get_positions()[:,2].argsort()][-16:-12]
        self.sortedfl = self.fourthlayer[self.fourthlayer[:,0].argsort()]

        return


    def OnSurfAdsorptionFinder(self):
        self.dx = self.sortedtl[2][0] - self.sortedtl[1][0]
        self.dy = self.sortedtl[1][1] - self.sortedtl[0][1]

        # determine the initial adsorption site, with offset = [0,0]
        self.SB_init = (self.sortedtl[0] + self.sortedtl[3]) / 2
        self.LB_init = (self.sortedtl[0] + self.sortedtl[1]) / 2
        self.OT_init = self.sortedtl[0] + np.array([0,0,2]) # with a z-distance with 2 Angstrom
        self.HW_init = (self.sortedtl[0] + self.sortedtl[2]) / 2 + np.array([0,0,1.5])

        # find all periodic adsorption sites within the unit cell, set as (2 x 2)
        # order of the offset is (1,0),(1,1),(0,1),(0,0)
        self.SB = np.array([self.SB_init, 
                            self.SB_init + np.array([self.dx,0,0]),
                            self.SB_init + np.array([self.dx,self.dy,0]),
                            self.SB_init + np.array([0,self.dy,0])
                            ])
        
        self.LB = np.array([self.LB_init, 
                            self.LB_init + np.array([self.dx,0,0]),
                            self.LB_init + np.array([self.dx,self.dy,0]),
                            self.LB_init + np.array([0,self.dy,0])
                            ])

        self.OT = np.array([self.OT_init, 
                            self.OT_init + np.array([self.dx,0,0]),
                            self.OT_init + np.array([self.dx,self.dy,0]),
                            self.OT_init + np.array([0,self.dy,0])
                            ])
        
        self.HW = np.array([self.HW_init, 
                            self.HW_init + np.array([self.dx,0,0]),
                            self.HW_init + np.array([self.dx,self.dy,0]),
                            self.HW_init + np.array([0,self.dy,0])
                            ])

        return

    def SubSurfAdsorptionFinder(self,bulklatticeA,bulklatticeC):
        self.tetradist = bulklatticeC/2 * (1/2 + 2/3 * ((bulklatticeA/bulklatticeC)**2 )) 

        # Determine Octahedral(A/B) Sites
        self.OCTa_init = (self.sortedtl[0] + self.sortedfl[1]) / 2

        self.OCTa = np.array([self.OCTa_init, 
                              self.OCTa_init + np.array([self.dx,0,0]),
                              self.OCTa_init + np.array([self.dx,self.dy,0]),
                              self.OCTa_init + np.array([0,self.dy,0])
                             ])

        self.OCTb_init = (self.sortedtl[0] + self.sortedfl[0]) / 2

        self.OCTb = np.array([self.OCTb_init, 
                              self.OCTb_init + np.array([self.dx,0,0]),
                              self.OCTb_init + np.array([self.dx,self.dy,0]),
                              self.OCTb_init + np.array([0,self.dy,0])
                             ])

        # Determine Tetrahedral-I (A/B) Sites
        self.TetraIa_init = self.sortedsl[1] + np.array([0,self.tetradist,0])

        self.TetraIa = np.array([self.TetraIa_init, 
                                 self.TetraIa_init + np.array([self.dx,0,0]),
                                 self.TetraIa_init + np.array([self.dx,self.dy,0]),
                                 self.TetraIa_init + np.array([0,self.dy,0])
                                ])

        self.TetraIb_init = self.sortedsl[0] - np.array([0,self.tetradist,0])

        self.TetraIb = np.array([self.TetraIb_init, 
                                 self.TetraIb_init + np.array([self.dx,0,0]),
                                 self.TetraIb_init + np.array([self.dx,self.dy,0]),
                                 self.TetraIb_init + np.array([0,self.dy,0])
                                ])

        # Determine Tetrahedral-II (A/B) Sites
        self.TetraIIa_init = self.sortedthl[0] - np.array([0,self.tetradist,0])

        self.TetraIIa = np.array([self.TetraIIa_init, 
                                  self.TetraIIa_init + np.array([self.dx,0,0]),
                                  self.TetraIIa_init + np.array([self.dx,self.dy,0]),
                                  self.TetraIIa_init + np.array([0,self.dy,0])
                                 ])

        self.TetraIIb_init = self.sortedthl[0] + np.array([0,self.tetradist,0])

        self.TetraIIb = np.array([self.TetraIIb_init, 
                                  self.TetraIIb_init + np.array([self.dx,0,0]),
                                  self.TetraIIb_init + np.array([self.dx,self.dy,0]),
                                  self.TetraIIb_init + np.array([0,self.dy,0])
                                 ])

        return