import numpy as np
import copy as cp
from ase.build import bulk, add_adsorbate, surface, make_supercell

class Mg_10m10_AdsorptionAdder():
    """
    Class adding the O/OH/H on the Mg (10m10) surface

    Attributes:
            @SymmetricSiteFinder(): Find the corresponding symmetric site of the input adsorption site
            @PureAdsorption(): Adsorbing "Single Type" radicals/atoms on/in the Mg (10m10) surface
            @HybridAdsorption(): Adsorbing multiple types of radicals/atoms on the Mg (10m10) surface

    """ 
    def __init__(self,slab,site):
        '''
        @slab: input the slab generated by the ASE modules
        @site: input the position of the adsorption sites
        '''
        self.slab = cp.deepcopy(slab)
        self.site = site
        self.toppose = self.slab.get_positions()[self.slab.get_positions()[:,2].argsort()][-1][2]
        self.zoffset = self.site[0][2] - self.toppose
        self.OHbond = 0.964

        return

    def SymmetricSiteFinder(self,adsite):
        self.dx = self.slab.get_cell()[0][0] / 2
        self.dy = self.slab.get_cell()[1][1] / 2
        self.Z = self.slab.get_cell()[2][2]

        self.movingVector = np.array([0,self.dy/2,self.Z-2*adsite[2]])
        self.mirrorsite = adsite + self.movingVector
        self.mirrorzoffset = self.mirrorsite[2] - self.toppose

        return

    def PureAdsorption(self,coverage,adsorbate):
        # need consider OH adsorption and symmetric adsorption
        self.coverage = int(coverage * 4) # coverage is the fractional coverage, 0.25, 0.50, 0.75, 1.0.
        self.adsorbate = adsorbate

        if self.coverage !=2:
            slab_tmp = cp.deepcopy(self.slab)
            name = 'slab_Pure'

            for i in range(self.coverage):
                
                add_adsorbate(slab_tmp, self.adsorbate[0], self.zoffset, (self.site[i][0],self.site[i][1]))
                self.SymmetricSiteFinder(self.site[i])
                add_adsorbate(slab_tmp, self.adsorbate[0], self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))

                if self.adsorbate == 'OH':
                    add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

            setattr(self,name,slab_tmp)

        # for coverage = 0.5 ML(2 X 2 slab), there are 3 different adsorption pattern
        elif self.coverage == 2:
            
            for i in range(1,4):
                name = 'slab_Pure2_p' + str(i)
                slab_tmp = cp.deepcopy(self.slab)
                
                add_adsorbate(slab_tmp, self.adsorbate[0], self.zoffset, (self.site[0][0],self.site[0][1]))
                self.SymmetricSiteFinder(self.site[0])
                add_adsorbate(slab_tmp, self.adsorbate[0], self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))
                
                add_adsorbate(slab_tmp, self.adsorbate[0], self.zoffset, (self.site[i][0],self.site[i][1]))
                self.SymmetricSiteFinder(self.site[i])
                add_adsorbate(slab_tmp, self.adsorbate[0], self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))
                
                if self.adsorbate == 'OH':
                    add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[0][0],self.site[0][1]))
                    self.SymmetricSiteFinder(self.site[0])
                    add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))
                    
                    add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    self.SymmetricSiteFinder(self.site[i])
                    add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                setattr(self,name,slab_tmp)
        return


    def HybridAdsorption(self,N_o,N_oh):
        self.totN = N_o + N_oh
        self.N_o = N_o
        self.N_oh = N_oh

        if self.totN == 2:
            # first adsorb O
            self.slab_1O_1OH = cp.deepcopy(self.slab)
            add_adsorbate(self.slab_1O_1OH, 'O', self.zoffset, (self.site[0][0],self.site[0][1]))
            self.SymmetricSiteFinder(self.site[0])
            add_adsorbate(self.slab_1O_1OH, 'O', self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))
            
            # then adsorb OH 
            for i in range(1,4):
                name = 'slab_1OH_1O_p' + str(i)
                slab_tmp = cp.deepcopy(self.slab_1O_1OH)

                add_adsorbate(slab_tmp, 'O', self.zoffset, (self.site[i][0],self.site[i][1]))
                self.SymmetricSiteFinder(self.site[i])
                add_adsorbate(slab_tmp, 'O', self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))

                add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                setattr(self,name,slab_tmp)


        elif self.totN == 3:
            # No matter how O/OH adsorption patterns are, the O could be adsorbed firstly.
            self.slab_O_OH_c3 = cp.deepcopy(self.slab)

            for i in range(3):
                add_adsorbate(self.slab_O_OH_c3, 'O', self.zoffset, (self.site[i][0],self.site[i][1]))
                self.SymmetricSiteFinder(self.site[i])
                add_adsorbate(self.slab_O_OH_c3, 'O', self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))

            if self.N_o == 2:
                for i in range(3):
                    name = 'slab_1OH_2O_p' + str(i+1)
                    slab_tmp = cp.deepcopy(self.slab_O_OH_c3)

                    add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    self.SymmetricSiteFinder(self.site[i])
                    add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))
                    
                    setattr(self,name,slab_tmp)

            if self.N_o == 1:
                index = 1
                for i in range(2):
                   
                    slab_tmp1 = cp.deepcopy(self.slab_O_OH_c3)

                    add_adsorbate(slab_tmp1, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    self.SymmetricSiteFinder(self.site[i])
                    add_adsorbate(slab_tmp1, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                    for j in range(i+1,3):
                        name = 'slab_2OH_1O_p' + str(index)

                        slab_tmp2 = cp.deepcopy(slab_tmp1)

                        add_adsorbate(slab_tmp2, 'H', self.zoffset+self.OHbond, (self.site[j][0],self.site[j][1]))
                        self.SymmetricSiteFinder(self.site[j])
                        add_adsorbate(slab_tmp2, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                        setattr(self,name,slab_tmp2)
                        index +=1
                        

        elif self.totN == 4:
            # No matter how O/OH adsorption patterns are, the O could be adsorbed firstly.
            self.slab_O_OH_c4 = cp.deepcopy(self.slab)

            for i in range(4):
                add_adsorbate(self.slab_O_OH_c4, 'O', self.zoffset, (self.site[i][0],self.site[i][1]))
                self.SymmetricSiteFinder(self.site[i])
                add_adsorbate(self.slab_O_OH_c4, 'O', self.mirrorzoffset, (self.mirrorsite[0],self.mirrorsite[1]))

            if self.N_oh == 1:
                name = 'slab_1OH_3O_p1'

                slab_tmp = cp.deepcopy(self.slab_O_OH_c4)

                add_adsorbate(slab_tmp, 'H', self.zoffset+self.OHbond, (self.site[0][0],self.site[0][1]))
                self.SymmetricSiteFinder(self.site[0])
                add_adsorbate(slab_tmp, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                setattr(self,name,slab_tmp)

            if self.N_oh == 2:
                slab_tmp1 = cp.deepcopy(self.slab_O_OH_c4)

                add_adsorbate(slab_tmp1, 'H', self.zoffset+self.OHbond, (self.site[0][0],self.site[0][1]))
                self.SymmetricSiteFinder(self.site[0])
                add_adsorbate(slab_tmp1, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                for i in range(1,4):
                    name = 'slab_2OH_2O_p' + str(i)
                    slab_tmp2 = cp.deepcopy(slab_tmp1)
                    add_adsorbate(slab_tmp2, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    self.SymmetricSiteFinder(self.site[i])
                    add_adsorbate(slab_tmp2, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                    setattr(self,name,slab_tmp2)

            if self.N_oh == 3:
                name = 'slab_3OH_1O_p1'

                slab_tmp1 = cp.deepcopy(self.slab_O_OH_c4)

                for i in range(3):
                    add_adsorbate(slab_tmp1, 'H', self.zoffset+self.OHbond, (self.site[i][0],self.site[i][1]))
                    self.SymmetricSiteFinder(self.site[i])
                    add_adsorbate(slab_tmp1, 'H', self.mirrorzoffset-self.OHbond, (self.mirrorsite[0],self.mirrorsite[1]))

                setattr(self,name,slab_tmp1)

        return
