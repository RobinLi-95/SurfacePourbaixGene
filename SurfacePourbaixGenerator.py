import numpy as np

class SurfacePourbaixDiagramGenerator():
    """
    Class for generating 1-D (fixed pH or potential) and 2-D Surface Pourbaix Diagram Data

    Attributes:
            @TwoDGibbsEnergy(pH,U_she): determine the excess Gibbs free surface energy of different surfaces at different pH and U
            @TwoDPhaseDetermination(lowpH,highpH,lowU,highU,pHgrids,Ugrids): determine the most stable surface phase at certain pH and U
        
    """
    def __init__(self,SurfSet):
        self.surfaceset = SurfSet
        self.k = 8.617333262e-5
        self.T = 298.15
        self.const = self.k * self.T * np.log(10)

        return

    def TwoDGibbsEnergy(self,pH,U_she):
        G = []
        surface = locals()

        for i in range(len(self.surfaceset)):

            name = str(self.surfaceset[i]["Name"])
            surface[name] = self.surfaceset[i]["G_constant"] - self.surfaceset[i]["k_index"] * U_she - self.surfaceset[i]["k_index"] * pH * self.const
            G.append([i+1,surface[name]])
        
        Gmin = 0
        index = 0 

        for k in range(len(G)):
            if G[k][1] <= Gmin:
                Gmin = G[k][1]
                index = G[k][0]

        return(index)

    def TwoDPhaseDetermination(self,lowpH,highpH,lowU,highU,pHgrids,Ugrids):
        self.pH = np.linspace(lowpH,highpH,pHgrids)
        self.U_she = np.linspace(lowU,highU,Ugrids)

        self.Phase = []

        for i in range(len(self.pH)):
            for j in range(len(self.U_she)):
                self.Phase.append([self.pH[i],self.U_she[i],self.TwoDGibbsEnergy(self.pH[i],self.U_she[j])])

        return




