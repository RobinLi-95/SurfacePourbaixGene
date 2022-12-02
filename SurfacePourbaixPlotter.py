import matplotlib.pyplot as plt
import numpy as np

class SurfacePourbaixDiagramPlotter():
    """
    Class for generating 1-D (fixed pH or potential) and 2-D Surface Pourbaix Diagram Figures

    Attributes:
            @TwoDPlotter(pHRange,URange): Plot 2D Surface Pourbaix Diagram
        
    """

    def __init__(self,SurfPhase,pHgrids,Ugrids):
        SurfacePhase_tmp = np.array(SurfPhase)
        self.phase = SurfacePhase_tmp.reshape([pHgrids,Ugrids,3])

        return

    def TwoDPlotter(self,pHRange,URange):

        pHlow = pHRange[0]
        pHhigh = pHRange[1]
        Ulow = URange[0]
        Uhigh = URange[1]

        self.ratio = (pHhigh - pHlow)/(Uhigh - Ulow)

        fig, ax = plt.subplots(figsize=(6,6))
        ax.set_title('Pourbaix Diagram of Mg (0001) Surface')
        surf = ax.imshow(self.phase[:,:,2].T, extent=[pHlow,pHhigh,Ulow,Uhigh], origin='lower', cmap='magma')
        ax.set_aspect(self.ratio)
        ax.set_xlabel('pH')
        ax.set_ylabel('$U_{SHE} (V)$')
        # ax.hlines(-2.36, xmin=0, xmax=14,color='white',linestyles='dashed')
        
        cbar = fig.colorbar(surf, ax=ax, shrink=0.8,label='Surface Type')
        cbar.solids.set_edgecolor("face")
        fig.savefig('full_figure.png')

        return