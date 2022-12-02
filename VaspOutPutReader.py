from ase.io import read
import subprocess

class VaspOutputReader():
    """
    Class reading the VASP output file/data and initially generating essential parameter for the Pourbaix Diagram constructing

    Attributes:
            @SurfaceDetermine(): determine the total energy and chemical symbols of the surface
    """

    def __init__(self,filename):
        self.filename = filename
        self.structure = read(self.filename,format='vasp')
        command = "grep -A 4 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' OUTCAR | grep 'energy(sigma->0)' | tail -1 | awk '{print $7}'"
        self.toten = float(subprocess.getoutput(command))
        self.chemicalsymbol = self.structure.get_chemical_formula()

        return
  
    def SurfaceDetermine(self):
        self.symbols = self.structure.get_chemical_symbols()
        self.N_Mg = 0
        self.N_O = 0
        self.N_H = 0 
        self.E_H_ref = -9.87332962 / 2 # reference H total energy
        self.E_O_ref = -6.77273526 / 2
        self.Mg_10m10_clean = -91.27406327
        self.G_H2O = -2.46

        for i in range(len(self.symbols)):
            if self.symbols[i] == 'Mg':
                self.N_Mg += 1
            elif self.symbols[i] == 'O':
                self.N_O += 1
            elif self.symbols[i] == 'H':
                self.N_H += 1

        self.N_O  = int(self.N_O/2)
        self.N_H  = int(self.N_H/2)

        self.k_factor = 2 * self.N_O - self.N_H
        self.G0 = (self.toten - self.Mg_10m10_clean) / 2  - self.N_O * self.E_O_ref - self.N_H * self.E_H_ref - self.N_O * self.G_H2O

        return
