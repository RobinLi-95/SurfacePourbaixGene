# SurfacePourbaixGene
### This is a set of scripts for generating Surface Pourbaix Diagram of Mg surfaces

1. <span style="color:red">*SlabBuilder*</span> is a tool for constructing different Mg surfaces based on relaxed bulk Mg structure.
2. <span style="color:red">*Mg10m10SiteFinder*</span> is a tool for determining the typical adsorption sites (onsurface:OT/HW/LB/SB and subsurface:Octa/Tetra-I/Tetra-II) in Mg (10m10) surface.
3. <span style="color:red">*Mg10m10AdsorptionHelper*</span> helps adsorbing -OH/-O/-H in the Mg (10m10) surfaces symmetrically.
4. <span style="color:red">*StructureOutput*</span> outputs the constructed surfaces in the VASP POSCAR format.
5. <span style="color:red">*VaspOutPutReader*</span> reads the data from VASP OUTCAR/POSCAR and determines important parameters for constructing surface Pourbaix diagram.
6. <span style="color:red">*SurfacePourbaixGenerator*</span> finds the most stable phase among a set of surfaces at certain pH and potential(U).
7. <span style="color:red">*SurfacePourbaixPlotter*</span> plots the 2D/1D(to be implemented) surface Pourbaix Diagram.
8. The <span style="color:blue">*Jupyter Notebooks*</span> (end with *ipynb) illustrate how to use these scripts. 
