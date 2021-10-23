# tau-calc
A Python3 script that calculates the geometry indices τ<sub>4</sub> and τ<sub>5</sub> of selected atoms from crystallographic information file (CIF). The script saves you the tedious checking of the two largest angles and the calculation of the τ-values.
If you use the τ<sub>4</sub> or τ<sub>5</sub> index to describe the coordination geometry of your compounds, please cite one or more of the following articles:

τ<sub>4</sub>:
> "Structural variation in copper(i) complexes with pyridylmethylamide ligands: 
>  structural analysis with a new four-coordinate geometry index, τ<sub>4</sub>"
>  
> Lei Yang, Douglas R. Powell, Robert P. Houser,
> *Dalton Trans.* **2007**, 955-964
> 
> DOI: https://doi.org/10.1039/B617136B

τ<sub>4</sub>' (τ<sub>4</sub> improved):
> "Coordination polymers and molecular structures among complexes of 
>  mercury(II) halides with selected 1-benzoylthioureas"
> 
> Andrzej Okuniewski, Damian Rosiak, Jarosław Chojnacki, Barbara Becker,
> *Polyhedron* **2015**, *90*, 47–57
> 
> DOI: https://doi.org/10.1016/j.poly.2015.01.035

τ<sub>5</sub>:
> "Synthesis, structure, and spectroscopic properties of copper(II) compounds containing 
>  nitrogen–sulphur donor ligands; the crystal and molecular structure of 
>  aqua[1,7-bis(N-methylbenzimidazol-2′-yl)-2,6-dithiaheptane]copper(II) perchlorate"
>  
> Anthony W. Addison, T. Nageswara Rao, Jan Reedijk, Jacobus van Rijn, Gerrit C. Verschoor, 
> *J. Chem. Soc., Dalton Trans.* **1984**, 1349-1356.
> 
> DOI: https://doi.org/10.1039/DT9840001349

The script uses the Gemmi library for CIF processing:

https://gemmi.readthedocs.io/en/latest/

https://github.com/project-gemmi/gemmi

## Introduction:
The geometry indices τ<sub>4</sub> and τ<sub>5</sub> help you in assigning a coordination geometry for four-coordinated (tetrahedral, trigonal pyramidal, square planar or seesaw geometry) or five-coordinates compounds (square pyramidal or trigonal bipyramidal geometry). Only the two largest angles enclosing the central atom are needed for the assignment. Please check the concise [Wikipedia article](https://en.wikipedia.org/wiki/Geometry_index) or the papers above for more information.

## External modules
 `re`,  `gemmi`
 
## Quick start
 Start the script with:
```console
python3 tau-calc.py filename.cif atom
```
The following output will be printed:
	Cu1 binds to:
	-----------------------------------------------------------------
	Cu1-O1 1.907(2) Å .
	Cu1-N1 1.911(3) Å .
	Cu1-N3 2.154(3) Å .
	Cu1-N5 2.170(3) Å .
	Cu1-N4 2.187(3) Å .
	
	The predicted coordination number for Cu1 is 5.
	
	Cu1 angles are:
	-----------------------------------------------------------------
	O1-Cu1-N1 98.68(12)° . .
	O1-Cu1-N3 171.94(12)° . .
	N1-Cu1-N3 80.09(13)° . .
	O1-Cu1-N5 103.92(12)° . .
	N1-Cu1-N5 137.78(13)° . .
	N3-Cu1-N5 81.87(12)° . .
	O1-Cu1-N4 91.77(11)° . .
	N1-Cu1-N4 130.58(13)° . .
	N3-Cu1-N4 83.17(12)° . .
	N5-Cu1-N4 84.14(13)° . .
	
	The two largest angles are beta = 171.94° and alpha = 137.78°.
	
	Cu1 geometry indices:
	-----------------------------------------------------------------
	tau_4  = 0.36 
	tau_4' = 0.25 
	tau_5  = 0.57 <--
	
	Table of typical geometries and their corresponding tau_x values
	-----------------------------------------------------------------
	Coordination number 4:
	Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00
	Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85
	Seesaw               : tau_4 = 0.42       tau_4' = 0.24
	Square planar        : tau_4 = 0.00       tau_4' = 0.00
	
	Coordination number 5:
	Trigonal bipyramidal : tau_5 = 1.00                     
	Square pyramidal     : tau_5 = 0.00 


