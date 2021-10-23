# tau-calc
A Python3 script that calculates the geometry indices τ<sub>4</sub> and τ<sub>5</sub> of selected atoms from crystallographic information file (CIF). The script saves you the tedious checking of the two largest angles and the calculation of the τ-values.
If you use the τ<sub>4</sub> or τ<sub>5</sub> index to describe the coordination geometry of your compounds, please cite one or more of the following articles:

**τ<sub>4</sub>**:
> "Structural variation in copper(i) complexes with pyridylmethylamide ligands: 
>  structural analysis with a new four-coordinate geometry index, τ<sub>4</sub>"
>  
> Lei Yang, Douglas R. Powell, Robert P. Houser,
> *Dalton Trans.* **2007**, 955-964.
> 
> DOI: https://doi.org/10.1039/B617136B

**τ<sub>4</sub>' (τ<sub>4</sub> improved)**:
> "Coordination polymers and molecular structures among complexes of 
>  mercury(II) halides with selected 1-benzoylthioureas"
> 
> Andrzej Okuniewski, Damian Rosiak, Jarosław Chojnacki, Barbara Becker,
> *Polyhedron* **2015**, *90*, 47–57.
> 
> DOI: https://doi.org/10.1016/j.poly.2015.01.035

**τ<sub>5</sub>**:
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
The input is case sensitive.

The following output will be printed:

	atom binds to:
	-----------------------------------------------------------------
	atom1-atom2 distance Å site symmetry
	...
	
	The predicted coordination number for atom is x.
	
	atom angles are:
	-----------------------------------------------------------------
	atom2-atom1-atom3 angle°  site symmetry  site symmetry
	...
	
	The two largest angles are beta = angle° and alpha = angle°.
	
	atom geometry indices:
	-----------------------------------------------------------------
	tau_4  = value 
	tau_4' = value 
	tau_5  = value <--
	
	Table of typical geometries and their corresponding tau_x values
	-----------------------------------------------------------------
	...
	
The suggested τ parameter is marked with an arrow (<--).

## Command-line options
- `filename` , required: the CIF (crystallographic information file)
- `atom_name`, required: the central atom, input is case sensitive, e.g. `Co1` calculates τ for Co1 
- `-e` `atom(s)`, optional: exclude atoms, e.g. `-e N1 N3` excludes bonds and angles to N1 and N3 from calculation 
- `-d` `N`, optional: excludes atoms outside `d = N Å` from calculation, e.g. `-d 2.1` excludes atoms with bond lengths larger than 2.1 Å from the central atom from calculation

## Remarks
- If the predicted coordination number is larger than 2, τ will be calculated independently of the real coordination geometry. 
- The suggestion τ<sub>4</sub> or τ<sub>5</sub> (<--) is based on the number of angles (6 for τ<sub>4</sub>, 10 for τ<sub>5</sub>).

## Know Issues
- The script is not very well tested with symmetry generated atom positions. However, this is rarely the case with small molecule structures.
- All flavors of τ are calculated as soon as two angles are present. So you have to check if it makes sense.

## Examples

### Example 1:
```console
python3 tau-calc.py test.cif Cu1
```
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

### Example 2:
```console
python3 tau-calc.py test2.cif Co1 -e N12 -d 2.0
```
	Excluded atoms: ['N12']
	
	Excluded atoms (distance larger than 2.0 Å): ['N11']
	
	Co1 binds to:
	-----------------------------------------------------------------
	Co1-N7 1.8860(15) Å .
	Co1-N8 1.8900(15) Å 3_666
	Co1-N14 1.9404(15) Å .
	Co1-N13 1.9502(15) Å .
	
	The predicted coordination number for Co1 is 4.
	
	Co1 angles are:
	-----------------------------------------------------------------
	N7-Co1-N8 92.29(6)° . 3_666
	N7-Co1-N14 91.55(6)° . .
	N8-Co1-N14 87.16(6)° 3_666 .
	N7-Co1-N13 85.35(6)° . .
	N8-Co1-N13 92.05(6)° 3_666 .
	N14-Co1-N13 176.77(6)° . .
	
	The two largest angles are beta = 176.77° and alpha = 92.29°.
	
	Co1 geometry indices:
	-----------------------------------------------------------------
	tau_4  = 0.64 <--
	tau_4' = 0.38 <--
	tau_5  = 1.41 
	
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
