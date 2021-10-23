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
- `atom_name`, required: the central atom e.g. `Co1` calculates τ for Co1
- `-e` `atom(s)`, optional: exclude atoms, e.g. `-e N1 N3` excludes bonds and angles to N1 and N3 from calculation 
- `-d` `N`, optional: excludes atoms outside `d = N Å` from calculation, e.g. `-d 2.1` excludes atoms with bond lengths larger than 2.1 Å from the central atom from calculation

## Remarks
- If the predicted coordination number is larger than 2, τ will be calculated independently of the real coordination geometry. 
- The suggestion τ<sub>4</sub> or τ<sub>5</sub> (<--) is based on the number of angles (6 for τ<sub>4</sub>, 10 for τ<sub>5</sub>).

## Know Issues
- The script is not very well tested with symmetry generated atom positions. However, this is rarely the case with small molecule structures.
- All flavors of τ are calculated as soon as two angles are present. So you have to check if it makes sense.
