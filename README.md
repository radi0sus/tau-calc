# tau-calc
A Python 3 script that calculates geometry indices, including τ<sub>4</sub>, τ<sub>5</sub>, O (octahedricity), and CShM (Continuous Shape Measures) values for various shapes of selected atoms from a crystallographic information file (CIF). The script saves you the tedious task of checking the two largest angles and calculating the τ-values. It also computes octahedricity and polyhedral volume and can print or save XYZ coordinates of the central atom and its neighboring atoms.

## Introduction:
The geometry indices τ<sub>4</sub> and τ<sub>5</sub> helps in assigning a coordination geometry for four-coordinated (tetrahedral, trigonal pyramidal, square planar or seesaw geometry) or five-coordinates compounds (square pyramidal or trigonal bipyramidal geometry). Only the two largest angles enclosing the central atom are needed for the assignment. Please check the concise [Wikipedia article](https://en.wikipedia.org/wiki/Geometry_index) or the papers above for more information. 
The octahedricity *O* is calculated with the following equation (see paper for more information):

$$ O = \sqrt{\frac{1}{15}\sum_{i=1}^{15}(\hat{\theta_i} - \theta)^2} $$ 

$\hat{\theta_i}$ = 180° for *trans* X-M-X angles and 90° for *cis* X-M-X angles\
$\theta$ = experimental X-M-X angles 

And is close to zero for an almost ideal octahedron

The CShM (Continuous Shape Measures) value approaches zero when a shape closely resembles the chosen ideal shape. For further details, refer to the paper and related articles on CShM.

The polyhedral volume is calculated using the Convex Hull Algorithm (implemented via SciPy and the [Qhull Library](http://www.qhull.org).

## External modules
 `gemmi`, `numpy`, `scipy`
 
## Quick start
 Start the script with:
```console
python3 tau-calc.py filename.cif atom
```
The input is case sensitive.

The following output will be printed:

    The predicted coordination number for atom is value.
    
    The two largest angles are beta = value and alpha = value.
    Note: Angles for the calculation of tau_4, tau_4' and tau_5.
    
    Number of cis angles ~ 90° (< 135°)    = value
    Number of trans angles ~ 180° (> 135°) = value
    Note: First value should be 12, second 3 for an octahedron.
    
    --------------------------------------------------------------------------------
    atom geometry indices  ("<--" indicates the likely structural parameter):
    --------------------------------------------------------------------------------
    tau_4  =  value 
    tau_4' =  value 
    tau_5  =  value <--
    O      =  value
    
    --------------------------------------------------------------------------------
    Continuous shape measure (CShM):
    --------------------------------------------------------------------------------
    CShM S(Shape 1) = value ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    CShM S(Shape 2) = value ░░░░
    CShM S(Shape 3) = value ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    --------------------------------------------------------------------------------
    Polyhedral volume = value A³
    --------------------------------------------------------------------------------
    
    --------------------------------------------------------------------------------
    Table of typical geometries and their corresponding tau_x and O values: 
    --------------------------------------------------------------------------------
    Coordination number 4:
    Tetrahedral          :  tau_4 = 1.00       tau_4' = 1.00
    Trigonal pyramidal   :  tau_4 = 0.85       tau_4' = 0.85
    Seesaw               :  tau_4 = 0.43       tau_4' = 0.24
    Square planar        :  tau_4 = 0.00       tau_4' = 0.00
    
    Coordination number 5:
    Trigonal bipyramidal :  tau_5 = 1.00                     
    Square pyramidal     :  tau_5 = 0.00                    
    
    Coordination number 6:
    Ideal octahedron     :      O = 0.00 
	
The likely structural parameter is marked with an arrow (<--).

## Command-line options
- `filename` , required: the CIF (crystallographic information file)
- `atom_name`, required: the central atom, input is case sensitive, e.g. `Co1` calculates τ for Co1 
- `-e` `atom(s)`, optional: exclude atoms, e.g. `-e N1 N3` excludes bonds and angles to N1 and N3 from calculation 
- `-d` `N`, optional: excludes atoms outside `d = N Å` from calculation, e.g. `-d 2.1` excludes atoms with bond lengths larger than 2.1 Å from the central atom from calculation
- `-v` optional: verbose output, prints all bond lengths and angles of the central atom to the neighboring atoms and the XYZ coordinates
- `-sxyz` optional: save the XYZ coordinates of the central atom and its neighboring atoms (filename: `cif_name-atom_name.xyz`)
  
## Remarks
- If the predicted coordination number is larger than 2, τ will be calculated independently of the real coordination geometry. 
- The suggestion τ<sub>4</sub>, τ<sub>5</sub> or *O* (<--) is based on the number of angles (6 for τ<sub>4</sub>, 10 for τ<sub>5</sub>, 15 for *O*).
- In octahedral coordination, all angles less than 135° are considered "90°" *cis* angles, while angles greater than 135° are considered "180°" *trans* angles.
- The XYZ coordinates of the neighboring atoms are given relative to the central atom, which is positioned at [0, 0, 0].
- The XYZ file (option: `-sxyz`) is useful for further studies of coordination geometry, such as with the Continuous Shape Measures (CShM) method.
- The polyhedral volume should match the value calculated by [Olex2](https://www.olexsys.org/olex2/).
- More reference shapes for CShM can be defined in the code. You can find definitions of several ideal shapes  
[here (cosymlib)](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml) or
[here (shape)](https://github.com/continuous-symmetry-measure/shape). Ideal shapes that include `'-'` are from the first source, while other shapes come from the second.
- The script's results should match those of the [online calculator](https://csm.ouproj.org.il/molecule) or the [Shape program](https://www.ee.ub.edu/downloads/) when using default options.
- Although a strictly defined τ<sub>3</sub> value does not exist, shape measures for 3-coordinated atoms are still calculated.
- If the calculation is not fast enough for your needs, you can uncomment the optimization process using the Hungarian algorithm. 
  Ensure that the number of trials (`num_trials`) is set high enough to avoid missing the global minimum.
  
## Known Issues
- The script is not very well tested with symmetry generated atom positions. However, this is rarely the case with small molecule structures.
- All flavors of τ and *O* are calculated as soon as two angles are present. So you have to check if it makes sense.
- The script can only remove atoms from the coordination sphere, not add atoms. Therefore, make sure that the connectivity list is appropriate.
- For the generation of the XYZ file, hydrogen atoms are generally ignored. As a result, in metal hydrides, hydrogen atoms will not be included in the XYZ file.
- The CShM method is rewritten from the [C++ code](https://github.com/continuous-symmetry-measure/shape) and may still contain errors.
- As the number of vertices increases, both computation time and memory usage grow rapidly ($n!$).
  To address this issue, enable the optimization process using the Hungarian algorithm.
  Ensure that the number of trials (`num_trials`) is set sufficiently high to avoid missing the global minimum.
  
## Examples

### Example 1:
```console
python3 tau-calc.py test.cif Hg1
```
    The predicted coordination number for Hg1 is 4.
    
    The two largest angles are beta = 135.94° and alpha = 124.44°.
    Note: Angles for the calculation of tau_4, tau_4' and tau_5.
     
    Number of cis angles ~ 90° (< 135°)    = 5
    Number of trans angles ~ 180° (> 135°) = 1
    Note: First value should be 12, second 3 for an octahedron.
     
    --------------------------------------------------------------------------------
    Hg1 geometry indices ("<--" indicates the likely structural parameter):
    --------------------------------------------------------------------------------
    tau_4  =   0.71 <--
    tau_4' =   0.67 <--
    tau_5  =   0.19 
    O      =  23.49 
     
    --------------------------------------------------------------------------------
    Continuous shape measure (CShM):
    --------------------------------------------------------------------------------
    S(SP-4, Square)                               =  32.0294 ░░░░░░░░░░░░░░░░░░░░░░░
    S(T-4, Tetrahedron)                           =   2.8099 ░
    S(SS-4, Seesaw or sawhorse)                   =   5.8584 ░░░
    S(vTBPY-4, Axially vacant trigonal bipyramid) =   2.8078 ░
     
    --------------------------------------------------------------------------------
    Polyhedral volume = 9.7478 A³
    --------------------------------------------------------------------------------
     
    --------------------------------------------------------------------------------
    Table of typical geometries and their corresponding tau_x and O values: 
    --------------------------------------------------------------------------------
    Coordination number 4:
    Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00
    Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85
    Seesaw               : tau_4 = 0.43       tau_4' = 0.24
    Square planar        : tau_4 = 0.00       tau_4' = 0.00
    
    Coordination number 5:
    Trigonal bipyramidal : tau_5 = 1.00                     
    Square pyramidal     : tau_5 = 0.00                    
    
    Coordination number 6:
    Ideal octahedron     :     O = 0.00  
	
### Example 2:
```console
python3 tau-calc.py test2.cif Ru1 -v
```
    --------------------------------------------------------------------------------
    Ru1 binds to:
    --------------------------------------------------------------------------------
    N1-Ru1 2.065(3) Å .
    Ru1-N1 2.065(3) Å 3_665
    Ru1-N1 2.065(3) Å 9_655
    Ru1-N1 2.065(3) Å 6_655
    Ru1-N1 2.065(3) Å 12_665
    Ru1-N1 2.065(3) Å 7
     
    The predicted coordination number for Ru1 is 6.
    
    --------------------------------------------------------------------------------
    Ru1 angles are:
    --------------------------------------------------------------------------------
    N1-Ru1-N1 78.98(19)° . 3_665
    N1-Ru1-N1 94.66(13)° . 9_655
    N1-Ru1-N1 92.37(18)° 3_665 9_655
    N1-Ru1-N1 92.37(18)° . 6_655
    N1-Ru1-N1 94.66(13)° 3_665 6_655
    N1-Ru1-N1 170.89(18)° 9_655 6_655
    N1-Ru1-N1 94.66(13)° . 12_665
    N1-Ru1-N1 170.89(18)° 3_665 12_665
    N1-Ru1-N1 94.66(13)° 9_655 12_665
    N1-Ru1-N1 78.98(19)° 6_655 12_665
    N1-Ru1-N1 170.89(18)° . 7
    N1-Ru1-N1 94.66(13)° 3_665 7
    N1-Ru1-N1 78.98(19)° 9_655 7
    N1-Ru1-N1 94.66(13)° 6_655 7
    N1-Ru1-N1 92.37(18)° 12_665 7
     
    The two largest angles are beta = 170.89° and alpha = 170.89°.
    Note: Angles for the calculation of tau_4, tau_4' and tau_5.
     
    Number of cis angles ~ 90° (< 135°)    = 12
    Number of trans angles ~ 180° (> 135°) = 3
    Note: First value should be 12, second 3 for an octahedron.
     
    --------------------------------------------------------------------------------
    Ru1 geometry indices ("<--" indicates the likely structural parameter):
    --------------------------------------------------------------------------------
    tau_4  =   0.13 
    tau_4' =   0.13 
    tau_5  =   0.00 
    O      =   7.12 <--
     
    --------------------------------------------------------------------------------
    Continuous shape measure (CShM):
    --------------------------------------------------------------------------------
    S(HP-6, Hexagon)                           =  28.5605 ░░░░░░░░░░░░░░░░░░░░░░░░
    S(PPY-6, Pentagonal pyramid)               =  27.0546 ░░░░░░░░░░░░░░░░░░░░░░
    S(OC-6, Octahedron)                        =   0.9516 ░
    S(TPR-6, Trigonal prism)                   =  13.6312 ░░░░░░░░░░░
    S(JPPY-6, Johnson pentagonal pyramid (J2)) =  30.8553 ░░░░░░░░░░░░░░░░░░░░░░░░░░
     
    --------------------------------------------------------------------------------
    Polyhedral volume = 11.5171 A³
    --------------------------------------------------------------------------------
     
    --------------------------------------------------------------------------------
    Table of typical geometries and their corresponding tau_x and O values: 
    --------------------------------------------------------------------------------
    Coordination number 4:
    Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00
    Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85
    Seesaw               : tau_4 = 0.43       tau_4' = 0.24
    Square planar        : tau_4 = 0.00       tau_4' = 0.00
    
    Coordination number 5:
    Trigonal bipyramidal : tau_5 = 1.00                     
    Square pyramidal     : tau_5 = 0.00                    
    
    Coordination number 6:
    Ideal octahedron     :     O = 0.00                      
     
    --------------------------------------------------------------------------------
    XYZ coordinates of the central atom and its neighbors: 
    --------------------------------------------------------------------------------
    7
    879418.cif Ru1
    Ru  0.00000000  0.00000000  0.00000000
    N   1.01427676 -1.42908956 -1.09132010
    N   0.72983256  1.59407203 -1.09132010
    N  -1.74608007 -0.16384466 -1.09132010
    N   1.74476624 -0.16384466  1.09132010
    N  -1.01559059 -1.42908956  1.09132010
    N  -0.73114639  1.59407203  1.09132010

### Example 3:
```console
python3 tau-calc.py test3.cif Co1 -e N12
```
    Excluded atoms: ['N12']
     
    The predicted coordination number for Co1 is 5.
    
    The two largest angles are beta = 176.77° and alpha = 173.52°.
    Note: Angles for the calculation of tau_4, tau_4' and tau_5.
     
    Number of cis angles ~ 90° (< 135°)    = 8
    Number of trans angles ~ 180° (> 135°) = 2
    Note: First value should be 12, second 3 for an octahedron.
     
    --------------------------------------------------------------------------------
    Co1 geometry indices ("<--" indicates the likely structural parameter):
    --------------------------------------------------------------------------------
    tau_4  =   0.07 
    tau_4' =   0.06 
    tau_5  =   0.05 <--
    O      =   4.65 
     
    --------------------------------------------------------------------------------
    Continuous shape measure (CShM):
    --------------------------------------------------------------------------------
    S(PP-5, Pentagon)                            =  30.9225 ░░░░░░░░░░░░░░░░░░░░░░░░
    S(vOC-5, Vacant octahedron (J1))             =   0.3515 ░
    S(TBPY-5, Trigonal bipyramid)                =   6.6922 ░░░░░
    S(SPY-5, Square pyramid)                     =   2.2513 ░░
    S(JTBPY-5, Johnson trigonal bipyramid (J12)) =   7.9883 ░░░░░░
     
    --------------------------------------------------------------------------------
    Polyhedral volume = 4.8918 A³
    --------------------------------------------------------------------------------
     
    --------------------------------------------------------------------------------
    Table of typical geometries and their corresponding tau_x and O values: 
    --------------------------------------------------------------------------------
    Coordination number 4:
    Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00
    Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85
    Seesaw               : tau_4 = 0.43       tau_4' = 0.24
    Square planar        : tau_4 = 0.00       tau_4' = 0.00
    
    Coordination number 5:
    Trigonal bipyramidal : tau_5 = 1.00                     
    Square pyramidal     : tau_5 = 0.00                    
    
    Coordination number 6:
    Ideal octahedron     :     O = 0.00

## References
If you use τ<sub>4</sub>, τ<sub>5</sub>, the *O* index or CShM to describe the coordination geometry of your compounds, please cite one or more of the following articles:

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

***O***:
> "Structural, electrochemical and photophysical behavior of Ru(II) complexes with 
>  large bite angle sulfur-bridged terpyridyl ligands"
>  
> Christopher M. Brown, Nicole E. Arsenault, Trevor N. K. Cross, Duane Hean, Zhen Xu, 
> Michael O. Wolf, 
> *Inorg. Chem. Front.* **2020**, *7*, 117-127.
> 
> DOI: https://doi.org/10.1039/C9QI01009B

**CShM**:
> "Continuous Symmetry Measures. 5. The Classical Polyhedra"
>  
> Mark Pinsky, David Avnir, 
> *Inorg. Chem.* **1998**, *37*, 5575–5582.
> 
> DOI: https://doi.org/10.1021/ic9804925
> 
> "Shape maps and polyhedral interconversion paths in transition metal chemistry"
>  
> Santiago Alvarez, Pere Alemany, David Casanova, Jordi Cirera, Miquel Llunell, David Avnir,
> *Coord. Chem. Rev.*, **2005**, *249*, 1693–1708.
> 
> DOI: https://doi.org/10.1016/j.ccr.2005.03.031

The script uses the **Gemmi** library for CIF processing:
> "GEMMI: A library for structural biology"
> 
> Marcin Wojdyr,
> *Journal of Open Source Software* **2022**, *7 (73)*, 4200.
>
> DOI: https://doi.org/10.1021/ic9804925

https://gemmi.readthedocs.io/en/latest/

https://github.com/project-gemmi/gemmi
