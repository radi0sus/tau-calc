#!/usr/bin/env python3
#
# Calculation of tau_4, tau_4', tau_5, or O, several CShM
# for 3-, 4-, 5- or 6-coordinated atoms and the polyhedral volume.
# For a deeper explanation of tau_4 and tau_5 have a look at Wikipedia:
#
# https://en.wikipedia.org/wiki/Geometry_index
# 
# and the following articles:
#
# For tau_4, please cite:
#
# "Structural variation in copper(i) complexes with pyridylmethylamide ligands: 
#  structural analysis with a new four-coordinate geometry index, τ4"
#
# Lei Yang, Douglas R. Powell, Robert P. Houser
# Dalton Trans. 2007, 955-964.
# DOI: https://doi.org/10.1039/B617136B
#
# For the improved tau_4', please cite:
#
# "Coordination polymers and molecular structures among complexes of 
#  mercury(II) halides with selected 1-benzoylthioureas"
#
# Andrzej Okuniewski, Damian Rosiak, Jarosław Chojnacki, Barbara Becker
# Polyhedron 2015, 90, 47–57.
# DOI: https://doi.org/10.1016/j.poly.2015.01.035
#
# For tau_5, please cite:
#
# "Synthesis, structure, and spectroscopic properties of copper(II) compounds containing 
#  nitrogen–sulphur donor ligands; the crystal and molecular structure of 
#  aqua[1,7-bis(N-methylbenzimidazol-2′-yl)-2,6-dithiaheptane]copper(II) perchlorate"
#
# Anthony W. Addison, T. Nageswara Rao, Jan Reedijk, Jacobus van Rijn, Gerrit C. Verschoor  
# J. Chem. Soc., Dalton Trans. 1984, 1349-1356.
# DOI: https://doi.org/10.1039/DT9840001349
#
# For octahedricity O, please cite:
#
# "Structural, electrochemical and photophysical behavior of Ru(ii) complexes with 
# large bite angle sulfur-bridged terpyridyl ligands"
# 
# Christopher M. Brown, Nicole E. Arsenault, Trevor N. K. Cross, Duane Hean, Zhen Xu, 
# Michael O. Wolf
# Inorg. Chem. Front. 2020, 7, 117-127. 
# DOI: https://doi.org/10.1039/C9QI01009B 
# 
# For CShM (Continuous Shape Measures) please cite:
#
# "Continuous Symmetry Measures. 5. The Classical Polyhedra"
# 
# Mark Pinsky, David Avnir
# Inorg. Chem. 1998, 37, 21, 5575–5582.
# DOI: https://doi.org/10.1021/ic9804925
#
# "Shape maps and polyhedral interconversion paths in transition metal chemistry"
#
# Santiago Alvarez, Pere Alemany, David Casanova, Jordi Cirera, Miquel Llunell, 
# David Avnir,
# Coord. Chem. Rev., 2005, 249, 1693–1708.
# DOI: https://doi.org/10.1016/j.ccr.2005.03.031
# 
# For Gemmi please cite: 
#
# "GEMMI: A library for structural biology"
# 
# Marcin Wojdyr
# Journal of Open Source Software 2022, 7 (73), 4200
# DOI: https://doi.org/10.21105/joss.04200
#
# https://gemmi.readthedocs.io/en/latest/
# https://github.com/project-gemmi/gemmi
#

import sys                                            # sys
import os                                             # for save xyz
import argparse                                       # argument parser
import re                                             # regular expressions
import gemmi                                          # CIF processing, coordinates
import numpy as np                                    # all sort of math
from scipy.linalg import svd                          # SVD for CShM
from itertools import permutations                    # permutations of vertices for CShM
from scipy.spatial import ConvexHull                  # polyhedral volume

#regex for bonds and angles --> string to float, no esd
ang_bond_val = re.compile(r'\d{1,}[\.]?\d{0,}')

#list for the exclusion of atoms in the angle table
list_of_atoms_with_large_bonds=[]

# Definitions for several Shapes START ##################
# from: https://github.com/continuous-symmetry-measure/shape/blob/master/src/shape.cpp
# for CShM (Continuous Shape Measures)
# AB6
# define the ideal octahedron with center
#ao = 1.0 / np.sqrt(2.0)
#IDEAL_AB6 = np.array([
#    [  0.0, 0.0,   1.0],
#    [  ao,   ao,   0.0],
#    [ -ao,   ao,   0.0],
#    [ -ao,  -ao,   0.0],
#    [  ao,  -ao,   0.0],
#    [ 0.0,  0.0,  -1.0],
#    [ 0.0,  0.0,   0.0]
#])

# for CShM (Continuous Shape Measures)
# APR
# define the ideal trigonal prism with center 
#ap = 1.0 / np.sqrt(2.0)
#bp = 1.0 / np.sqrt(3.0)
#cp = 1.0 / np.sqrt(6.0)
#dp = 2.0 * cp

#IDEAL_APR = np.array([
#    [ 0.0, -dp,  bp],
#    [-ap,   cp,  bp],
#    [ ap,   cp,  bp],
#    [ 0.0, -dp, -bp],
#    [-ap,   cp, -bp],
#    [ ap,   cp, -bp],
#    [ 0.0, 0.0, 0.0]  
#])

# APR_EQ
# define the ideal trigonal equilateral prism with center 
# SAME as TPR-6
#IDEAL_APR_EQ = np.array([
#    [ 0.0, -dp,  ap],
#    [-ap,   cp,  ap],
#    [ ap,   cp,  ap],
#    [ 0.0, -dp, -ap],
#    [-ap,   cp, -ap],
#    [ ap,   cp, -ap],
#    [ 0.0, 0.0, 0.0]  
#])

# for CShM (Continuous Shape Measures)
# OC-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Octahedron
# same as AB6
IDEAL_OC6 = np.array([
    [ 0.0           ,  0.0           , -1.080123449735],
    [ 1.080123449735,  0.0           ,  0.0           ],
    [ 0.0           ,  1.080123449735,  0.0           ],
    [-1.080123449735,  0.0,             0.0           ],
    [ 0.0           , -1.080123449735,  0.0           ],
    [ 0.0           ,  0.0           ,  1.080123449735],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# TPR-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal trigonal prism
# same as APR_EQ
IDEAL_TPR6 = np.array([
    [ 0.816496580928,  0.0           , -0.707106781187],
    [-0.408248290464,  0.707106781187, -0.707106781187],
    [-0.408248290464, -0.707106781187, -0.707106781187],
    [ 0.816496580928,  0.0           ,  0.707106781187],
    [-0.408248290464,  0.707106781187,  0.707106781187],
    [-0.408248290464, -0.707106781187,  0.707106781187],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# JJPY-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Johnson pentagonal pyramid (J2)
IDEAL_JPPY6 = np.array([
    [ 1.146281780821,  0.0           ,  0.101205871605],
    [ 0.354220550616,  1.090178757161,  0.101205871605],
    [-0.927361441027,  0.673767525738,  0.101205871605],
    [-0.927361441027, -0.673767525738,  0.101205871605],
    [ 0.354220550616, -1.090178757161,  0.101205871605],
    [ 0.0           ,  0.0           , -0.607235229628],
    [ 0.0           ,  0.0           ,  0.101205871605]
])

# for CShM (Continuous Shape Measures)
# HP-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Hexagon
IDEAL_HP6 = np.array([
    [ 1.080123449735,  0.0           , 0.0],
    [ 0.540061724867,  0.935414346693, 0.0],
    [-0.540061724867,  0.935414346693, 0.0],
    [-1.080123449735,  0.0           , 0.0],
    [-0.540061724867, -0.935414346693, 0.0],
    [ 0.540061724867, -0.935414346693, 0.0],
    [ 0.0           ,  0.0           , 0.0],
])

# for CShM (Continuous Shape Measures)
# PPY-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Pentagonal pyramid
IDEAL_PPY6 =  np.array([
    [ 0.0           ,  0.0           , -0.937042571332],
    [ 1.093216333220,  0.0           ,  0.156173761889],
    [ 0.337822425493,  1.039710517429,  0.156173761889],
    [-0.884430592103,  0.642576438232,  0.156173761889],
    [-0.884430592103, -0.642576438232,  0.156173761889],
    [ 0.337822425493, -1.039710517429,  0.156173761889],
    [ 0.0           ,  0.0           ,  0.156173761889]
])
# for CShM (Continuous Shape Measures)
# AB5
# define the ideal bipyramid with center 
#ab = np.sqrt(3.0/8.0)
#bb = 1 / np.sqrt(8.0)
#cb = 1 / np.sqrt(2.0)
#
#IDEAL_AB5 = np.array([
#    [ 0.0, 0.0,  1.0],
#    [ -ab, -bb,  0.0],
#    [  ab, -bb,  0.0],
#    [ 0.0,  cb,  0.0],
#    [ 0.0, 0.0, -1.0],
#    [ 0.0, 0.0,  0.0]  
#])

# for CShM (Continuous Shape Measures)
# AB5_
# define the ideal bipyramid with center 
# same as TBPY-5
#ab_ = np.sqrt(3.0) / 2.0
#IDEAL_AB5_ = np.array([
#    [ 0.0,  0.0,  1.0],
#    [-ab_, -0.5,  0.0],
#    [ ab_, -0.5,  0.0],
#    [ 0.0,  1.0,  0.0],
#    [ 0.0,  0.0, -1.0],
#    [ 0.0,  0.0,  0.0]  
#])

# for CShM (Continuous Shape Measures)
# SPY-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal square pyramid
IDEAL_SPY5 = np.array([
    [ 0.0           ,  0.0           ,      1.095445115010],
    [ 1.060660171780,  0.0           ,     -0.273861278753],
    [ 0.0           ,  1.060660171780,     -0.273861278753],
    [-1.060660171780,  0.0           ,     -0.273861278753],
    [ 0.0           , -1.060660171780,     -0.273861278753],
    [ 0.0           ,  0.0           ,      0.0           ]
])

# for CShM (Continuous Shape Measures)
# TBPY-5 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal trigonal bipyramid
# same as AB5_
IDEAL_TBPY5 = np.array([
    [ 0.0,             0.0,            -1.095445115010],
    [ 1.095445115010,  0.0,             0.0           ],
    [-0.547722557505,  0.948683298051,  0.0           ],
    [-0.547722557505, -0.948683298051,  0.0           ],
    [ 0.0,             0.0,             1.095445115010],
    [ 0.0,             0.0,             0.0           ],
])

# for CShM (Continuous Shape Measures)
# vOC-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal vacant octahedron (Johnson square pyramid, J1)
IDEAL_vOC5 = np.array([
    [ 0.0           ,  0.0           , -0.928476690885],
    [ 1.114172029062,  0.0           ,  0.185695338177],
    [ 0.0           ,  1.114172029062,  0.185695338177],
    [-1.114172029062,  0.0           ,  0.185695338177],
    [ 0.0           , -1.114172029062,  0.185695338177],
    [ 0.0           ,  0.0           ,  0.185695338177],
])

# for CShM (Continuous Shape Measures)
# TP-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Pentagon
IDEAL_PP5 = np.array([
    [ 1.095445115010,  0.0           , 0.0],
    [ 0.338511156943,  1.041830214874, 0.0],
    [-0.886233714448,  0.643886483299, 0.0],
    [-0.886233714448, -0.643886483299, 0.0],
    [ 0.338511156943, -1.041830214874, 0.0],
    [ 0.0           ,  0.0           , 0.0]
])

# for CShM (Continuous Shape Measures)
# JTBPY-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Johnson trigonal bipyramid (J12)
IDEAL_JTBPY5 = np.array([
    [ 0.925820099773,  0.0           ,  0.0           ],
    [-0.462910049886,  0.801783725737,  0.0           ],
    [-0.462910049886, -0.801783725737,  0.0           ],
    [ 0.0           ,  0.0           ,  1.309307341416],
    [ 0.0           ,  0.0           , -1.309307341416],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# AB4
# define the ideal tetrahedron with center 
#at = np.sqrt(8.0) / 3.0
#bt = 1.0 / 3.0
#ct = np.sqrt(2.0 / 3.0)
#dt = np.sqrt(2.0) / 3.0
#IDEAL_AB4 = np.array([
#    [ 0.0,  0.0,  1.0],
#    [ 0.0,   at,  -bt],
#    [ ct,   -dt,  -bt],
#    [-ct,   -dt,  -bt],
#    [ 0.0,  0.0,  0.0]  
#])

# for CShM (Continuous Shape Measures)
# SQ5
# define the ideal square with center 
#asq = 1 / np.sqrt(2.0)
#
#IDEAL_SQ5 = np.array([
#    [ asq,  asq,  0.0     ],
#    [ asq, -asq,  0.0     ],
#    [-asq, -asq,  0.0     ],
#    [-asq,  asq, -0.000001],
#    [ 0.0,  0.0,  0.000001]  
#])

# for CShM (Continuous Shape Measures)
# T-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal tetrahedron 
# same as AB4
IDEAL_T4 = np.array([
    [ 0.0           ,  0.912870929175, -0.645497224368],
    [ 0.0           , -0.912870929175, -0.645497224368],
    [ 0.912870929175,  0.0           ,  0.645497224368],
    [-0.912870929175,  0.0           ,  0.645497224368],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# SP-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal square 
# same as SQ5
IDEAL_SP4 = np.array([
    [ 1.118033988750,  0.0           , 0.0],
    [ 0.0           ,  1.118033988750, 0.0],
    [-1.118033988750,  0.0           , 0.0],
    [ 0.0           , -1.118033988750, 0.0],
    [ 0.0           ,  0.0           , 0.0],
])

# for CShM (Continuous Shape Measures)
# SS-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal seesaw with center 
IDEAL_SS4 = np.array([
    [-0.235702260396, -0.235702260396, -1.178511301978],
    [ 0.942809041582, -0.235702260396,  0.0           ],
    [-0.235702260396,  0.942809041582,  0.0           ],
    [-0.235702260396, -0.235702260396,  1.178511301978],
    [-0.235702260396, -0.235702260396,  0.0           ]#,
   #[ 0.0,             0.0,             0.0           ]
])

# for CShM (Continuous Shape Measures)
# vTBPY-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the axially vacant trigonal bipyramid with center (trigonal pyramidal)
IDEAL_vTBPY4 = np.array([
    [ 0.0,             0.0,            -0.917662935482],
    [ 1.147078669353,  0.0,             0.229415733871],
    [-0.573539334676,  0.993399267799,  0.229415733871],
    [-0.573539334676, -0.993399267799,  0.229415733871],
    [ 0.0,             0.0,             0.229415733871]#,
   #[ 0.0,             0.0,             0.0           ]
])

# for CShM (Continuous Shape Measures)
# TP-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for trigonal planar
IDEAL_TP3 = np.array([
    [ 1.154700538379,  0.0, 0.0],
    [-0.577350269190,  1.0, 0.0],
    [-0.577350269190, -1.0, 0.0],
    [ 0.0           ,  0.0, 0.0]
])

# for CShM (Continuous Shape Measures)
# vT-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for Pyramid (vacant tetrahedron)
IDEAL_vT3 = np.array([
    [ 1.137070487230,  0.0,             0.100503781526],
    [-0.568535243615,  0.984731927835,  0.100503781526],
    [-0.568535243615, -0.984731927835,  0.100503781526],
    [ 0.0,             0.0,            -0.301511344578]
])

# for CShM (Continuous Shape Measures)
# fac-vOC-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for fac-Trivacant octahedron
IDEAL_facvOC3 = np.array([
    [ 1.0,            -0.333333333333, -0.333333333333],
    [-0.333333333333,  1.0,            -0.333333333333],
    [-0.333333333333, -0.333333333333,  1.0           ],
    [-0.333333333333, -0.333333333333, -0.333333333333]
])

# for CShM (Continuous Shape Measures)
# mer-vOC-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for mer-Trivacant octahedron (T-shape)
IDEAL_mvOC3 = np.array([
    [ 1.206045378311, -0.301511344578, 0.0],
    [ 0.0,             0.904534033733, 0.0],
    [-1.206045378311, -0.301511344578, 0.0],
    [ 0.0,            -0.301511344578, 0.0]
])
# Definitions for several Shapes END ##################

#calculation of tau5
def calc_tau5(beta, alpha):
    tau5 = (beta - alpha) / 60
    return tau5

#calculation of tau4 classic
def calc_tau4(beta, alpha):
    tau4 = (360 - (alpha + beta )) / (360 - 2*109.5)
    return tau4

#calculation of tau4' improved
def calc_tau4impr(beta, alpha):
    tau4impr = (beta - alpha) / (360 - 109.5) + (180 - beta) / (180 - 109.5)
    return tau4impr

#calulation of the octahedricity O
def calc_octahedricity(measured_angles):
    # O (rms angle deviation) = sqrt ((1/15)*sum(omega_measured - omega_ideal)^2)
    # ideal angle is 90 deg if angle is < 135 deg and 180 deg if angle is > 135 deg
    ideal_angles = [90 if angle <= 135 else 180 for angle in measured_angles]
    # deviation of the measured angles from ideal angles: delta omega
    # deviations = [(measured - ideal) for measured, ideal in zip(measured_angles, ideal_angles)]
    deviations = [(ideal - measured) for ideal, measured in zip(ideal_angles, measured_angles)]
    # squared deviations: delta omega^2
    squared_deviations = [dev**2 for dev in deviations]
    # octahedricity
    octahedricity = np.sqrt(sum(squared_deviations) / len(squared_deviations))
    return octahedricity

# center and normalize the structure for CShM calculations
def normalize_structure(coordinates):
    centered_coords = coordinates - np.mean(coordinates, axis=0)
    norm = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))
    return centered_coords / norm

# calculation the continuous shape measures (CShM) parameter for a given structure
# from the c++ code with some help of AI
# https://github.com/continuous-symmetry-measure/shape
# all permutations considered
def calc_cshm(coordinates, ideal_shape):
    permut_list = list(permutations(range(len(coordinates))))
    ideal_sq_norms = np.sum(ideal_shape**2)
    input_structure = normalize_structure(coordinates)
    min_cshm = float('inf')
    for permuted_ideal in map(lambda p: ideal_shape[list(p)], permut_list):
        H = np.dot(input_structure.T, permuted_ideal)
        U, _, Vt = svd(H)
        R = np.dot(Vt.T, U.T)

        rotated_ideal = np.dot(permuted_ideal, R)
        scale = np.sum(input_structure * rotated_ideal) / ideal_sq_norms
        cshm = np.mean(np.sum((input_structure - scale * rotated_ideal)**2, axis=1))
        
        min_cshm = min(min_cshm, cshm)
        
    return min_cshm * 100

# faster Hungarian algorithm optimization
# check number of trials, if it is to low, it calculates the
# local and not the global minimum
#def calc_cshm(coordinates, ideal_shape, num_trials = 24):
#    from scipy.optimize import linear_sum_assignment
#    input_structure = normalize_structure(coordinates)
#    ideal_sq_norms = np.sum(ideal_shape**2)
#    
#    min_cshm = float('inf')
#    
#    # Generate some initial rotations to avoid local minima
#    # This simulates checking multiple permutations like in the exhaustive approach
#    for trial in range(num_trials):
#       if trial == 0:
#            # First trial with identity rotation
#            R_init = np.eye(3)
#        else:
#            # Random rotation matrix for subsequent trials
#            # Generate a random rotation matrix using QR decomposition
#            A = np.random.randn(3, 3)
#            Q, _ = np.linalg.qr(A)
#            R_init = Q
#            
#            # Ensure it's a proper rotation (det=1)
#            if np.linalg.det(R_init) < 0:
#                R_init[:, 0] *= -1
#                
#        # Apply initial rotation to ideal shape
#        rotated_ideal_init = np.dot(ideal_shape, R_init)
#        
#        # Compute cost matrix based on squared Euclidean distances
#        cost_matrix = np.linalg.norm(input_structure[:, None, :] - rotated_ideal_init[None, :, :], axis=2)
#        
#        # Solve assignment problem (Hungarian algorithm)
#        row_ind, col_ind = linear_sum_assignment(cost_matrix)
#        
#        # Rearrange ideal_shape based on optimal assignment
#        permuted_ideal = ideal_shape[col_ind]
#
#        # Compute optimal rotation using SVD
#        H = np.dot(input_structure.T, permuted_ideal)
#        U, _, Vt = svd(H)
#        R = np.dot(Vt.T, U.T)
#
#        rotated_ideal = np.dot(permuted_ideal, R)
#        scale = np.sum(input_structure * rotated_ideal) / ideal_sq_norms
#        cshm = np.mean(np.sum((input_structure - scale * rotated_ideal) ** 2, axis=1))
#        
#        min_cshm = min(min_cshm, cshm)
#
#    return min_cshm * 100

# nicely print the CShM with bars
def print_cshm_with_bars(cshm_values):
    # find the minimum and maximum cshm values
    min_cshm = min(cshm_values, key=lambda x: x[1])[1]
    max_cshm = max(cshm_values, key=lambda x: x[1])[1]
    
    # find the length of the longest label for alignment
    max_label_length = max(len(label) for label, _ in cshm_values)
    
    bar_space =  80  - (max_label_length + 12)
    
    # calculate the scale factor (highest value corresponds to remaining space to 80 chars,
    # lowest to 1 bar)
    # bar_space - 1 because 1 bar is reserved for the lowest value
    scale_factor = (bar_space - 1) / (max_cshm - min_cshm)  
    
    # print the results with the proportional bars
    for label, cshm_value in cshm_values:
        # scale the value to the range [1, to end] bars
        if max_cshm != min_cshm:
            bar_length = int((cshm_value - min_cshm) * scale_factor) + 1
        else:
            bar_length = 1  # handle the case where all values are the same
        
        bar = '░' * bar_length  # the bar is represented by '░' characters
        
        # print the label with dynamic width based on the longest label
        print(f'{label:<{max_label_length}} = {cshm_value:8.4f} {bar}')

#argument parser
parser = argparse.ArgumentParser(prog='tau-calc', 
        description = "Calculation of tau_4, tau_4', tau_5, O geometry indices and CShM.")

#filename is required
parser.add_argument('filename', 
    help = 'filename, CIF; e.g. mystructure.cif')

#atom for tau_x calculation
parser.add_argument('atom_name', 
    type = str,
    help = 'atom name; e.g. Co1')

#exclude atoms
parser.add_argument('-e','--exclude',
    nargs = '+',
    type = str,
    help = 'exclude bonds to specified atoms; e.g. -e N1 or -e N1 N2')

#exclude atoms by distance
parser.add_argument('-d','--distance',
    type = float,
    help = 'exclude atoms with distances larger than d in Å; e.g. -d 2.2')

#verbode printing
parser.add_argument('-v','--verbose',
   action = 'store_true',
     help = 'verbose print')

#save xyz coordinates of the central atom and the neighboring atoms
parser.add_argument('-sxyz','--savexyz',
   action = 'store_true',
     help = 'save the XYZ coordinates of the central atom and its neighboring atoms')

#parse arguments
args = parser.parse_args()

#load cif
try:
    doc = gemmi.cif.read_file(args.filename)
#file not found
except IOError:
    print(f"'{args.filename}'" + " not found")
    sys.exit(1)
#not a valid cif
except ValueError:
    print(f"'{args.filename}'" + " is not a valid CIF. Exit.\n")
    sys.exit(1)

#more than one or no data_block --> exit
if len(doc) != 1:
    print("CIF should contain one structure or one 'data_' block. Exit.")
    sys.exit(1)

#get the block
block = doc.sole_block()

#check if selected atom is in the CIF
if args.atom_name not in list(block.find_loop('_atom_site_label')):
    print("The atom is not part of the CIF. Exit.\n")
    sys.exit(1)

#check if excluded atoms are in the CIF
if args.exclude:
    if set(args.exclude).issubset(list(block.find_loop('_atom_site_label'))) is False:
        print("One or more excluded atoms are not part of the CIF. Exit.\n")
        sys.exit(1)

#build a bond table atom1-atom2 bond-length symmetry
bond_table = block.find(['_geom_bond_atom_site_label_1', 
                         '_geom_bond_atom_site_label_2',
                         '_geom_bond_distance',
                         '_geom_bond_site_symmetry_2'])

#delete bonds not containing the atom from input
for i in range(len(bond_table)-1, -1, -1):
    if args.atom_name not in bond_table[i]:
        del bond_table[i]

#exclude symmetry eq. atoms, e.g. for polymeric compounds M1-X-M1'
#comment in case of problems
for i in range(len(bond_table)-1, -1, -1):
    if args.atom_name in bond_table[i][1] and "." not in bond_table[i][3]:
        del bond_table[i]

#exclude bonds to specified atoms
if args.exclude:
    for atom in args.exclude:
        for i in range(len(bond_table)-1, -1, -1):
            if atom in bond_table[i]:
                del bond_table[i]
    #exit if to many excluded atoms or print excluded atoms
    if len(bond_table) == 0:
        print("Warning! To many excluded atoms. Exit.\n")
        sys.exit(1)
    else:
        print(' ')
        print("Excluded atoms: " + str(args.exclude))

#exclude bonds outside larger than a specific distance d
if args.distance:
    for i in range(len(bond_table)-1, -1, -1):
        if float(ang_bond_val.match(bond_table[i][2]).group()) > abs(args.distance):
            #add atom names of excluded atoms to list --> for angles
            if args.atom_name != bond_table[i][0]: list_of_atoms_with_large_bonds.append(bond_table[i][0]) 
            if args.atom_name != bond_table[i][1]: list_of_atoms_with_large_bonds.append(bond_table[i][1])
            del bond_table[i]
    #exit if to many excluded atoms or print excluded atoms
    if len(bond_table) == 0:
        print("Warning! To many excluded atoms. Exit.\n")
        sys.exit(1)
    #print only if the list was created
    elif list_of_atoms_with_large_bonds:
        print(' ')
        print("Excluded atoms (distance larger than " + str(abs(args.distance)) +  " Å): "
            + str(list_of_atoms_with_large_bonds))


#build an angle table atom2-atom1-atom3 angle symmetry symmetry
angle_table=block.find(['_geom_angle_atom_site_label_1',
                        '_geom_angle_atom_site_label_2',
                        '_geom_angle_atom_site_label_3',
                        '_geom_angle',
                        '_geom_angle_site_symmetry_1',
                        '_geom_angle_site_symmetry_3'])

#delete angles not containing the atom from input in the center X-M-X
for i in range(len(angle_table)-1, -1, -1):
    if args.atom_name not in angle_table[i][1].split():
        del angle_table[i]

#exclude symmetry eq. atoms, e.g. for polymeric compounds X-M1-M1'
#comment in case of problems
for i in range(len(angle_table)-1, -1, -1):
    if args.atom_name in angle_table[i][2] and "." not in angle_table[i][5]:
        del angle_table[i]

#exclude angles with specified atoms
if args.exclude:
    for atom in args.exclude:
        for i in range(len(angle_table)-1, -1, -1):
            if atom in angle_table[i]:
                del angle_table[i]

#exclude angles with specified atoms from the list of excluded atoms outside a specific distance d 
if args.distance:
    for atom in list_of_atoms_with_large_bonds:
        for i in range(len(angle_table)-1, -1, -1):
            if atom in angle_table[i]:
                del angle_table[i]

# generate list with bond lengths and sort
list_of_bonds = []
for row in bond_table:
    list_of_bonds.append(float(ang_bond_val.match(row[2]).group()))
list_of_bonds.sort()

#print bonds on request
if args.verbose:
    print(' ')
    print('--------------------------------------------------------------------------------')
    print(args.atom_name + " binds to:")
    print('--------------------------------------------------------------------------------')
    for row in bond_table:
        print(f'{row[0]}-{row[1]} {row[2]} Å {row[3]}') 

#calculate coordination number from occurance of atom in list
cnd = (list(block.find_loop('_geom_bond_atom_site_label_1')).count(args.atom_name) + 
      list(block.find_loop('_geom_bond_atom_site_label_2')).count(args.atom_name))
    
print(' ')
print(f'The predicted coordination number for {args.atom_name} is {cnd}.\n')

#exit if cn is < 3
if cnd < 3:
    print("Warning! the predicted coordination number is < 3. Exit.\n")
    sys.exit(1)
    
# generate list with angles
list_of_angles = []
for row in angle_table:
    list_of_angles.append(float(ang_bond_val.match(row[3]).group()))

# the two largest angles alpha and beta for tau_4, tau4' and tau_5
list_of_angles.sort(reverse=True)
beta = list_of_angles[0]
alpha = list_of_angles[1]

#print angles on request
if args.verbose:
    print('--------------------------------------------------------------------------------')
    print(args.atom_name + " angles are:")
    print('--------------------------------------------------------------------------------')
    for row in angle_table:
        print(f'{row[0]}-{row[1]}-{row[2]} {row[3]}° {row[4]} {row[5]}') 
    print(' ')
    
#if there is only one angle, exit
if len(list_of_angles) < 2:
    print(' ')
    print("Warning! Number of angles is too low. Exit.\n")
    sys.exit(1)

print('The two largest angles are beta = ' +
        str(beta) + '°' + ' and alpha = ' +
        str(alpha) + '°.')
print("Note: Angles for the calculation of tau_4, tau_4' and tau_5.")

# cis angles around 90 deg (< 135 deg) and trans angles around 180 deg (> 135 deg) for O
cis_ang = sum(1 for angle in list_of_angles if angle < 135)
trans_ang = sum(1 for angle in list_of_angles if angle > 135)
print(' ')
print(f'Number of cis angles ~ 90° (< 135°)    =', cis_ang)
print(f'Number of trans angles ~ 180° (> 135°) =', trans_ang)
print('Note: First value should be 12, second 3 for an octahedron.')

#calculate the cn from number of angles cn=4 number of angles = 6, 
#cn=5 number of angles=10
#cn=6 number of angles=15
cn = list(block.find_loop('_geom_angle_atom_site_label_2')).count(args.atom_name) 

printmark4 = ('')
printmark5 = ('')
printmark6 = ('')

#print '<--' arrow for most probale tau_x assignement
if cn == 6:
    printmark4=('<--')
if cn == 10:
    printmark5=('<--')
if cn == 15:
    printmark6=('<--')
    
#print warning if cn is not 4, 5 or 6
elif cn != 6 and cn != 10 and cn !=15:
    print("")
    print("Warning! The predicted coordination number seems not suitable\n" +
            "for the calculation of tau_4, tau_5 or O.")
    print("Calculated values are probably meaningless.")    

# get the XYZ coordinates for the central atom and neighbors atoms 
# gemmi neighbor search to also include symmetry equivalent positions
st = gemmi.make_small_structure_from_block(doc.sole_block())
ns = gemmi.NeighborSearch(st, 4).populate(include_h = False)

# after reaching the central atom, break
# find neighbors within min bond distances -0.01 Å and min bond distances +0.01 Å
# bond list was created above
for site in st.sites:
    if site.label == args.atom_name:
        marks = ns.find_site_neighbors(site, min_dist = list_of_bonds[0] - 0.01, 
                                             max_dist = list_of_bonds[-1] + 0.01)
        # orthogonalize the fractional coordinates of the central atom (ca)
        cart_coord_ca = st.cell.orthogonalize(site.fract)
        break

# for the calculation of the CShM values
# the coordinates must be in a numpy array
if cn == 3 or cn == 6 or cn==10 or cn == 15:
    coordinates = np.array([0.0, 0.0, 0.0]) # the central atom is at 0, 0, 0
    for mark in marks:
        # important: mark.pos gives position in unit cell, not outside
        # to_site and fract is useless in case of symmetry equivalents
        # pbc_position is the way to go 
        # check if the atom was excluded
        label = mark.to_site(st).label
        if label in list(bond_table.column(0)) or label in list(bond_table.column(1)):
            real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
            neighbor_coordinate = np.array([
                               [real_pos.x - cart_coord_ca.x, 
                                real_pos.y - cart_coord_ca.y,
                                real_pos.z - cart_coord_ca.z]
                                ])
            # add coordinates of neighbors
            coordinates = np.vstack([coordinates, neighbor_coordinate])

# calculate and print tau_x and O
print(' ')
print('--------------------------------------------------------------------------------')
print(f'{args.atom_name} geometry index ("<--" indicates the parameter for coordination number {cnd}):')
print('--------------------------------------------------------------------------------')
print(f'tau_4  = {calc_tau4(beta, alpha):6.2f} {printmark4}')
print(f"tau_4' = {calc_tau4impr(beta, alpha):6.2f} {printmark4}")
print(f'tau_5  = {calc_tau5(beta, alpha):6.2f} {printmark5}')
print(f'O      = {calc_octahedricity(list_of_angles):6.2f} {printmark6}')
# calculate and print CShM
print(' ')
print('--------------------------------------------------------------------------------')
print(f'Continuous shape measure (CShM) for coordination number {cnd}:')
print('--------------------------------------------------------------------------------')
if cn == 3 and cnd == 3:
    # calculate cshm 
    cshm_values = [
        ('S(TP-3, Trigonal planar)', calc_cshm(coordinates, IDEAL_TP3)),
        ('S(vT-3, Pyramid (vacant tetrahedron))', calc_cshm(coordinates, IDEAL_vT3)),
        ('S(fac-vOC-3, fac-Trivacant octahedron)', calc_cshm(coordinates, IDEAL_facvOC3)),
        ('S(mer-vOC-3, mer-Trivacant octahedron (T-shape))', calc_cshm(coordinates, IDEAL_mvOC3)),
    ]
    
    print_cshm_with_bars(cshm_values)
            
elif cn == 6 and cnd == 4:
    # calculate cshm
    cshm_values = [
        ('S(SP-4, Square)', calc_cshm(coordinates, IDEAL_SP4)),
        ('S(T-4, Tetrahedron)', calc_cshm(coordinates, IDEAL_T4)),
        ('S(SS-4, Seesaw or sawhorse)', calc_cshm(coordinates, IDEAL_SS4)),
        ('S(vTBPY-4, Axially vacant trigonal bipyramid)', calc_cshm(coordinates, IDEAL_vTBPY4)),
    ]
    
    print_cshm_with_bars(cshm_values)
    
elif cn == 10 and cnd == 5:
    # calculate cshm 
    cshm_values = [
        ('S(PP-5, Pentagon)', calc_cshm(coordinates, IDEAL_PP5)),
        ('S(vOC-5, Vacant octahedron (J1))', calc_cshm(coordinates, IDEAL_vOC5)),
        ('S(TBPY-5, Trigonal bipyramid)', calc_cshm(coordinates, IDEAL_TBPY5)),
        ('S(SPY-5, Square pyramid)', calc_cshm(coordinates, IDEAL_SPY5)),
        ('S(JTBPY-5, Johnson trigonal bipyramid (J12))', calc_cshm(coordinates, IDEAL_JTBPY5)),
    ]
    
    print_cshm_with_bars(cshm_values)
            
elif cn == 15 and cnd == 6:
    # calculate cshm 
    cshm_values = [
        ('S(HP-6, Hexagon)', calc_cshm(coordinates, IDEAL_HP6)),
        ('S(PPY-6, Pentagonal pyramid)', calc_cshm(coordinates, IDEAL_PPY6)),
        ('S(OC-6, Octahedron)', calc_cshm(coordinates, IDEAL_OC6)),
        ('S(TPR-6, Trigonal prism)', calc_cshm(coordinates, IDEAL_TPR6)),
        ('S(JPPY-6, Johnson pentagonal pyramid (J2))', calc_cshm(coordinates, IDEAL_JPPY6)),
    ]
    
    print_cshm_with_bars(cshm_values)

else:
    print('CShM not calculated.\n'
          'The coordination number differs from 3, 4, 5, or 6, or there is a mismatch \n'
          'between the predicted coordination number and the coordination geometry.'
          )
# calculate and print the polyhedral volume
print(' ')
print('--------------------------------------------------------------------------------')
print(f'Polyhedral volume (coordination number {cnd}) = {ConvexHull(coordinates).volume:.4f} A³')
print('--------------------------------------------------------------------------------')
#print a table of typical tau_x values
#values different from 0 or 1 and the corresponding geometries have been taken
#from an internet source - don't take it too seriously
print(' ')
print('--------------------------------------------------------------------------------')
print('Table of typical geometries and their corresponding tau_4 or tau_5 values: ')
print('--------------------------------------------------------------------------------')
print(f"Coordination number 4:")
print(f"Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00")
print(f"Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85")
print(f"Seesaw               : tau_4 = 0.43       tau_4' = 0.24")
print(f"Square planar        : tau_4 = 0.00       tau_4' = 0.00\n")
print(f"Coordination number 5:")
print(f"Trigonal bipyramidal : tau_5 = 1.00                     ")
print(f"Square pyramidal     : tau_5 = 0.00                    \n")

# print XYZ coordinates 
# set in relation to the central atom (ca) at 0, 0, 0      
if args.verbose:
    print(' ')
    print(f"--------------------------------------------------------------------------------")
    print(f"XYZ coordinates of the central atom and its neighbors: ")
    print(f"--------------------------------------------------------------------------------")
    print(f'{len(bond_table) + 1}')
    print(f'{args.filename} {args.atom_name}')
    print(f'{site.element.name:<2} {0:>11.8f} {0:>11.8f} {0:>11.8f}')
    # print orthogonalized cartesian coordinates (.pos)
    # set in relation to the central atom (ca) at 0, 0, 0
    for mark in marks:
        el_label = mark.to_site(st).element
        # important: mark.pos gives position in unit cell, not outside
        # to_site and fract is useless in case of symmetry equivalents
        # pbc_position is the way to go 
        # check if the atom was excluded
        label = mark.to_site(st).label
        if label in list(bond_table.column(0)) or label in list(bond_table.column(1)):
            real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
            print(f'{el_label.name:<2} {(real_pos.x - cart_coord_ca.x) :>11.8f} '
                  f'{(real_pos.y - cart_coord_ca.y):>11.8f} ' 
                  f'{(real_pos.z - cart_coord_ca.z):>11.8f}')

# save XYZ coordinates
# set in relation to the central atom (ca) at 0, 0, 0
if args.savexyz:
    file_name, file_extension = os.path.splitext(args.filename)
    try:
        with open(file_name + "-" + args.atom_name + '.xyz', 'w') as output_file:
            output_file.write(f'{len(bond_table) + 1}\n')
            output_file.write(f'{args.filename} {args.atom_name}\n')
            output_file.write(f'{site.element.name:<2} {0:>11.8f} {0:>11.8f} {0:>11.8f}\n')
            for mark in marks:
                el_label = mark.to_site(st).element
                # important: mark.pos gives position in unit cell, not outside
                # to_site and fract is useless in case of symmetry equivalents
                # pbc_position is the way to go 
                # check if the atom was excluded
                label = mark.to_site(st).label
                if label in list(bond_table.column(0)) or label in list(bond_table.column(1)):
                    real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
                    output_file.write(f'{el_label.name:<2} {(real_pos.x - cart_coord_ca.x):>11.8f} ' 
                                      f'{(real_pos.y - cart_coord_ca.y):>11.8f} ' 
                                      f'{(real_pos.z - cart_coord_ca.z):>11.8f}\n')
    # file not found -> exit here
    except IOError:
        print("Write error. Exit.")
        sys.exit(1)    
