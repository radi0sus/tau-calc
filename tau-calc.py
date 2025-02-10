#!/usr/bin/env python3
#
# Calculation of tau_4, tau_4', tau5 and O for 4-, 5- or 6-coordinated atoms
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
# Gemmi: 
# https://gemmi.readthedocs.io/en/latest/
# https://github.com/project-gemmi/gemmi
#

import sys                                              #sys
import os                                               #for save xyz
import argparse                                         #argument parser
import re                                               #regular expressions
import math                                             #sqrt
import gemmi                                            #CIF processing, coordinates

#regex for bonds and angles --> string to float, no esd
ang_bond_val = re.compile('\d{1,}[\.]?\d{0,}')

#list for the exclusion of atoms in the angle table
list_of_atoms_with_large_bonds=[]

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
    octahedricity = math.sqrt(sum(squared_deviations) / len(squared_deviations))
    return octahedricity

#argument parser
parser = argparse.ArgumentParser(prog='tau-calc', 
        description = "Calculation of tau_4, tau_4', tau_5 and O geometry indices.")

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
    if args.atom_name not in angle_table[i][1]:
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
    print(args.atom_name + " binds to:")
    print('------------------------------------------------------------------------')
    for row in bond_table:
        print(f'{row[0]}-{row[1]} {row[2]} Å {row[3]}') 

#calculate coordination number from occurance of atom in list
cn = (list(block.find_loop('_geom_bond_atom_site_label_1')).count(args.atom_name) + 
      list(block.find_loop('_geom_bond_atom_site_label_2')).count(args.atom_name))
    
print(' ')
print(f'The predicted coordination number for {args.atom_name} is {cn}.\n')

#exit if cn is < 3
if cn < 3:
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
    print(args.atom_name + " angles are:")
    print('------------------------------------------------------------------------')
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

#calculate and print tau_x values
print(' ')
print(args.atom_name + ' geometry indices  ("<--" indicates the likely structural parameter):')
print('------------------------------------------------------------------------')
print(f'tau_4  = {calc_tau4(beta, alpha):6.2f} {printmark4}')
print(f"tau_4' = {calc_tau4impr(beta, alpha):6.2f} {printmark4}")
print(f'tau_5  = {calc_tau5(beta, alpha):6.2f} {printmark5}')
print(f'O      = {calc_octahedricity(list_of_angles):6.2f} {printmark6}')
    
#print a table of typical tau_x values
#values different from 0 or 1 and the corresponding geometries have been taken
#from an internet source - don't take it too seriously
print(' ')
print(f"Table of typical geometries and their corresponding tau_x and O values: ")
print(f"------------------------------------------------------------------------")
print(f"Coordination number 4:")
print(f"Tetrahedral          : tau_4 = 1.00       tau_4' = 1.00")
print(f"Trigonal pyramidal   : tau_4 = 0.85       tau_4' = 0.85")
print(f"Seesaw               : tau_4 = 0.43       tau_4' = 0.24")
print(f"Square planar        : tau_4 = 0.00       tau_4' = 0.00\n")
print(f"Coordination number 5:")
print(f"Trigonal bipyramidal : tau_5 = 1.00                     ")
print(f"Square pyramidal     : tau_5 = 0.00                    \n")
print(f"Coordination number 6:")
print(f"Ideal octahedron     :     O = 0.00                    \n")

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
        
if args.verbose:
    print(f"XYZ coordinates of the central atom and its neighbors: ")
    print(f"------------------------------------------------------------------------")
    print(f'{len(marks) + 1}')
    print(f'{args.filename} {args.atom_name}')
    print(f'{site.element.name:<2} {0:>11.8f} {0:>11.8f} {0:>11.8f}')
    # print orthogonalized cartesian coordinates (.pos)
    # set in relation to the central atom (ca) at 0, 0, 0
    for mark in marks:
        label = mark.to_site(st).element
        # important: mark.pos gives position in unit cell, not outside
        # to_site and fract is useless in case of symmetry equivalents
        # pbc_position is the way to go 
        real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
        print(f'{label.name:<2} {(real_pos.x - cart_coord_ca.x) :>11.8f} '
              f'{(real_pos.y - cart_coord_ca.y):>11.8f} ' 
              f'{(real_pos.z - cart_coord_ca.z):>11.8f}')

# save XYZ coordinates
# set in relation to the central atom (ca) at 0, 0, 0
if args.savexyz:
    file_name, file_extension = os.path.splitext(args.filename)
    try:
        with open(file_name + "-" + args.atom_name + '.xyz', 'w') as output_file:
            output_file.write(f'{len(marks) + 1}\n')
            output_file.write(f'{args.filename} {args.atom_name}\n')
            output_file.write(f'{site.element.name:<2} {0:>11.8f} {0:>11.8f} {0:>11.8f}\n')
            for mark in marks:
                label = mark.to_site(st).element
                # important: mark.pos gives position in unit cell, not outside
                # to_site and fract is useless in case of symmetry equivalents
                # pbc_position is the way to go 
                real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
                output_file.write(f'{label.name:<2} {(real_pos.x - cart_coord_ca.x):>11.8f} ' 
                                  f'{(real_pos.y - cart_coord_ca.y):>11.8f} ' 
                                  f'{(real_pos.z - cart_coord_ca.z):>11.8f}\n')
    # file not found -> exit here
    except IOError:
        print("Write error. Exit.")
        sys.exit(1)    
