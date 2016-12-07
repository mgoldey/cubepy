import numpy as np

#####################################################################
#                                                                   #
#      Constants taken from:  http://physics.nist.gov/constants     #
#                           (1/14/2014)                             #
#                                                                   #
#   Atomic masses taken from http://www.nist.gov/pml/data/comp.cfm  #
#                            (2/8/2014)                             #
#                                                                   #
#####################################################################

def info():
    print("Selection of constants and conversion factors")
    print("---------------------------------------------")
    print("")
    print("Fundamental constants are given in SI units. Factors for conversion from x to y are abbreviated as x_to_y.")
    print("Constants are taken from http://physics.nist.gov/constants (1/14/2014)")

# Fundamentals in SI units                                     Value                 Uncertainty           Unit
#---------------------------------------------------------------------------------------------------------------------------
Angstrom                                        =           1.00000000e-10      #  (exact)                  m
atomic_mass_constant                            =           1.660538921e-27     #  0.000 000 073 e-27       kg
atomic_mass_constant_energy_equivalent          =           1.492417954e-10     #  0.000 000 066 e-10       J
atomic_unit_of_1st_hyperpolarizability          =           3.206361449e-53     #  0.000 000 071 e-53       C^3 m^3 J^-2
atomic_unit_of_2nd_hyperpolarizability          =           6.23538054e-65      #  0.000 000 28 e-65        C^4 m^4 J^-3
atomic_unit_of_action                           =           1.054571726e-34     #  0.000 000 047 e-34       J s
atomic_unit_of_charge                           =           1.602176565e-19     #  0.000 000 035 e-19       C
atomic_unit_of_charge_density                   =           1.081202338e12      #  0.000 000 024 e12        C m^-3
atomic_unit_of_electric_dipole_moment           =           8.47835326e-30      #  0.000 000 19 e-30        C m
atomic_unit_of_electric_field                   =           5.14220652e11       #  0.000 000 11 e11         V m^-1
atomic_unit_of_electric_field_gradient          =           9.71736200e21       #  0.000 000 21 e21         V m^-2
atomic_unit_of_electric_polarizability          =           1.6487772754e-41    #  0.000 000 0016 e-41      C^2 m^2 J^-1
atomic_unit_of_electric_potential               =           27.21138505         #  0.000 000 60             V
atomic_unit_of_electricvquadrupole_moment       =           4.486551331e-40     #  0.000 000 099 e-40       C m^2
atomic_unit_of_energy                           =           4.35974434e-18      #  0.000 000 19 e-18        J
atomic_unit_of_force                            =           8.23872278e-8       #  0.000 000 36 e-8         N
atomic_unit_of_length                           =           0.52917721092e-10   #  0.000 000 000 17 e-10    m
atomic_unit_of_mag_dipole_moment                =           1.854801936e-23     #  0.000 000 041 e-23       J T^-1
atomic_unit_of_mag_flux_density                 =           2.350517464e5       #  0.000 000 052 e5         T
atomic_unit_of_magnetizability                  =           7.891036607e-29     #  0.000 000 013 e-29       J T^-2
atomic_unit_of_mass                             =           9.10938291e-31      #  0.000 000 40 e-31        kg
atomic_unit_of_momentum                         =           1.992851740e-24     #  0.000 000 088 e-24       kg m s^-1
atomic_unit_of_permittivity                     =           1.112650056e-10     #  (exact)                  F m^-1
atomic_unit_of_time                             =           2.418884326502e-17  #  0.000 000 000 012 e-17   s
atomic_unit_of_velocity                         =           2.18769126379e6     #  0.000 000 000 71 e6      m s^-1
Avogadro_constant                               =           6.02214129e23       #  0.000 000 27 e23         mol^-1
Bohr_magneton                                   =           927.400968e-26      #  0.000 020 e-26           J T^-1
Bohr_radius                                     =           0.52917721092e-10   #  0.000 000 000 17 e-10    m
Boltzmann_constant                              =           1.3806488e-23       #  0.000 0013 e-23          J K^-1
classical_electron_radius                       =           2.8179403267e-15    #  0.000 000 0027 e-15      m
Compton_wavelength                              =           2.4263102389e-12    #  0.000 000 0016 e-12      m
electric_constant                               =           8.854187817e-12     #  (exact)                  F m^-1
electron_charge_to_mass_quotient                =           -1.758820088e11     #  0.000 000 039 e11        C kg^-1
electron_g_factor                               =           -2.00231930436153   #  0.000 000 000 000 53
electron_mass                                   =           9.10938291e-31      #  0.000 000 40 e-31        kg
electron_volt                                   =           1.602176565e-19     #  0.000 000 035 e-19       J
elementary_charge                               =           1.602176565e-19     #  0.000 000 035 e-19       C
Fermi_coupling_constant                         =           1.166364e-5         #  0.000 005 e-5            GeV^-2
fine_structure_constant                         =           7.2973525698e-3     #  0.000 000 0024 e-3
lattice_parameter_of_silicon                    =           543.1020504e-12     #  0.000 0089 e-12          m
Loschmidt_constant_1_bar                        =           2.6516462e25        #  0.000 0024 e25           m^-3
Loschmidt_constant_1_atm                        =           2.6867805e25        #  0.000 0024 e25           m^-3
magnetic_constant                               =           12.566370614e-7     #  (exact)                  N A^-2
molar_gas_constant                              =           8.3144621           #  0.000 0075               J mol^-1 K^-1
molar_Planck_constant                           =           3.9903127176e-10    #  0.000 000 0028 e-10      J s mol^-1
Newtonian_constant_of_gravitation               =           6.67384e-11         #  0.000 80 e-11            m^3 kg^-1 s^-2
nuclear_magneton                                =           5.05078353e-27      #  0.000 000 11 e-27        J T^-1
Planck_constant                                 =           6.62606957e-34      #  0.000 000 29 e-34        J s
Planck_constant_over_2_pi                       =           1.054571726e-34     #  0.000 000 047 e-34       J s
Rydberg_constant                                =           10973731.568539     #  0.000 055                m^-1
speed_of_light_in_vacuum                        =           299792458           #  (exact)                  m s^-1
standard_acceleration_of_gravity                =           9.80665             #  (exact)                  m s^-2
standard_atmosphere                             =           101325              #  (exact)                  Pa
standard_state_pressure                         =           100000              #  (exact)                  Pa
Stefan_Boltzmann_constant                       =           5.670373e-8         #  0.000 021 e-8            W m^-2 K^-4
Thomson_cross_section                           =           0.6652458734e-28    #  0.000 000 0013 e-28      m^2
Wien_frequency_displacement_law_constant        =           5.8789254e10        #  0.000 0053 e10           Hz K^-1
Wien_wavelength_displacement_law_constant       =           2.8977721e-3        #  0.000 0026 e-3           m K



#  Conversion factors                  Value                 Uncertainty           Final Unit
#---------------------------------------------------------------------------------------------------------------------------
eV_to_amu                      =         1.073544150e-9      #  0.000 000 024 e-9        u
eV_to_hartree                  =         3.674932379e-2      #  0.000 000 081 e-2        E_h
eV_to_hertz                    =         2.417989348e14      #  0.000 000 053 e14        Hz
eV_to_inverse_m                =         8.06554429e5        #  0.000 000 18 e5          m^-1
eV_to_inverse_cm               =         8.06554429e3        #  0.000 000 18 e3          cm^-1
eV_to_joule                    =         1.602176565e-19     #  0.000 000 035 e-19       J
eV_to_kelvin                   =         1.1604519e4         #  0.000 0011 e4            K
eV_to_kilogram                 =         1.782661845e-36     #  0.000 000 039 e-36       kg
hartree_to_amu                 =         2.9212623246e-8     #  0.000 000 0021 e-8       u
hartree_to_eV                  =        27.21138505          #  0.000 000 60             eV
hartree_to_hertz               =         6.579683920729e15   #  0.000 000 000 033 e15    Hz
hartree_to_inverse_m           =         2.194746313708e7    #  0.000 000 000 011 e7     m^-1
hartree_to_inverse_cm          =         2.194746313708e5    #  0.000 000 000 011 e5     cm^-1
hartree_to_joule               =         4.35974434e-18      #  0.000 000 19 e-18        J
hartree_to_kelvin              =         3.1577504e5         #  0.000 0029 e5            K
hartree_to_kilogram            =         4.85086979e-35      #  0.000 000 21 e-35        kg
hertz_to_amu                   =         4.4398216689e-24    #  0.000 000 0031 e-24      u
hertz_to_eV                    =         4.135667516e-15     #  0.000 000 091 e-15       eV
hertz_to_hartree               =         1.5198298460045e-16 #  0.000 000 000 0076 e-16  E_h
hertz_to_inverse_m             =         3.335640951e-9      #  (exact)                  m^-1
hertz_to_inverse_cm            =         3.335640951e-11     #  (exact)                  cm^-1
hertz_to_joule                 =         6.62606957e-34      #  0.000 000 29 e-34        J
hertz_to_kelvin                =         4.7992434e-11       #  0.000 0044 e-11          K
hertz_to_kilogram              =         7.37249668e-51      #  0.000 000 33 e-51        kg
joule_to_amu                   =         6.70053585e9        #  0.000 000 30 e9          u
joule_to_eV                    =         6.24150934e18       #  0.000 000 14 e18         eV
joule_to_hartree               =         2.29371248e17       #  0.000 000 10 e17         E_h
joule_to_hertz                 =         1.509190311e33      #  0.000 000 067 e33        Hz
joule_to_inverse_m             =         5.03411701e24       #  0.000 000 22 e24         m^-1
joule_to_inverse_cm            =         5.03411701e22       #  0.000 000 22 e22         cm^-1
joule_to_kelvin                =         7.2429716e22        #  0.000 0066 e22           K
joule_to_kilogram              =         1.112650056e-17     #  (exact)                  kg
kelvin_to_amu                  =         9.2510868e-14       #  0.000 0084 e-14          u
kelvin_to_eV                   =         8.6173324e-5        #  0.000 0078 e-5           eV
kelvin_to_hartree              =         3.1668114e-6        #  0.000 0029 e-6           E_h
kelvin_to_hertz                =         2.0836618e10        #  0.000 0019 e10           Hz
kelvin_to_inverse_m            =        69.503476            #  0.00000063               m^-1
kelvin_to_inverse_cm           =         0.69503476          #  0.0000000063             cm^-1
kelvin_to_joule                =         1.3806488e-23       #  0.000 0013 e-23          J
kelvin_to_kilogram             =         1.5361790e-40       #  0.000 0014 e-40          kg
kilogram_to_amu                =         6.02214129e26       #  0.000 000 27 e26         u
kilogram_to_eV                 =         5.60958885e35       #  0.000 000 12 e35         eV
kilogram_to_hartree            =         2.061485968e34      #  0.000 000 091 e34        E_h
kilogram_to_hertz              =         1.356392608e50      #  0.000 000 060 e50        Hz
kilogram_to_inverse_m          =         4.52443873e41       #  0.000 000 20 e41         m^-1
kilogram_to_inverse_cm         =         4.52443873e39       #  0.000 000 20 e39         m^-1
kilogram_to_joule              =         8.987551787e16      #  (exact)                  J
kilogram_to_kelvin             =         6.5096582e39        #  0.000 0059 e39           K



# Conversions added to NIST selection
#---------------------------------------------------------------------------
bohr_to_angstrom                =         0.52917721092
angstrom_to_bohr                =         1.889726125

joule_to_cal                    =         0.2390057361
cal_to_joule                    =         4.184
joule_to_kcal                   =         0.2390057361e-3
kcal_to_joule                   =         4184

inverse_cm_to_hartree           =         1/2.194746313708e5
inverse_cm_to_hertz             =         1/3.335640951e-11
inverse_cm_to_joule             =         1/5.03411701e22
inverse_cm_to_kelvin            =         1/0.69503476
inverse_cm_to_kilogram          =         1/4.52443873e39
inverse_cm_to_eV                =         1/8.06554429e3


hartree_to_kcal_pro_mole        =         hartree_to_joule*joule_to_kcal*Avogadro_constant
hartree_to_kJ_pro_mole          =         hartree_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_hartree        =         1/hartree_to_kcal_pro_mole
kJ_pro_mole_to_hartree          =         1/hartree_to_kJ_pro_mole


kelvin_to_kcal_pro_mole         =         kelvin_to_joule*joule_to_kcal*Avogadro_constant
kelvin_to_kJ_pro_mole           =         kelvin_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_kelvin         =         1/kelvin_to_kcal_pro_mole
kJ_pro_mole_to_kelvin           =         1/kelvin_to_kJ_pro_mole

hertz_to_kcal_pro_mole          =         hertz_to_joule*joule_to_kcal*Avogadro_constant
hertz_to_kJ_pro_mole            =         hertz_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_hertz          =         1/hertz_to_kcal_pro_mole
kJ_pro_mole_to_hertz            =         1/hertz_to_kJ_pro_mole

inverse_cm_to_kcal_pro_mole     =         inverse_cm_to_joule*joule_to_kcal*Avogadro_constant
inverse_cm_to_kJ_pro_mole       =         inverse_cm_to_joule*Avogadro_constant/1000
kcal_pro_mole_to_inverse_cm     =         1/inverse_cm_to_kcal_pro_mole
kJ_pro_mole_to_inverse_cm       =         1/inverse_cm_to_kJ_pro_mole

atomic_unit_of_time_to_picosec  =        1e-12/atomic_unit_of_time
atomic_unit_of_time_to_femtosec =        1e-15/atomic_unit_of_time
atomic_unit_of_time_to_attosec  =        1e-18/atomic_unit_of_time

picosec_to_atomic_unit_of_time  =        1/atomic_unit_of_time_to_picosec
femtosec_to_atomic_unit_of_time =        1/atomic_unit_of_time_to_femtosec
attosec_to_atomic_unit_of_time  =        1/atomic_unit_of_time_to_attosec

picosec_to_femtosec             =        1e3
picosec_to_attosec              =        1e6
femtosec_to_picosec             =        1e-3
femtosec_to_attosec             =        1e3
attosec_to_femtosec             =        1e-3
attosec_to_picosec              =        1e-6



# Relative atomic mass per most common isotope
#---------------------------------------------
# Using http://www.nist.gov/pml/data/comp.cfm

dict_of_atomic_masses = {'Ru': 95.907598, 'Yb': 167.933897, 'Db': 268.12545, 'Re': 184.952955, 'Rf': 265.1167, 'Ra': 223.0185022, 'Rb': 84.911789738, '90': 230.0331338, 'Rn': 210.990601, 'Rh': 102.905504, '24': 49.9460442, 'Be': 9.0121822, '26': 53.9396105, '27': 58.933195, '20': 39.96259098, '21': 44.9559119, '22': 45.9526316, '23': 49.9471585, 'Sg': 271.13347, '28': 57.9353429, '29': 62.9295975, 'Bk': 247.070307, '4': 9.0121822, 'Br': 78.9183371, '8': 15.99491461956, '96': 243.0613891, '87': 223.0197359, 'H': 1.00782503207, 'P': 30.97376163, '15': 30.97376163, 'Os': 183.9524891, '91': 231.035884, 'Es': 252.08298, '59': 140.9076528, '58': 135.907172, '55': 132.905451933, '54': 123.905893, '57': 137.907112, 'Cl': 34.96885268, '51': 120.9038157, '50': 111.904818, '53': 126.904473, '52': 119.90402, 'Ge': 69.9242474, 'Gd': 151.919791, 'Ga': 68.9255736, 'Pr': 140.9076528, '56': 129.9063208, 'Pt': 189.959932, 'Pu': 238.0495599, 'C': 12.0, 'Pb': 203.9730436, 'Np': 236.04657, '88': 223.0185022, '89': 227.0277521, 'Pd': 101.905609, '82': 203.9730436, '83': 208.9803987, '80': 195.965833, 'Xe': 123.905893, '86': 210.990601, 'Po': 208.9824304, '84': 208.9824304, 'Pm': 144.912749, '108': 270.13465, 'Hs': 270.13465, '3': 6.015122795, 'Ho': 164.9303221, '7': 14.0030740048, 'Hf': 173.940046, 'K': 38.96370668, 'He': 4.00260325415, 'Md': 258.098431, 'Mg': 23.9850417, '25': 54.9380451, 'Mo': 91.906811, 'Mn': 54.9380451, 'O': 15.99491461956, 'Mt': 276.15116, 'S': 31.972071, '109': 276.15116, 'W': 179.946704, '102': 259.10103, '103': 262.10963, '100': 257.095105, '101': 258.098431, '106': 271.13347, '107': 272.13803, '104': 265.1167, 'Ba': 129.9063208, '39': 88.9058483, '38': 83.913425, '33': 74.9215965, '32': 69.9242474, '31': 68.9255736, '30': 63.9291422, '37': 84.911789738, '36': 77.9203648, '35': 78.9183371, '34': 73.9224764, 'Eu': 150.9198502, 'U': 233.0396352, 'Zr': 89.9047044, 'Ni': 57.9353429, 'No': 259.10103, 'Na': 22.9897692809, 'Nb': 92.9063781, 'Nd': 141.9077233, 'Ne': 19.9924401754, 'Bi': 208.9803987, '60': 141.9077233, '61': 144.912749, '62': 143.911999, '63': 150.9198502, '64': 151.919791, '65': 158.9253468, '66': 155.924283, '67': 164.9303221, '68': 161.928778, '69': 168.9342133, 'Bh': 272.13803, 'Fr': 223.0197359, '2': 4.00260325415, 'Fe': 53.9396105, '6': 12.0, '105': 268.12545, 'Pa': 231.035884, 'Fm': 257.095105, 'B': 10.012937, '97': 247.070307, 'F': 18.99840322, 'Sr': 83.913425, 'Zn': 63.9291422, 'N': 14.0030740048, '99': 252.08298, 'Kr': 77.9203648, 'Si': 27.9769265325, 'Sn': 111.904818, 'Sm': 143.911999, 'V': 49.9471585, 'Sc': 44.9559119, 'Sb': 120.9038157, '93': 236.04657, '92': 233.0396352, '95': 241.0568291, '94': 238.0495599, 'Se': 73.9224764, 'Hg': 195.965833, '11': 22.9897692809, '10': 19.9924401754, '13': 26.98153863, '12': 23.9850417, 'Co': 58.933195, '14': 27.9769265325, '17': 34.96885268, '16': 31.972071, '19': 38.96370668, '18': 35.967545106, 'Ca': 39.96259098, 'Cf': 249.0748535, 'Ce': 135.907172, 'Cd': 105.906459, 'Lu': 174.9407718, 'Cs': 132.905451933, 'Cr': 49.9460442, 'Cu': 62.9295975, 'La': 137.907112, 'Li': 6.015122795, 'Tl': 202.9723442, 'Tm': 168.9342133, 'Lr': 262.10963, 'Th': 230.0331338, 'Ti': 45.9526316, 'Te': 119.90402, 'Tb': 158.9253468, 'Tc': 96.906365, 'Ta': 179.9474648, '81': 202.9723442, '48': 105.906459, '49': 112.904058, '46': 101.905609, '47': 106.905097, '44': 95.907598, '45': 102.905504, '42': 91.906811, '43': 96.906365, '40': 89.9047044, '41': 92.9063781, '1': 1.00782503207, '5': 10.012937, 'Dy': 155.924283, '9': 18.99840322, '85': 209.987148, 'I': 126.904473, 'Cm': 243.0613891, '77': 190.960594, '76': 183.9524891, '75': 184.952955, '74': 179.946704, '73': 179.9474648, '72': 175.9414086, '71': 174.9407718, '70': 167.933897, 'Y': 88.9058483, '79': 196.9665687, '78': 189.959932, 'Ac': 227.0277521, 'Ag': 106.905097, '98': 249.0748535, 'Ir': 190.960594, 'Am': 241.0568291, 'Al': 26.98153863, 'As': 74.9215965, 'Ar': 35.967545106, 'Au': 196.9665687, 'At': 209.987148, 'In': 112.904058}
for i in range(109):
    dict_of_atomic_masses[i+1]=dict_of_atomic_masses[str(i+1)]

dict_of_atomic_abbr = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Uub', 113: 'Uut', 114: 'Uuq', 115: 'Uup', 116: 'Uuh', 117: 'Uus', 118: 'Uuo'}
for i in range(118):
    dict_of_atomic_abbr[str(i+1)]=dict_of_atomic_abbr[i+1]

dict_of_atomic_numbers= {'Ru': 44, 'Re': 75, 'Rf': 104, 'Rg': 111, 'Ra': 88, 'Rb': 37, 'Rn': 86, 'Rh': 45, 'Be': 4, 'Ba': 56, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'Uuh': 116, 'H': 1, 'P': 15, 'Os': 76, 'Es': 99, 'Hg': 80, 'Ge': 32, 'Gd': 64, 'Ga': 31, 'Uub': 112, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'C': 6, 'Pb': 82, 'Pa': 91, 'Pd': 46, 'Cd': 48, 'Po': 84, 'Pm': 61, 'Hs': 108, 'Uuq': 114, 'Uup': 115, 'Uus': 117, 'Uuo': 118, 'Ho': 67, 'Hf': 72, 'K': 19, 'He': 2, 'Md': 101, 'Mg': 12, 'Mo': 42, 'Mn': 25, 'O': 8, 'Mt': 109, 'S': 16, 'W': 74, 'Zn': 30, 'Eu': 63, 'Zr': 40, 'Er': 68, 'Ni': 28, 'No': 102, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10, 'Np': 93, 'Fr': 87, 'Fe': 26, 'Fm': 100, 'B': 5, 'F': 9, 'Sr': 38, 'N': 7, 'Kr': 36, 'Si': 14, 'Sn': 50, 'Sm': 62, 'V': 23, 'Sc': 21, 'Sb': 51, 'Sg': 106, 'Se': 34, 'Co': 27, 'Cm': 96, 'Cl': 17, 'Ca': 20, 'Cf': 98, 'Ce': 58, 'Xe': 54, 'Lu': 71, 'Cs': 55, 'Cr': 24, 'Cu': 29, 'La': 57, 'Li': 3, 'Tl': 81, 'Tm': 69, 'Lr': 103, 'Th': 90, 'Ti': 22, 'Te': 52, 'Tb': 65, 'Tc': 43, 'Ta': 73, 'Yb': 70, 'Db': 105, 'Dy': 66, 'Ds': 110, 'I': 53, 'U': 92, 'Y': 39, 'Ac': 89, 'Ag': 47, 'Uut': 113, 'Ir': 77, 'Am': 95, 'Al': 13, 'As': 33, 'Ar': 18, 'Au': 79, 'At': 85, 'In': 49}

dict_of_atomic_names = {1: 'Hydrogen', 2: 'Helium', 3: 'Lithium', 4: 'Beryllium', 5: 'Boron', 6: 'Carbon', 7: 'Nitrogen', 8: 'Oxygen', 9: 'Fluorine', 10: 'Neon', 11: 'Sodium', 12: 'Magnesium', 13: 'Aluminum', 14: 'Silicon', 15: 'Phosphorus', 16: 'Sulfur', 17: 'Chlorine', 18: 'Argon', 19: 'Potassium', 20: 'Calcium', 21: 'Scandium', 22: 'Titanium', 23: 'Vanadium', 24: 'Chromium', 25: 'Manganese', 26: 'Iron', 27: 'Cobalt', 28: 'Nickel', 29: 'Copper', 30: 'Zinc', 31: 'Gallium', 32: 'Germanium', 33: 'Arsenic', 34: 'Selenium', 35: 'Bromine', 36: 'Krypton', 37: 'Rubidium', 38: 'Strontium', 39: 'Yttrium', 40: 'Zirconium', 41: 'Niobium', 42: 'Molybdenum', 43: 'Technetium', 44: 'Ruthenium', 45: 'Rhodium', 46: 'Palladium', 47: 'Silver', 48: 'Cadmium', 49: 'Indium', 50: 'Tin', 51: 'Antimony', 52: 'Tellurium', 53: 'Iodine', 54: 'Xenon', 55: 'Cesium', 56: 'Barium', 57: 'Lanthanum', 58: 'Cerium', 59: 'Praseodymium', 60: 'Neodymium', 61: 'Promethium', 62: 'Samarium', 63: 'Europium', 64: 'Gadolinium', 65: 'Terbium', 66: 'Dysprosium', 67: 'Holmium', 68: 'Erbium', 69: 'Thulium', 70: 'Ytterbium', 71: 'Lutetium', 72: 'Hafnium', 73: 'Tantalum', 74: 'Tungsten', 75: 'Rhenium', 76: 'Osmium', 77: 'Iridium', 78: 'Platinum', 79: 'Gold', 80: 'Mercury', 81: 'Thallium', 82: 'Lead', 83: 'Bismuth', 84: 'Polonium', 85: 'Astatine', 86: 'Radon', 87: 'Francium', 88: 'Radium', 89: 'Actinium', 90: 'Thorium', 91: 'Protactinium', 92: 'Uranium', 93: 'Neptunium', 94: 'Plutonium', 95: 'Americium', 96: 'Curium', 97: 'Berkelium', 98: 'Californium', 99: 'Einsteinium', 100: 'Fermium', 101: 'Mendelevium', 102: 'Nobelium', 103: 'Lawrencium', 104: 'Rutherfordium', 105: 'Dubnium', 106: 'Seaborgium', 107: 'Bohrium', 108: 'Hassium', 109: 'Meitnerium', 110: 'Darmstadtium', 111: 'Roentgenium', 112: 'Ununbium', 113: 'Ununtrium', 114: 'Ununquadium', 115: 'Ununpentium', 116: 'Ununhexium', 117: 'Ununseptium', 118: 'Ununoctium'}
for i in range(118):
    dict_of_atomic_names[str(i+1)]=dict_of_atomic_names[i+1]

dict_abbr_to_name={'Ru': 'Ruthenium', 'Re': 'Rhenium', 'Rf': 'Rutherfordium', 'Rg': 'Roentgenium', 'Ra': 'Radium', 'Rb': 'Rubidium', 'Rn': 'Radon', 'Rh': 'Rhodium', 'Be': 'Beryllium', 'Ba': 'Barium', 'Bh': 'Bohrium', 'Bi': 'Bismuth', 'Bk': 'Berkelium', 'Br': 'Bromine', 'Uuh': 'Ununhexium', 'H': 'Hydrogen', 'P': 'Phosphorus', 'Os': 'Osmium', 'Es': 'Einsteinium', 'Hg': 'Mercury', 'Ge': 'Germanium', 'Gd': 'Gadolinium', 'Ga': 'Gallium', 'Uub': 'Ununbium', 'Pr': 'Praseodymium', 'Pt': 'Platinum', 'Pu': 'Plutonium', 'C': 'Carbon', 'Pb': 'Lead', 'Pa': 'Protactinium', 'Pd': 'Palladium', 'Cd': 'Cadmium', 'Po': 'Polonium', 'Pm': 'Promethium', 'Hs': 'Hassium', 'Uuq': 'Ununquadium', 'Uup': 'Ununpentium', 'Uus': 'Ununseptium', 'Ho': 'Holmium', 'Hf': 'Hafnium', 'K': 'Potassium', 'He': 'Helium', 'Md': 'Mendelevium', 'Mg': 'Magnesium', 'Mo': 'Molybdenum', 'Mn': 'Manganese', 'O': 'Oxygen', 'Mt': 'Meitnerium', 'S': 'Sulfur', 'W': 'Tungsten', 'Zn': 'Zinc', 'Eu': 'Europium', 'Zr': 'Zirconium', 'Er': 'Erbium', 'Ni': 'Nickel', 'No': 'Nobelium', 'Na': 'Sodium', 'Nb': 'Niobium', 'Nd': 'Neodymium', 'Ne': 'Neon', 'Np': 'Neptunium', 'Fr': 'Francium', 'Fe': 'Iron', 'Fm': 'Fermium', 'B': 'Boron', 'F': 'Fluorine', 'Sr': 'Strontium', 'N': 'Nitrogen', 'Kr': 'Krypton', 'Si': 'Silicon', 'Sn': 'Tin', 'Sm': 'Samarium', 'V': 'Vanadium', 'Sc': 'Scandium', 'Sb': 'Antimony', 'Sg': 'Seaborgium', 'Se': 'Selenium', 'Co': 'Cobalt', 'Cm': 'Curium', 'Cl': 'Chlorine', 'Ca': 'Calcium', 'Cf': 'Californium', 'Ce': 'Cerium', 'Xe': 'Xenon', 'Lu': 'Lutetium', 'Cs': 'Cesium', 'Cr': 'Chromium', 'Cu': 'Copper', 'La': 'Lanthanum', 'Li': 'Lithium', 'Tl': 'Thallium', 'Tm': 'Thulium', 'Lr': 'Lawrencium', 'Th': 'Thorium', 'Ti': 'Titanium', 'Te': 'Tellurium', 'Tb': 'Terbium', 'Tc': 'Technetium', 'Ta': 'Tantalum', 'Yb': 'Ytterbium', 'Db': 'Dubnium', 'Dy': 'Dysprosium', 'Ds': 'Darmstadtium', 'I': 'Iodine', 'U': 'Uranium', 'Y': 'Yttrium', 'Ac': 'Actinium', 'Ag': 'Silver', 'Uut': 'Ununtrium', 'Ir': 'Iridium', 'Am': 'Americium', 'Al': 'Aluminum', 'As': 'Arsenic', 'Ar': 'Argon', 'Au': 'Gold', 'At': 'Astatine', 'In': 'Indium'}


#from http://en.wikipedia.org/wiki/Covalent_radius
covalent_radii=np.array([0.31,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,
                0.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,
                2.03,1.76,1.7,1.6,1.53,1.39,1.61,1.52,1.5,
                1.24,1.32,1.22,1.22,1.2,1.19,1.2,1.2,1.16,
                2.2,1.95,1.9,1.75,1.64,1.54,1.47,1.46,1.42,
                1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.4,
                2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,
                1.96,1.94,1.92,1.92,1.89,1.9,1.87,1.87,1.75,
                1.7,1.62,1.51,1.44,1.41,1.36,1.36,1.32,1.45,
                1.46,1.48,1.4,1.5,1.5,2.6,2.21,2.15,2.06,
                2.,1.96,1.9,1.87,1.8,1.69])


#from http://toc.uni-muenster.de/DFTD3/index.html, times 1.1, then converted to Angstroms
vdw_radii=np.array([1.892,1.912,1.559,2.661,2.806,2.744,2.640,2.536,2.432,2.349,
           2.162,2.578,3.097,3.243,3.222,3.180,3.097,3.014,2.806,2.785,
           2.952,2.952,2.952,2.952,2.952,2.952,2.952,2.952,2.952,2.952,
           3.118,3.264,3.326,3.347,3.305,3.264,3.076,3.035,3.097,3.097,
           3.097,3.097,3.097,3.097,3.097,3.097,3.097,3.097,3.160,3.409,
           3.555,3.575,3.575,3.555,3.405,3.330,3.251,3.313,3.313,3.313,
           3.313,3.313,3.313,3.313,3.313,3.313,3.313,3.313,3.313,3.313,
           3.313,3.378,3.349,3.349,3.349,3.349,3.349,3.349,3.349,3.322,
           3.752,3.673,3.586,3.789,3.762,3.636])


# from http://life.nthu.edu.tw/~fmhsu/rasframe/CPKCLRS.HTM
atomic_colors={1: (1.0, 1.0, 1.0), 6: (0.0,0.0,0.0), 7: (0.5607843137254902, 0.5607843137254902, 1.0), 8: (0.9411764705882353, 0.0, 0.0), 11: (0.0, 0.0, 1.0), 12: (0.16470588235294117, 0.5019607843137255, 0.16470588235294117), 15: (1.0, 0.6470588235294118, 0.0), 16: (1.0, 0.7843137254901961, 0.19607843137254902), 17: (0.0, 1.0, 0.0), 20: (0.5019607843137255, 0.5019607843137255, 0.5019607843137255), 26: (1.0, 0.6470588235294118, 0.0), 30: (0.6470588235294118, 0.16470588235294117, 0.16470588235294117), 35: (0.6470588235294118, 0.16470588235294117, 0.16470588235294117)}

 # orig 6: (0.7843137254901961, 0.7843137254901961, 0.7843137254901961),
#Unknown  [255,20,147]   
