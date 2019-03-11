# -- Third-party modules

import os
import subprocess
import sys
import time
from itertools import product
from math import pi
from subprocess import check_output
import numpy as np
import xlsxwriter
import paramiko
import threading
import matplotlib.pyplot as plt
import copy as copy
import cluster

class point_charges:
    '''classe qui étudie la disposition des chages ponctuelles dans un cluster'''
    
    def __init__(self, cluster_new):
        
        self.cluster_1 = cluster_new
        self.name = "analysor_test"
        return


    def optimal_ponctual_charge(self, atomindex):
        atom, index = atomindex
        nb1 = self.cluster_1.get_global_atom_index(atom, index)
        pop = self.cluster_1.lowdin_charge(nb1) - self.cluster_1[atomindex][0]
        coords = self.cluster_1[atomindex][1]

        project_temp = cluster.cluster(self.cluster_1.wfn)

        project_temp.atoms_coords = self.cluster_1.atoms_coords
        project_temp.atoms_charges = self.cluster_1.atoms_charges
        project_temp.atoms_names = self.cluster_1.atoms_names
        project_temp.atoms_indexes = self.cluster_1.atoms_indexes
        project_temp.atoms_nb = self.cluster_1.atoms_nb

        project_temp['BQ', None] = (-self.cluster_1[atomindex][0], coords)
        del project_temp[atomindex]

        g = [pop * (1 / 3)]
        d = [pop * (2 / 3)]

        project_temp['MP', None] = (g[0], coords)
        nb2 = len(np.where(self.cluster_1.env_names == 'BQ')[0])
        g.append(pop - project_temp.lowdin_charge(nb1 + nb2))
        index_mp = len(np.where(self.cluster_1.env_names == 'MP')[0])
        print(index_mp)
        project_temp['MP', index_mp] = (d[0], coords)
        d.append(pop - project_temp.lowdin_charge(nb2))
        a = (d[1] - g[1]) / (d[0] - g[0])
        b = d[1] - d[0] * a
        print(nb1)
        print(coords)
        print(pop)
        print(nb2)
        print(-b / a)
        return -b / a


    
