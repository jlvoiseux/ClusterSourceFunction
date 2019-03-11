# -- Third-party modules

import os
import re
import subprocess
import sys
import time
from itertools import product
from math import pi
from subprocess import check_output
import numpy as np
import xlsxwriter
import paramiko
import csv
import fileinput
import threading
import matplotlib.pyplot as plt
import shutil
import copy 
import http.server
import socketserver
import os
import webbrowser

# -- class cluster

class cluster:
    """
    Cette classe permet de créer et de manipuler un cluster d'atomes. Elle utilise
    les logiciels externes Gamess-US et MultiWFN, desquels on peut extraire des données
    telles que la densité éléctronique, le laplacien de la densité, ou les bassins. Elle
    effectue aussi des calculs de Source Function et des calculs de matrice d'influence.

    """
    
    # SECTION DEDIEE À L'INITIALISATION ET À LA MODIFICATION DU CLUSTER

    periodictable = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"]

    def __init__(self, wfn_path):

        # Informations concernant la molécule
        
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        self.already_compute = [False for i in range(4)]

        self.atoms_coords = np.empty((0,3),float)
        self.atoms_charges = np.empty(0,int)
        self.atoms_names = []
        self.atoms_indexes = np.empty(0,int)
        self.atoms_nb = 0
        
        self.env_coords = np.empty((0,3),float)
        self.env_charges = np.empty(0,int)
        self.env_names = []
        self.env_indexes = np.empty(0,int)
        self.env_nb = 0

        # Informations concernant la grille en position
        self.grid_origin = [0,0,0] # Coordonnée de l'origine de la grille
        self.grid_tot_pt = 0 # Nombre de point total dans la grille
        self.grid_dir_pt = [0,0,0]  # Nombre de points selon X, Y et Z
        self.grid_incr = [0,0,0]  # Distance entre deux points selon X, Y et Z
        self.grid_elem_volume = 0  # Volume d'un cube de la grille
        self.grid_all_coords = np.zeros((0, 3))  # Liste des coordonnées (x,y,z) de tous les points de la grille

        # Informations concernant les données calculées en position
        self.density = np.array([])
        self.laplacian = np.array([])
        self.basins = np.array([])
        self.basins_indicator = np.array([])
        self.criticalpoints_coords = np.array([])
        self.criticalpoints_type = []
        self.criticalpoints_nb = 0
        self.criticalpoints_density = []

        # Informations concernant les points d'intérêt
        self.interestpoints_coords1 = np.empty((0, 4), float)
        self.interestpoints_coords2 = np.empty((0, 7), float)
        self.interestpoints_indexes = []
        self.interestpoints_names = []
        self.interestpoints_types = []
        self.interestpoints_nb = 0
        self.interestpoints_densities = []
        
        # Matrice d'influence et bassins sélectionnés
        self.AIM = []
        self.selectedbasins = []
        self.selectedbasins_coords = []
        self.selectedbasins_clusterindexes = []

        self.wfn = wfn_path

        return

    def __setitem__(self, atomindex, chargecoordinates):
        """
        Cette méthode permet d'ajouter un atome au cluster selon cluster[("A",n)] = (C, (x,y,z)). 
        Les coordonnées sont en Bohr. Mettre None à la place de l'indice si l'on ne souhaite pas préciser.
        ex : cluster[("Na",0)] = (11,(0,0,0))

        """
        
        if not(atomindex[0] in ['MP','DP','BQ']):
            (atom, index) = atomindex
            (charge, coords) = chargecoordinates
            if index == None:
                self.add_atom(atom, coords, charge)
            elif index in self.atoms_indexes[np.where(self.atoms_names == atom)[0]]:
                global_index = self.get_global_atom_index(atom, index)
                self.atoms_coords[global_index] = coords
                self.atoms_charges[global_index] = charge
            else:
                raise ValueError("Indice non conforme")
        else:
            (elm, index) = atomindex
            (charge, coords) = chargecoordinates
            if index == None:
                self.add_atom(elm, coords, charge)
            elif index in self.env_indexes[np.where(self.env_names == elm)[0]]:
                global_index = self.get_global_atom_index(elm, index)
                self.env_coords[global_index] = coords
                self.env_charges[global_index] = charge
            else:
                raise ValueError("Indice non conforme")
            
        # no longer updated
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        
        return

    def __delitem__(self, atomindex):
        """
        Cette méthode permet de supprimer un atome du cluster.

        """
        if not(atomindex[0] in ['MP','DP','BQ']): 
            (atom, index) = atomindex
            global_index = self.get_global_atom_index(atom,index)
            self.atoms_coords = np.delete(self.atoms_coords, global_index, axis = 0)
            self.atoms_charges = np.delete(self.atoms_charges, global_index)
            self.atoms_indexes = np.delete(self.atoms_indexes, global_index)
            self.atoms_names = np.delete(self.atoms_names, global_index)
            self.atoms_nb -=1
            
            # update indexes
            if index <= self.count_atom(atom):
                for i in np.where(self.atoms_names == atom)[0]:
                    if self.atoms_indexes[i] > index:
                        self.atoms_indexes[i] -= 1
        else:
            (env, index) = atomindex
            global_index = self.get_global_atom_index(env,index)
            self.env_coords = np.delete(self.env_coords, global_index, axis = 0)
            self.env_charges = np.delete(self.env_charges, global_index)
            self.env_indexes = np.delete(self.env_indexes, global_index)
            self.env_names = np.delete(self.env_names, global_index)
            self.env_nb -=1

            # update indexes
            if index <= self.count_atom(env):
                for i in np.where([self.env_names[i] == env for i in range(len(self.env_names))]):
                    if self.env_indexes[i] > index:
                        self.env_indexes[i] -= 1
        
        # no longer updated
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        
        return

    def __str__(self):
        """
        Cette méthode permet d'afficher la composition du cluster.
        ex : print(cluster)

        """
        s = ""
        for i in range(self.atoms_nb):
            s += "{}({})    {}    {}\n".format(self.atoms_names[i], self.atoms_indexes[i], self.atoms_charges[i], "    ".join(format(f, ".3f") for f in self.atoms_coords[i]))
        return s

    def __getitem__(self, atomindex):
        """
        Cette méthode renvoie un tuple de type atomindex avec la charge et les coordonnées 
        de l'atome ("A",n).

        """
        (atom, index) = atomindex
        global_index = self.get_global_atom_index(atom,index)
        
        if not(atom in ['MP','DP','BQ']): 
            return (self.atoms_charges[global_index], self.atoms_coords[global_index])
        else:
            return (self.env_charges[global_index], self.env_coords[global_index])

    def get_global_atom_index(self, atom, index):
        """
        Cette méthode renvoie l'index général du ième atome "atom".

        """
        if not(atom in ['MP','DP','BQ']): 
            return np.where(self.atoms_names == atom)[0][index]
        else:
            return np.where(self.env_names == atom)[0][index]

    def get_total_charge(self):
        """
        Cette méthode renvoie la charge totale du cluster.

        """
        return np.sum(self.atoms_charges)

    def count_atom(self,atom):
        """
        Cette méthode renvoie le nombre d'atome "atom" dans le cluster.

        """
        if not(atom in ['MP','DP','BQ']):
            try:
                return len(np.where(self.atoms_names == atom)[0])
            except:
                return 0
        
        else:
            try:
                return len(np.where(self.env_names == atom)[0])
            except:
                return 0

    def add_atom(self, atom, coordinates, charge = None):
        """
        Cette méthode ajoute un atome au cluster. Si la charge n'est pas précisée,
        elle est déterminée en fonction du symbole de l'atome et du tableau périodique.
        
        """
        if not(atom in ['MP','DP','BQ']):
            if charge == None:
                charge = self.periodictable.index(atom)+1
            self.atoms_coords = np.append(self.atoms_coords, [coordinates], axis=0)
            self.atoms_charges = np.append(self.atoms_charges, charge)
            index = self.count_atom(atom)
            self.atoms_indexes = np.append(self.atoms_indexes, index)
            self.atoms_names = np.append(self.atoms_names,atom)
            self.atoms_nb += 1
        
        else:
            self.env_coords = np.append(self.env_coords, [coordinates], axis=0)
            self.env_charges = np.append(self.env_charges, charge)
            index = self.count_atom(atom)
            self.env_indexes = np.append(self.env_indexes, index)
            self.env_names = np.append(self.env_names,atom)
            self.env_nb += 1
        
        # no longer updated
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        
        return

    def reset_composition(self):
        """
        Cette méthode supprime tous les atomes du cluster.

        """
        self.atoms_coords = np.empty((0,3),float)
        self.atoms_nb = 0
        self.atoms_charges = np.empty(0,int)
        self.atoms_names = []
        self.atoms_indexes = np.empty(0,int)
        
        # no longer updated
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        
        return

    def read_txt_atoms(self, txt_name):
        """
        Cette méthode importe les données du cluster depuis un fichier txt d'input gamess.

        """
        with open(txt_name, 'r') as txt_file:   
            for buf in txt_file.readlines():
                atom = buf.split()[0]
                charge = float(buf.split()[1])
                coords = list(map(float,buf.split()[2:5]))
                self[(atom,None)] = (charge, coords)
        self.already_compute = [False for i in range(4)]
        
        # no longer updated
        self.updated = {'input_gamess':False, 'output_gamess':False, 'input_wfn':False, 'output_wfn':False}
        
        return

    
    # SECTION DEDIEE À L'USAGE GENERAL DE GAMESS ET MULTIWFN (LOCAL + DEPORTE)
    
    def build_gamess_input(self, output_name, optimize=False, basis = "6-31G**"):
        """
        Cette méthode construit un fichier d'input pour Gamess-US.

        """
        with open(output_name + ".inp", "w") as inp:
            inp.write("! Gamess-US input file generated by our cluster building tool\n")

            # On compte la charge éléctronique totale pour déterminer la multiplicité de spin
            total_charge = self.get_total_charge()

            inp.write(" $CONTRL {} RUNTYP={} UNITS=BOHR AIMPAC=.TRUE. $END\n".format("MULT=2 SCFTYP=ROHF" \
                if total_charge%2 != 0 else "SCFTYP=RHF", "OPTIMIZE" if optimize else "ENERGY"))
            inp.write(" $SYSTEM MWORDS=5 $END\n")
            inp.write(" $SCF DIRSCF=.TRUE. $END\n")

            # Paramètres qui rendent Gamess un peu plus laxiste sur la vitesse de convergence et permettent d'augmenter
            # sa tolérance aux clusters "fantaisistes"
            inp.write(" $STATPT NSTEP=50 OPTTOL=0.0005  $END\n")

            if basis == "STO3G":
                inp.write(" $BASIS  GBASIS=STO NGAUSS=3 $END\n")
            elif basis == "6-31G**":
                inp.write(" $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1 $END\n")

            inp.write(" $GUESS  GUESS=HUCKEL $END\n")
            inp.write(" $DATA\n")
            inp.write("Title\n")
            inp.write("C1\n")

            # Encodage des atomes du cluster
            monopole, dipole = [],[]
            
            for i in range(self.atoms_nb):
                coord = self.atoms_coords[i]
                inp.write("{}      {}      {}      {}      {}\n".format(self.atoms_names[i], self.atoms_charges[i], coord[0], coord[1], coord[2]))
            
            for i in range(self.env_nb):
                if self.env_names[i] == 'BQ':
                    coord = self.env_coords[i]
                    inp.write("{}      {}      {}      {}      {}\n".format(self.env_names[i], self.env_charges[i], coord[0], coord[1], coord[2]))
                elif self.env_names[i] == 'MP':
                    monopole.append(i)
                else:
                    dipole.append(i)
                
                
            inp.write(" $END")               
                 
            
            if len(monopole)+len(dipole) != 0:
                inp.write('\n' + ' ' + '$EFRAG' + '\n' + '\n')
                for i in monopole:
                    inp.write(' ' + 'FRAGNAME= ' +self.env_names[i]+str(self.env_indexes[i]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'1  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]) + '  ' + str(self.env_coords[i][2]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'2  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]+0.1) + '  ' + str(self.env_coords[i][2]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'3  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]) + '  ' + str(self.env_coords[i][2]+0.1) + '\n')
                for i in dipole:
                    inp.write(' ' + 'FRAGNAME= ' +self.env_names[i]+str(self.env_indexes[i]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'1  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]) + '  ' + str(self.env_coords[i][2]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'2  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]+0.1) + '  ' + str(self.env_coords[i][2]) + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'3  ' + str(self.env_coords[i][0]) + '  ' + str(self.env_coords[i][1]) + '  ' + str(self.env_coords[i][2]+0.1) + '\n')
                inp.write(' $END' + '\n')
                
                for i in monopole:
                    inp.write(' '  + '$' + self.env_names[i]+str(self.env_indexes[i]) + '\n')
                    inp.write(' ' + 'Title' + '\n')
                    inp.write('COORDINATES' + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'1  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'2  ' + '0.0' + '  ' + '0.1' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'3  ' + '0.0' + '  ' + '0.0' + '  ' + '0.1' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write('STOP' + '\n')
                    inp.write('MONOPOLES' + '\n')
                    inp.write( self.env_names[i] + str(self.env_indexes[i]) +'1  ' + str(self.env_charges[i]) + '  0' +'\n')
                    inp.write('STOP' + '\n')
                    inp.write('$END')
                inp.write('\n')
                
                for i in dipole:
                    inp.write(' '  + '$' + self.env_names[i]+str(self.env_indexes[i]) + '\n')
                    inp.write(' ' + 'Title' + '\n')
                    inp.write('COORDINATES' + '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'1  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'2  ' + '0.0' + '  ' + '0.1' + '  ' + '0.0' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write(' ' + self.env_names[i] + str(self.env_indexes[i]) +'3  ' + '0.0' + '  ' + '0.0' + '  ' + '0.1' + '  ' + '0.0' + '  ' + '0.0' +  '\n')
                    inp.write('STOP' + '\n')
                    inp.write('MONOPOLES' + '\n')
                    inp.write( self.env_names[i] + str(self.env_indexes[i]) +'1  ' + self.env_charges[i][0] + '  ' + self.env_charges[i][1] + '  ' + self.env_charges[i][2] + '  ' + '\n')
                    inp.write('STOP' + '\n')
                    inp.write(' $END')
        
        
        # input_gamess updated 
        self.updated['input_gamess'] = True
        
        return

    def process_gamess_input(self, gamess_input_name, output_name):
        """
        Cette méthode lance Gamess-US avec gamess_input_name.inp et extrait un fichier 
        output (format .log) et un fichier d'input pour MultiWFN (format .wfn).

        """
        with open(gamess_input_name + ".inp", "r") as gamess_inp:
            gamess_inp.readline()

            optimize = "OPTIMIZE" in gamess_inp.readline()

            if sys.platform == "win32":
                with open("rungms.gms", "r") as paths:
                    paths.readline()
                    paths.readline()  

                    # Suppression du fichier dans le dossier restart
                    restart_folder = paths.readline().split("=")[1].rstrip("\n")
                    print(check_output("del {}\\{}.dat".format(restart_folder,gamess_input_name), shell = True).decode())

                    # Lancement de Gamess-US avec gamess_input_name.inp
                    job_path = "rungms.bat {}.inp 2016-pgi-linux-mkl 4 0 {}.log".format(gamess_input_name,output_name)
                    #subprocess.call("start cmd /K " + job_path) # Pour s"assurer que Gamess tourne bien

            elif sys.platform == "linux":
                with open("rungms", "r") as paths:
                    restart_folder = ""
                    for line in paths:
                        if line.split("=")[0] == "set GMSPATH":
                            restart_folder = line.split("=")[1].rstrip("\n")
                            break
                            #print(check_output("rm " + restart_folder + "/" + name + ".dat", shell=True).decode())

                    # Lancement de Gamess-US avec gamess_input_name.inp
                    job_path = "./rungms {}.inp 00 1 &> {}.log".format(gamess_input_name,output_name)

            p = subprocess.Popen(job_path, shell = True)
            while p.poll()==None:
                time.sleep(1)
                print("Gamess working ...")
        self.already_compute[1] = True
        
        # output_gamess updated 
        self.updated['output_gamess'] = True
        
        return

    def update_optimized_composition(self, log_name):
        """
        Cette méthode extrait les nouvelles coordonnées optimisées à partir du fichier .log de Gamess-US.

        """
        with open(log_name + ".log", "r") as log_file:
            for line in log_file:
                if "COORDINATES OF ALL ATOMS ARE (ANGS)" in line:
                    log_file.readline()
                    log_file.readline()

                    ex_atoms_nb = self.atoms_nb
                    self.reset_composition()
                    for _ in range(ex_atoms_nb):
                        line = log_file.readline().split()
                        self.add_atom(line[0], (float(line[2]), float(line[3]), float(line[4])))
                break
        return

    def build_multiwfn_input(self, gamess_input_name, output_name):
        """
        Cette méthode construit un fichier d'input pour MultiWFN.

        """
        if sys.platform == "win32":
            with open("rungms.gms", "r") as paths:
                paths.readline()
                paths.readline()
                restart_folder = paths.readline().split("=")[1].rstrip("\n")
                multiwfn_input_path = "{}\\{}.dat".format(restart_folder,gamess_input_name)

        elif sys.platform == "linux":
            with open("rungms", "r") as paths:
                restart_folder = ""
                for line in paths:
                    if line.split("=")[0] == "set GMSPATH":
                        restart_folder = line.split("=")[1].rstrip("\n")
                        break
                multiwfn_input_path = "{}/scratch/{}.dat".format(restart_folder, gamess_input_name)

        with open(output_name + ".wfn", "w") as wfn_file_out, open(multiwfn_input_path, "r") as multiwfn_input_file :
            topreached = False
            
            lines = multiwfn_input_file.readlines()
            i = 0
            while not("----- TOP" in lines[i]):
                i += 1
            print(lines[i+2])
            nb_elements = int(lines[i+2].split()[6])
            wfn_file_out.write(lines[i+1])
            wfn_file_out.write(lines[i+2])
            
            for j in range(nb_elements):
                if len(lines[i + j + 3].split()) == 10:
                    if lines[i + j + 3].split()[9] == '0.0':
                        buf = lines[i + j + 3].split()
                        wfn_file_out.write(
                            '  H    ' + buf[1] + '    ' + buf[2] + '  ' + buf[3] + '   ' + buf[4] + ' ' + buf[5] + ' ' +
                            buf[6] + '  ' + buf[7] + ' ' + buf[8] + '  1.0' + '\n')
                    else:
                        wfn_file_out.write(lines[i + j + 3])
                else:
                    wfn_file_out.write(lines[i + j + 3])
            
            for line in lines[i+nb_elements+3:]:
                if "----- END" in line:
                    break
                else:
                    wfn_file_out.write(line)                    
        
        self.already_compute[2] = True
        
        # input_wfn updated 
        self.updated['input_wfn'] = True
        
        return

    def process_multiwfn_commands(self, multiwfn_input_name, commands):
        """
        Cette méthode permets de lancer MultiWFN avec un input et une série de commandes.

        """
        def stdout_printer(p):
            for line in p.stdout:
                print(line.rstrip())

        curr_dir = os.getcwd()
        if sys.platform == "win32":
            job_path = "{} {}\\{}.wfn".format(self.wfn,curr_dir,multiwfn_input_name)
            commands.insert(0, job_path)
            p = subprocess.Popen("cmd.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        elif sys.platform == "linux":
            job_path = "{} {}/{}".format(self.wfn, curr_dir, multiwfn_input_name)
            p = subprocess.Popen(job_path, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True)

        #t = threading.Thread(target=stdout_printer, args=(p,))
        #t.start()

        for command in commands:
            p.stdin.write(command + "\n")
            print("Multiwfn working ..." + " (Last command :"+command +")")
        p.stdin.close()

        with open("wfn_exec.log", "w+") as wfn_log:
            wfn_log.write(p.stdout.read())

        #t.join()

        self.already_compute[3] = True
        
        # output_gamess updated 
        self.updated['output_wfn'] = True
        
        return

    def generate_multiwfn_commands(self, quality = ["2"], density = True, basins = True, criticalpoints = True):
        """
        Cette méthode génère la suite de commande MultiWFN (le laplacien est toujours exporté).

        """
        commands = ["5","3"] + quality + ["2","0"]
        if density:
            commands += ["5","1","8","laplacian.cub","2","0"]
        if basins:
            commands += ["17","1","1","1","9","laplacian.cub","3","0","-5","0,0","-10"]
        if criticalpoints:
            commands += ["2", "2", "3", "-4", "4","0","-10"]
        return(commands)

    def ssh_transfer_compute_old(self, address, username, password, remote_dir, output_names, gamess_deported, multiwfn_deported, sleep_time, quality = ["2"], density = True, basins = True, criticalpoints = True):
        """
        Cette méthode exporte les inputs de Gamess et de Multiwfn sur le mésocentre et y 
        lance les calculs correspondants. Si seul "gamess_deported" est vrai : export de 
        l'input gamess, lancement du calcul et téléchargement de l'output.
        Si seul "multiwfn_deported" est vrai : export de l'input multiwfn, lancement du 
        calcul et téléchargement de l'output. Si "gamess_deported" et "multiwfn_deported" 
        sont tous les deux vrais : export de l'input gamess, lancement du calcul, extraction
        du fichier .wfn sur le mésocentre, lancement du calcul multiwfn, téléchargement de l'output.

        """
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(address, username=username, password=password)
        sftp_client = ssh.open_sftp()

        if gamess_deported and multiwfn_deported:
            sftp_client.put(os.getcwd() + '\\' + output_names + '.inp', remote_dir + output_names + '.inp')
            print('python3 gamess_multiwfn_meso.py' + ' ' + output_names + ' ' + output_names + ' ' + output_names + ' ' + str(sleep_time) + ' ' + str(quality) + ' ' + str(density) + ' ' + str(basins) + ' ' + str(criticalpoints) )
            stdin, stdout, stderr = ssh.exec_command('module load python/3.5.1 ; module load intel-mpi/5.1.2; cd ' + remote_dir +'; python3 gamess_multiwfn_meso.py' + ' ' + output_names + ' ' + output_names + ' ' + output_names + ' ' + str(sleep_time) + ' ' + str(quality) + ' ' + str(density) + ' ' + str(basins) + ' ' + str(criticalpoints) )
            exit_status = stdout.channel.recv_exit_status()
            print('Gamess and Multiwfn computations done')
            if exit_status == 0:
                sftp_client.get(remote_dir + output_names + '.wfn', os.getcwd() + '\\' + output_names + '.wfn')
                print('Wavefunction downloaded')
                sftp_client.get(remote_dir + 'laplacian.cub', os.getcwd() + '\\' + 'laplacian.cub')
                print('Laplacian downloaded')
                if density:
                    sftp_client.get(remote_dir + 'density.cub', os.getcwd() + '\\' + 'density.cub')
                    print('Density downloaded')
                if basins:
                    sftp_client.get(remote_dir + 'basin.cub', os.getcwd() + '\\' + 'basin.cub')
                    print('Basins downloaded')
                if criticalpoints:
                    sftp_client.get(remote_dir + 'CPs.txt', os.getcwd() + '\\' + 'CPs.txt')
                    print('CP downloaded')
            else:
                print("Task aborted")

        elif gamess_deported:
            sftp_client.put(os.getcwd() + '\\' + output_names + '.inp', remote_dir + output_names + '.inp')
            print('python3 gamess_meso.py' + ' ' + output_names + ' ' + output_names + ' ' + str(sleep_time))
            stdin, stdout, stderr = ssh.exec_command('module load python/3.5.1 ; module load intel-mpi/5.1.2; cd ' + remote_dir +'; python3 gamess_meso.py' + ' ' + output_names + ' ' + output_names + ' ' + str(sleep_time))
            exit_status = stdout.channel.recv_exit_status()
            print('Gamess computations done')
            if exit_status == 0:
                sftp_client.get(remote_dir + output_names + '.wfn', os.getcwd() + '\\' + output_names + '.wfn')
                print('Wavefunction downloaded')
            else:
                print("Task aborted")

        elif multiwfn_deported:
            sftp_client.put(os.getcwd() + '\\' + output_names + '.wfn', remote_dir + output_names + '.wfn')
            print('python3 multiwfn_meso.py' + ' ' + output_names + ' ' + str(quality) + ' ' + str(density) + ' ' + str(basins) + ' ' + str(criticalpoints) )
            stdin, stdout, stderr = ssh.exec_command('module load python/3.5.1 ; module load intel-mpi/5.1.2; cd ' + remote_dir +'; python3 multiwfn_meso.py' + ' ' + output_names + ' ' + str(quality) + ' ' + str(density) + ' ' + str(basins) + ' ' + str(criticalpoints) )
            exit_status = stdout.channel.recv_exit_status()
            print('Multiwfn computations done')
            if exit_status == 0:
                sftp_client.get(remote_dir + 'laplacian.cub', os.getcwd() + '\\' + 'laplacian.cub')
                print('Laplacian downloaded')
                if density:
                    sftp_client.get(remote_dir + 'density.cub', os.getcwd() + '\\' + 'density.cub')
                    print('Density downloaded')
                if basins:
                    sftp_client.get(remote_dir + 'basin.cub', os.getcwd() + '\\' + 'basin.cub')
                    print('Basins downloaded')
                if criticalpoints:
                    sftp_client.get(remote_dir + 'CPs.txt', os.getcwd() + '\\' + 'CPs.txt')
                    print('CP downloaded')
            else :
                print("Task aborted")

        sftp_client.close()
        ssh.close()

    def process_gamess_wfn(self, output_names, optimize = False, basis = "STO3G", quality = ["2"], density = True, basins = True, criticalpoints = True):
        """
        Cette méthode regroupe la totalité du processus de construction des inputs et de génération des outputs
        jusqu'à l'export des données liées au cluster (densité, laplacien, bassins), 
        à partir du cluster intialement créé.

        """
        self.build_gamess_input(output_names, optimize, basis)
        self.process_gamess_input(output_names, output_names)
        self.build_multiwfn_input(output_names, output_names)
        commands = self.generate_multiwfn_commands(quality,density,basins,criticalpoints)
        self.process_multiwfn_commands(output_names, commands)
        return

    
    # SECTION DEDIEE AU CALCUL DE LA MATRICE D'INFLUENCE APPROXIMEE (ALGO MAISON)
    
    def read_cub_header(self, cub_name = "laplacian.cub"):
        """
        Cette méthode lit les premières lignes d'un fichier .cub et met à jour les données correspondantes.

        """
        with open(cub_name,"r") as cub_file:

            # skip first line and read number of points
            cub_file.readline()
            buf = cub_file.readline()
            self.grid_tot_pt = int(buf.split()[1])

            # read grid origin
            buf = cub_file.readline()
            buf = buf.split()
            self.grid_origin = [float(buf[i]) for i in range(1, 4)]

            # read direction lines
            for i in range(3):
                buf = cub_file.readline()
                buf = buf.split()
                self.grid_dir_pt[i] = int(buf[0])
                self.grid_incr[i] = float(buf[1 + i])
            self.grid_elem_volume = self.grid_incr[0] * self.grid_incr[1] * self.grid_incr[2]
        return

    def read_cub_atoms(self, cub_name = "laplacian.cub"):
        """
        Cette méthode met à jour les coordonnées des atomes selon les données du fichier .cub.

        """
        self.reset_composition()
        with open(cub_name, "r") as cub_file:
            for _ in range(2):
                cub_file.readline()
            atoms_nb = int(cub_file.readline().split()[0])

            for _ in range(3):
                cub_file.readline()

            for i in range(atoms_nb):
                buf = cub_file.readline().split()
                charge = int(buf[0])
                self.add_atom(self.periodictable[charge-1], tuple([float(buf[i]) for i in range(2, 5)]))
        return

    def skip_cub_header(self, cub_file):
        """
        Cette méthode permet de passer les premières lignes du fichier .cub ouvert pour 
        arriver directement aux données.

        """
        for i in range(6+self.atoms_nb) :
            cub_file.readline()
        return

    def get_grid_all_coords(self):
        """
        Cette méthode remplit la variable self.grid_all_coords.

        """
        x = np.linspace(self.grid_origin[0], self.grid_origin[0] + self.grid_incr[0] * self.grid_dir_pt[0], self.grid_dir_pt[0], endpoint=False)
        y = np.linspace(self.grid_origin[1], self.grid_origin[1] + self.grid_incr[1] * self.grid_dir_pt[1], self.grid_dir_pt[1], endpoint=False)
        z = np.linspace(self.grid_origin[2], self.grid_origin[2] + self.grid_incr[2] * self.grid_dir_pt[2], self.grid_dir_pt[2], endpoint=False)
        self.grid_all_coords = np.array(list(product(x, y, z)))
        return

    def read_density(self, cub_name = "density.cub"):
        """
        Cette méthode remplit la variable self.density à partir du fichier .cub correspondant.

        """
        with open(cub_name, "r") as cub_file:
            self.skip_cub_header(cub_file)
            values = []
            buf = cub_file.readline()
            while buf != "":
                values += [float(s) for s in buf.split()]
                buf = cub_file.readline()
            self.density = np.array(values)
        return

    def read_laplacian(self, cub_name = "laplacian.cub"):
        """
        Cette méthode remplit la variable self.laplacian à partir du fichier .cub correspondant.

        """
        with open(cub_name, "r") as cub_file:
            self.skip_cub_header(cub_file)
            values = []
            buf = cub_file.readline()
            while buf != "":
                values += [float(s) for s in buf.split()]
                buf = cub_file.readline()
            self.laplacian = np.array(values)
        return

    def read_basins(self, cub_name = "basin.cub"):
        """
        Cette méthode remplit la variable self.basins à partir du fichier .cub correspondant.

        """
        with open(cub_name, "r") as cub_file:
            self.skip_cub_header(cub_file)
            values = []
            buf = cub_file.readline()
            while buf != "":
                values += [int(float(s)) - 1 if int(float(s)) > 0 else -1 for s in buf.split()]
                buf = cub_file.readline()
            self.basins = np.array(values)
        return

    def rearange_basins(self):
        """
        Cette méthode réordonne les bassins pour les mettre dans l'ordre des atomes correspondants.

        """
        count = 0
        temp_basins = np.zeros(self.grid_tot_pt)-1
        for i in range(self.atoms_nb):
            distances_to_nucleus = np.sum((self.grid_all_coords - self.atoms_coords[i])**2, axis=1)
            nucleus_basin = self.basins[np.argmin(distances_to_nucleus)]
            for pt in range(self.grid_tot_pt):
                if self.basins[pt] == nucleus_basin:
                    temp_basins[pt] = count
            count += 1
        self.basins = temp_basins
        return

    def generate_basins_indicator(self):
        """
        Cette méthode génère l'indicatrice des bassins.

        """
        self.basins_indicator = np.zeros((self.atoms_nb,self.grid_tot_pt))
        for basin in range(self.atoms_nb):
            self.basins_indicator[basin] = np.int64(self.basins == basin)
        return

    def calculate_pt_sourcefunction(self, pt_coord):
        """
         /! DEPRECATED\ Cette méthode calcule les contributions source de 
         tous les bassins sur le point spécifié.

        """
        pt_sourcefunction = np.zeros(self.atoms_nb)
        distances_to_pt = np.sqrt(np.sum((self.grid_all_coords - pt_coord) ** 2, axis=1))
        distances_to_pt[distances_to_pt == 0] = 1e-20
        for basin in range(self.atoms_nb):
            pt_sourcefunction[basin] = - self.grid_elem_volume / (4 * pi) * np.sum(self.laplacian * self.basins_indicator[basin] / distances_to_pt)
        return(pt_sourcefunction)

    def read_criticalpoints(self, txt_name = "CPs.txt"):
        """
        Cette méthode permet d'importer la liste des points critiques attracteurs et de liaison.

        """
        with open(txt_name, "r") as txt_file:
            buf = txt_file.readline() # Première ligne avec le nombre de CPs
            temp_criticalpoints_coords = []
            temp_criticalpoints_type = []
            buf = txt_file.readline()
            while buf != "":
                if buf.split()[4] == "1" or buf.split()[4] == "2":
                    temp_criticalpoints_coords += [[float(s) for s in buf.split()[1:4]]]
                    temp_criticalpoints_type += ["a" if buf.split()[4] == "1" else "b"]
                    self.criticalpoints_nb += 1
                buf = txt_file.readline()
            self.criticalpoints_type = temp_criticalpoints_type
            self.criticalpoints_coords = np.array(temp_criticalpoints_coords)
        return

    def update_criticalpoints_type_attractor(self):
        """
        Cette méthode assigne les points critiques attracteurs aux noyaux les plus proches et 
        les réordonne

        """
        temp_criticalpoints_coords = np.array(self.criticalpoints_coords, copy = True)
        temp_criticalpoints_type = self.criticalpoints_type
        for cp in [i for i, x in enumerate(self.criticalpoints_type) if x == "a"]:
            distances_to_cp = np.sum((self.atoms_coords - self.criticalpoints_coords[cp])**2, axis=1)
            global_index = np.argmin(distances_to_cp)
            atom = self.atoms_names[global_index]
            index = self.atoms_indexes[global_index]
            temp_criticalpoints_coords[global_index] = self.criticalpoints_coords[cp]
            temp_criticalpoints_type[global_index] = "ACP {}({})".format(atom,index)
            print(cp,global_index,atom,index)
            #self.criticalpoints_type[cp] = "ACP {}({})".format(atom,index)
        self.criticalpoints_coords = temp_criticalpoints_coords
        self.criticalpoints_type = temp_criticalpoints_type
        return

    def calculate_distance_line_point(self,pt, linept1, linept2):
        """
        Cette méthode calcule la distance entre un point et la droite formée par deux autres points.

        """
        return(np.sum(np.cross(linept2-linept1, linept1-pt)**2)/np.sum((linept2-linept1)**2))

    def update_criticalpoints_type_bond(self):
        """
        Cette méthode assigne les points critiques de liaison aux liaisons les plus proches.

        """
        for cp in [i for i, x in enumerate(self.criticalpoints_type) if x == "b"]:
            min_distance_to_line = float("inf")
            min_distance_to_atom1, min_distance_to_atom2 = float("inf"), float("inf")
            atom1, atom2, index1, index2 = "", "", -1, -1
            for n1 in range(self.atoms_nb):
                for n2 in range(n1+1, self.atoms_nb):
                    distance_to_line = self.calculate_distance_line_point(self.criticalpoints_coords[cp], self.atoms_coords[n1], self.atoms_coords[n2])
                    distance_to_atom1 = np.sum((self.criticalpoints_coords[cp]-self.atoms_coords[n1])**2)
                    distance_to_atom2 = np.sum((self.criticalpoints_coords[cp]-self.atoms_coords[n2])**2)
                    if (distance_to_line < min_distance_to_line) or \
                       (distance_to_line < min_distance_to_line+0.5 and (distance_to_atom1 < min_distance_to_atom1 or distance_to_atom2 < min_distance_to_atom2)):
                        min_distance_to_line = distance_to_line
                        min_distance_to_atom1, min_distance_to_atom2 = distance_to_atom1, distance_to_atom2
                        atom1, index1 = self.atoms_names[n1], self.atoms_indexes[n1]
                        atom2, index2 = self.atoms_names[n2], self.atoms_indexes[n2]
            self.criticalpoints_type[cp] = "BCP {}({})-{}({})".format(atom1, index1, atom2, index2)
        return

    def calculate_criticalpoints_sourcefunction(self, percentage = True):
        """
        Cette méthode calcule la sourcefunction sur tous les points critiques du système

        """
        approx_influence_matrix = np.zeros((self.criticalpoints_nb,self.atoms_nb))
        cp_density = np.zeros(self.criticalpoints_nb)
        for cp in range(self.criticalpoints_nb):
            approx_influence_matrix[cp] = self.calculate_pt_sourcefunction(self.criticalpoints_coords[cp])
            cp_density[cp] = np.sum(approx_influence_matrix[cp])
            if percentage:
                approx_influence_matrix[cp] = approx_influence_matrix[cp]*100/cp_density[cp]
        return(approx_influence_matrix, cp_density)

    def export_approx_influence_matrix(self,output_name, percentage = True):
        """
        Cette méthode exporte la matrice d'influence approximée.

        """
        (aim, cp_density) = self.calculate_criticalpoints_sourcefunction(percentage)
        maxlen = max(map(len,self.criticalpoints_type))
        with open(output_name+".out", "w") as output_file:
            # write first line
            line = " "*maxlen + " |||"
            for i in range(self.atoms_nb):
                line += "{:^11}".format("{}({})".format(self.atoms_names[i],self.atoms_indexes[i])) + "|"
            line += "|| TOTAL DENSITY"
            output_file.write(line + "\n")

            # write second line
            line = " "*(maxlen+1) + "-"*(19+12*self.atoms_nb)
            output_file.write(line + "\n")

            # write array
            for i, row in enumerate(aim):
                line = "%s ||| %s ||| %s" % ("{:<{width}}".format(self.criticalpoints_type[i], width = maxlen),
                                                " | ".join(map("{:< 9.3g}".format, row)),
                                                "{:< 9.3g}".format(cp_density[i]))
                output_file.write(line + "\n")
        return

    def process_AIM(self, output_name):
        """
        (/!\ DEPRECATED) Cette méthode regroupe la totalité du processus de calcul 
        de la matrice d'influence approximée, et l'exporte dans le fichier output_name.out.

        """
        self.read_cub_header()
        self.read_cub_atoms()
        self.get_grid_all_coords()
        self.read_laplacian()
        self.read_basins()
        self.rearange_basins()
        self.generate_basins_indicator()
        self.read_criticalpoints()
        self.update_criticalpoints_type_attractor()
        self.update_criticalpoints_type_bond()
        self.export_approx_influence_matrix(output_name)
        return

    
    # SECTION DEDIEE AUX POINTS D'INTERET ET AU CALCUL DES CONTRIBUTIONS SOURCE A L'AIDE DE MULTIWFN

    def add_interestpoint(self, atomindex, type, atomindex2=None):
        """
        Cette méthode ajoute un point d'intérêt (où les contributions des différents bassins sont calculées) 
        au cluster.

        """
        if type == 1:
            index = self.interestpoints_nb
            self.interestpoints_nb += 1
            coords = self.__getitem__(atomindex)[1]
            self.interestpoints_coords1 = np.append(self.interestpoints_coords1,
                                                    [[index, coords[0], coords[1], coords[2]]], axis=0)
            self.interestpoints_indexes.append(index)
            self.interestpoints_names.append("ACP " + atomindex[0] + "(" + str(atomindex[1]) + ")")
            self.interestpoints_types.append(type)

        elif type == 2:
            index = self.interestpoints_nb
            self.interestpoints_nb += 1
            coords1 = self.__getitem__(atomindex)[1]
            coords2 = self.__getitem__(atomindex2)[1]
            self.interestpoints_coords2 = np.append(self.interestpoints_coords2, [
                [index, coords1[0], coords1[1], coords1[2], coords2[0], coords2[1], coords2[2]]], axis=0)
            self.interestpoints_indexes.append(index)
            self.interestpoints_names.append(
                "BCP " + atomindex[0] + "(" + str(atomindex[1]) + ") " + atomindex2[0] + "(" + str(
                    atomindex2[1]) + ") ")
            self.interestpoints_types.append(type)

        else:
            print("Please enter a valid type (1 : ACP, 2 : BCP)")

        # no longer updated
        self.updated = {'input_gamess': False, 'output_gamess': False, 'input_wfn': False, 'output_wfn': False}

        return

    def add_selectedbasins(self, atomindex):
        """
            Cette méthode sélectionne les indices des atomes que l'on demande à Multiwfn 
            de prendre en compte pour le calcul des bassins.

        """
        (atom, index) = atomindex
        glob_index = self.get_global_atom_index(atom, index)
        self.selectedbasins.append(glob_index + 1)
        self.selectedbasins_clusterindexes.append(index)
        self.selectedbasins_coords.append(self.__getitem__(atomindex)[1])

    def update_interestpoints_names(self):
        """
                Cette méthode met à jour les noms des points d'intérêt en cas d'usage des 
                groupes de symétrie de Gamess.

        """
        for i in self.interestpoints_indexes:
            if self.interestpoints_types[i] == 1:
                count = 0
                for a in self.atoms_coords:
                    if self.find_interestpoint_coords(i)[0] == a[0] and self.find_interestpoint_coords(i)[1] == a[1] and \
                                    self.find_interestpoint_coords(i)[2] == a[2]:
                        self.interestpoints_names[i] = "ACP " + str(self.atoms_names[count]) + "(" + str(
                            self.atoms_indexes[count]) + ")"
                        break
                    count += 1
            elif self.interestpoints_types[i] == 2:
                count1 = 0
                for a1 in self.atoms_coords:
                    if self.find_interestpoint_coords(i)[0] == a1[0] and self.find_interestpoint_coords(i)[1] == a1[
                        1] and self.find_interestpoint_coords(i)[2] == a1[2]:
                        temp_at1 = str(self.atoms_names[count1]) + "(" + str(self.atoms_indexes[count1]) + ")"
                        count2 = 0
                        for a2 in self.atoms_coords:
                            if self.find_interestpoint_coords(i)[3] == a2[0] and self.find_interestpoint_coords(i)[4] == \
                                    a2[1] and self.find_interestpoint_coords(i)[5] == a2[2]:
                                temp_at2 = str(self.atoms_names[count2]) + "(" + str(self.atoms_indexes[count2]) + ")"
                                self.interestpoints_names[i] = "BCP " + temp_at1 + " " + temp_at2
                                break
                            count2 += 1
                    count1 += 1

    def find_interestpoint_coords(self, index):
        """
                Cette fonction renvoie les coordonnées d'un point d'intérêt à l'aide d'un indice.

        """
        for c1 in self.interestpoints_coords1:
            if c1[0] == index:
                return c1[1:]
        for c2 in self.interestpoints_coords2:
            if c2[0] == index:
                return c2[1:]

    def get_critical_points(self, multiwfn_input_name):
        """
                Cette méthode calcule les ACP/BCP et leur densité à l'aide de Multiwfn.

        """
        commands = ["2", "2", "3", "-4", "4", "0", "7", "0", "-10"]
        self.process_multiwfn_commands(multiwfn_input_name, commands)
        self.read_criticalpoints()
        self.read_critical_points_properties()

    def read_critical_points_properties(self, density=True):
        """
                Cette méthode doit être appelée après "get_critical_points" et permet d'en extraire 
                la densité aux points critiques.

        """
        if density:
            with open("CPprop.txt", "r") as cpprop:
                for line in cpprop:
                    if "Density of all electrons" in line.split(":")[0]:
                        self.criticalpoints_density.append(float(line.split(":")[1]))

    def compute_referencepoint(self, name, type, coordinates):
        """
                Cette méthode calcule le point critique correspondant au point d'intérêt 
                dont les coordonnées sont passées en argument.

        """
        count_cp = 0
        if type == 1:
            for c in self.criticalpoints_coords:
                if coordinates[0] - 0.1 <= c[0] <= coordinates[0] + 0.1 and coordinates[1] - 0.1 <= c[1] <= coordinates[
                    1] + 0.1 and coordinates[2] - 0.1 <= c[2] <= coordinates[2] + 0.1 and self.criticalpoints_type[
                    count_cp] == 'a':
                    self.interestpoints_densities.append(self.criticalpoints_density[count_cp])
                    return c
                count_cp += 1
            print("Critical point corresponding to " + name + " not found")

        elif type == 2:
            for c in self.criticalpoints_coords:
                if self.calculate_distance_line_point(c, np.array([coordinates[0], coordinates[1], coordinates[2]]),
                                                      np.array([coordinates[3], coordinates[4],
                                                                coordinates[5]])) < 0.1 and self.criticalpoints_type[
                    count_cp] == 'b':
                    if np.linalg.norm(c - np.array([coordinates[0], coordinates[1], coordinates[2]])) + np.linalg.norm(
                                    c - np.array([coordinates[3], coordinates[4], coordinates[5]])) < np.linalg.norm(
                                    np.array([coordinates[3], coordinates[4], coordinates[5]]) - np.array(
                                    [coordinates[0], coordinates[1], coordinates[2]])) + 0.1:
                        self.interestpoints_densities.append(self.criticalpoints_density[count_cp])
                        return c
                count_cp += 1
            print("Critical point corresponding to " + name + " not found")

    def update_cluster_symmetry(self, gamess_output_name):
        """
                Cette méthode met à jour les atomes du cluster en cas d'usage des groupes de symétrie de Gamess.

        """
        with open(gamess_output_name + ".log", "r") as gamess_log:
            temp_coord = open("effective_cluster.txt", "w+")
            write_temp = False
            for line in gamess_log:
                if write_temp:
                    if line == "\n":
                        break
                    else:
                        temp_coord.write(line)
                if "CHARGE" in line:
                    write_temp = True
            temp_coord.close()
        # output_gamess updated
        self.reset_composition()
        self.read_txt_atoms("effective_cluster.txt")
        #self.update_interestpoints_names()
        #self.update_selectedbasins()
        self.updated['input_gamess'] = True
        self.updated['output_gamess'] = True
        return

    def update_selectedbasins(self):
        """
                Cette méthode met à jour le nom des bassins choisis en cas d'usage des groupes de 
                symétrie de Gamess. Elle s'assure également que tous les bassins symétriques à 
                celui choisi sont également sélectionnés afin de minimiser la perte d'info.

        """
        self.selectedbasins = []
        self.selectedbasins_clusterindexes = []
        checked_norms = []
        for c in self.selectedbasins_coords:
            count_atoms = 0
            if np.linalg.norm(c) not in checked_norms:
                for a in self.atoms_coords:
                    count_atoms += 1
                    if np.linalg.norm(a) == np.linalg.norm(c):
                        self.selectedbasins.append(count_atoms)
                        self.selectedbasins_clusterindexes.append(self.atoms_indexes[count_atoms - 1])
                checked_norms.append(np.linalg.norm(c))

    def compute_interestpoints_SF_commands(self, quality):
        """
                Cette méthode calcule les contributions source en un point d'intérêt donné, à l'aide de Multiwfn.

        """
        commands = []
        if len(self.selectedbasins) != 0:
            selectedbasins_string = ""
            for b in self.selectedbasins:
                selectedbasins_string += (str(b) + ",")
            selectedbasins_string = selectedbasins_string[:-1]
            commands += ["-3", selectedbasins_string]
        for i in self.interestpoints_indexes:
            interestpoint_type = self.interestpoints_types[i]
            interestpoint_coords = self.find_interestpoint_coords(i)
            interestpoint_name = self.interestpoints_names[i]
            ref_coords = self.compute_referencepoint(interestpoint_name, interestpoint_type, interestpoint_coords)
            if i == 0:
                grid_settings = ["1"] + quality
            else:
                grid_settings = ["2"]
            commands += ["1000", "1", str(ref_coords[0]) + "," + str(ref_coords[1]) + "," + str(ref_coords[2])] + ["17",
                                                                                                                   "1"] + grid_settings + [
                            "7", "1", "19", "-10"]
        return commands

    def extract_interestpoint_SF(self):
        """
                Cette méthode lit le .log de Multiwfn et en extrait les infos relatives à la fonction source.

        """
        contrib = []
        basins = []
        basin_name_buffer = []
        global_array = []
        count_ip = 0
        with open("wfn_exec.log", "r") as wfn_log:
            source_reached = False
            for line in wfn_log:
                if source_reached:
                    if "Sum" in line:
                        source_reached = False
                        count_ip += 1
                        global_array.append((basins, contrib, line.split()[4]))
                        if count_ip == len(self.interestpoints_indexes):
                            return global_array
                    elif "Atom" not in line:
                        count_index = 0
                        temp_atom = re.split("\(|\)", line)[1]
                        temp_atom_num = line.split()[0]
                        for i in self.atoms_indexes:
                            if count_index == int(temp_atom_num) -1:
                                basin_name_buffer.append(temp_atom)
                                basins.append(temp_atom + " (" + str(i)+ ")")
                                contrib.append(line.split(")")[1].split()[1])
                                break
                            count_index+=1
                if "Total result" in line and source_reached == False:
                    source_reached = True
                    basins = []
                    contrib = []

    def compute_AIM_multiwfn(self, multiwfn_input_name, quality):
        """
                Cette méthode calcule les contributions source en chacun des points d'intérêt.

        """
        #self.update_cluster_symmetry(multiwfn_input_name)
        self.get_critical_points(multiwfn_input_name)
        basins_written = False
        self.AIM = []
        count_ip = 0
        self.process_multiwfn_commands(multiwfn_input_name, self.compute_interestpoints_SF_commands(quality))
        data_array = self.extract_interestpoint_SF()
        for d in data_array:
            (temp_basins, temp_contrib, temp_totaldensity) = d
            if basins_written == False:
                self.AIM.append(["///"] + temp_basins + ["Total density"])
                basins_written = True
            temp_result = [self.interestpoints_names[count_ip]] + [float(f) for f in temp_contrib] + [float(temp_totaldensity)]
            temp_contrib_perc = ["{0:.8f}".format(float(100 * float(c)) / float(self.interestpoints_densities[count_ip])) for c in
                                 temp_contrib]
            temp_result_perc = [self.interestpoints_names[count_ip] + " %"] + temp_contrib_perc + [
                self.interestpoints_densities[count_ip]]
            self.AIM.append(temp_result)
            self.AIM.append(temp_result_perc)
            self.AIM.append("//////")
            print("Density computed with the sum of the source contributions : " + temp_totaldensity)
            print("Density computed with a thorough topological analysis : " + str(
                self.interestpoints_densities[count_ip]))
            count_ip += 1
        self.AIM.append(["CAUTION : If you have selected specific basins, you might suffer a non-negligible loss of information"])
        self.AIM.append(["Check the difference between the total density resulting from the sum of source contributions (top) and the total density at the critical point computed through a thorough topological analysis."])
        self.AIM.append("//////")

        if len(self.selectedbasins) != 0:
            print("Done ! Source contributions computed for " + str(count_ip) + " interest points and " + str(
                len(self.selectedbasins)) + " selected basins.")
        else:
            print("Done ! Source contributions computed for " + str(count_ip) + " interest points and every basin.")

    def export_AIM_multiwfn_CSV(self, name):
        """
                Cette méthode met en forme les contributions source au sein d'un fichier .csv

        """
        try:
            os.remove(name + ".csv")
        except OSError:
            pass
        with open("source" + ".csv", "w+", newline='') as csvfile:
            csvw = csv.writer(csvfile, dialect='excel', delimiter='|')
            for line in self.AIM:
                formatted_line = ['{:<18}' for item in line]
                s = ','.join(formatted_line)
                line = s.format(*line)
                csvw.writerow([line])
        with fileinput.FileInput("source" + ".csv", inplace=True) as file:
            for line in file:
                print(line.replace(",", ""), end='')        

    def export_AIM_multiwfn_XSLX(self, name):
        """
                Cette méthode met en forme les contributions source au sein d'un fichier .xlsx

        """
        try:
            os.remove(name + ".xlsx")
        except OSError:
            pass
        workbook = xlsxwriter.Workbook(name + ".xlsx")
        worksheet = workbook.add_worksheet()
        row = 0
        for line in self.AIM:
            col = 0
            for entry in line:
                worksheet.write(row, col, entry)
                col += 1
            row += 1
        try:
            workbook.close()
        except:
            # Handle your exception here.
            print("Couldn't create xlsx file")

    def process_AIM_multiwfn(self, output_name, quality):
        """
                L'emploi de cette méthode permet de faire un calcul complet de contributions
                source à partir d'un input Gamess.

        """
        self.process_gamess_input(output_name, output_name)
        self.build_multiwfn_input(output_name, output_name)
        self.compute_bader_charges(output_name, quality)
        self.compute_AIM_multiwfn(output_name, quality)
        self.export_AIM_multiwfn_CSV(output_name)
        
    def ssh_transfer_compute(self, address, username, password, remote_dir, output_names, core_num, mem, quality, mode):
        """
        Cette méthode lance un calcul sur le mésocentre.
        /!\ : Seuls le mode "2" et core_num = "24" fonctionnent pour l'instant.

        """
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(address, username=username, password=password)
        sftp_client = ssh.open_sftp()

        with open("job.sh", "w+") as jobfile:
            jobfile.write("#PBS -S /bin/bash\n")
            jobfile.write("#PBS -N " + output_names + "\n")
            jobfile.write("#PBS -P project name\n")
            jobfile.write("#PBS -o output.txt\n")
            jobfile.write("#PBS -e error.txt\n")
            jobfile.write("#PBS -q ivyq\n")
            jobfile.write("#PBS -l walltime=04:00:00\n")
            jobfile.write("#PBS -M user email address\n")
            jobfile.write("#PBS -m e\n")
            jobfile.write("#PBS -l select=1:ncpus=" + str(core_num) + ":mem=" + str(mem) + "gb\n")
            jobfile.write("\n")
            jobfile.write("module load intel-mpi/5.1.2 anaconda3/4.3.1\n")
            jobfile.write("\n")
            jobfile.write("# Go to the directory where the job has been submitted\n")
            jobfile.write("cd $PBS_O_WORKDIR\n")
            jobfile.write("\n")
            jobfile.write("python3 computeRemote.py " + output_names + " " + core_num + " " + quality + " " + mode + "\n")

        with open(output_names + "_interestpoints" + '.txt', "w+") as ipfile:
            count = 0
            for at in self.interestpoints_names:
                ipfile.write(at + "\n")

        with open(output_names + "_selectedbasins" + '.txt', "w+") as sbfile:
            count = 0
            for bs in self.selectedbasins_names:
                sbfile.write(bs + "\n")

        sftp_client.put(os.getcwd() + '\\' + output_names + '.txt', remote_dir + output_names + '.txt')
        sftp_client.put(os.getcwd() + '\\' + output_names + '.inp', remote_dir + output_names + '.inp')
        sftp_client.put(os.getcwd() + '\\' + 'job.sh', remote_dir + '/job.sh')
        sftp_client.put(os.getcwd() + '\\' + output_names + "_interestpoints" + '.txt', remote_dir + output_names + "_interestpoints" + '.txt')
        sftp_client.put(os.getcwd() + '\\' + output_names + "_selectedbasins" + '.txt', remote_dir + output_names + "_selectedbasins" + '.txt')

        ssh.exec_command('cd '+remote_dir +'; rm *.xlsx ; rm *.csv ; dos2unix job.sh; chmod +x job.sh; qsub job.sh')
        ready = False
        wait_time = 0
        while not ready:
            time.sleep(1)
            print("Waiting for task completion. Time elapsed : " + str(wait_time) + " seconds")
            stdin, stdout, stderr = ssh.exec_command('cd '+remote_dir +'; if ( -f' + output_names + '.xlsx ) then echo Exists else if ( -f' + output_names + '.csv ) then echo Exists else echo Nope endif')
            try:
                filestat = sftp_client.stat(remote_dir + output_names + ".xlsx")
            except Exception:
                pass
            else:
                ready = True
            wait_time += 1

        sftp_client.get(remote_dir + output_names + '.wfn', os.getcwd() + "\\" + output_names + '.wfn')
        sftp_client.get(remote_dir + output_names + '.log', os.getcwd() + "\\" + output_names + '.log')
        sftp_client.get(remote_dir + output_names + '.xlsx', os.getcwd() + "\\" + output_names + '.xlsx')
        sftp_client.get(remote_dir + output_names + '.csv', os.getcwd() + "\\" + output_names + '.csv')
        sftp_client.close()
        ssh.close()
       
    def visualize(self):
        """
        Lance le module de visualisation pour le cluster courant
        /!\ À n'utiliser qu'après avoir calculé une première fonction source.
        
        """
        
        with open("cluster_viz.txt", "w+") as vizfile:
            vizfile.write(self.__str__())
        PORT = 8084
        Handler = http.server.SimpleHTTPRequestHandler
        httpd = socketserver.TCPServer(("", PORT), Handler)
        print("serving at port", PORT)
        webbrowser.open("http://localhost:8084/", 1, True)
        httpd.serve_forever()
        
    def compute_bader_charges(self, name, quality):
        """
        Méthode de calcul de charges de Bader
        Nécessite la création préalable d'un .wfn, et une qualité entre 1 et 4.
        """

        commands = ["17", "1", "1"] +  quality + ["2", "0"]
        self.process_multiwfn_commands(name, commands)
        basins_coords = []
        basins_charge = []
        with open("wfn_exec.log", "r") as wfn_log:
            coords_reached = False
            charge_reached = False
            bohr_radius = 0.52917721092

            for line in wfn_log:

                if "Detecting boundary grids..." in line:
                    coords_reached = False
                elif coords_reached:
                    basins_coords.append([str(float(line.split()[1])/bohr_radius),str(float(line.split()[2])/bohr_radius), str(float(line.split()[3])/bohr_radius)])
                elif "Attractor       X,Y,Z coordinate (Angstrom)" in line:
                    coords_reached = True

                if "Sum of above values" in line:
                    charge_reached = False
                elif charge_reached:
                    basins_charge.append(line.split()[1])
                elif "#Basin        Integral(a.u.)      Volume(a.u.^3)" in line:
                    charge_reached = True

        final_bas = [] # Liste des charges de bader, réordonnées pour coller à la liste d'atomes du cluster
        for at in self.atoms_coords:
            count_b = 0
            for bas in basins_coords:
                if abs(at[0] - float(bas[0])) < 0.5 and abs(at[1] - float(bas[1])) < 0.5 and abs(at[2] - float(bas[2])) < 0.5:
                    final_bas.append(basins_charge[count_b])
                    break
                count_b += 1

        with open(name + "_bader_charges.txt", "w+") as outputfile:
            count = 0
            for bas in final_bas:
                outputfile.write(self.atoms_names[count] +"("+ str(self.atoms_indexes[count]) +")"+ " " + bas +"\n")
                count +=1

    def mulliken_charge(self, nb):
        path1 = "mc_temp"
        path2 = "mc_temp"

        try:
            os.remove(path1)
            os.remove(path2)
        except:
            print("no files to remove")

        self.build_gamess_input(path1)
        self.process_gamess_input(path1, path2)

        file = open(path2, 'r')
        buf = file.readline().split()
        while len(buf) < 2 or buf[0] != "TOTAL" or buf[1] != "MULLIKEN":
            buf = file.readline().split()
        file.readline()
        for i in range(nb):
            file.readline()
        buf = file.readline().split()

        return float(buf[2])

    def lowdin_charge(self, nb):
        path1 = "mc_temp"
        path2 = "mc_temp"

        try:
            os.remove(path1)
            os.remove(path2)
        except:
            print("no files to remove")

        self.build_gamess_input(path1)
        self.process_gamess_input(path1, path2)

        file = open(path2 + ".log", 'r')
        buf = file.readline().split()
        while len(buf) < 2 or buf[0] != "TOTAL" or buf[1] != "MULLIKEN":
            buf = file.readline().split()
        file.readline()
        for i in range(nb):
            file.readline()
        buf = file.readline().split()

        return float(buf[4])




