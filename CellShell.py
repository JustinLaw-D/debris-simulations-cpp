# shell of Cell class in Python, for reading/writing purposes

import numpy as np
from ObjectsShell import *
import os
import csv

class Cell:
    
    def __init__(self, S_i, R_i, N_i, logL_edges, chi_edges, alt, dh, v=None):
        '''Constructor for Cell class
    
        Parameter(s):
        S_i : list of satellite types with initial values
        R_i : list of rocket body types with initial values
        N_i : initial array of number of debris by L and A/M
        logL_edges : bin edges in log10 of characteristic length (log10(m))
        chi_edges : bin edges in log10(A/M) (log10(m^2/kg))
        alt : altitude of the shell centre (km)
        dh : width of the shell (km)
        
        Keyword Parameter(s):
        v : relative collision speed (km/s, default 10km/s)

        Output(s):
        Cell instance
        '''

        # set default values as needed
        if v == None:
            v = 10

        # setup initial values for tracking live satallites, derelict satallites,
        # lethat debris, and non-lethal debris over time
        self.satellites = S_i
        self.rockets = R_i
        self.num_sat_types = len(self.satellites)
        self.num_rb_types = len(self.rockets)
        self.N_bins = [N_i]

        # setup other variables
        self.C_l = [0] # lethal collisions
        self.C_nl = [0] # non-lethal collisions
        self.alt = alt
        self.dh = dh
        self.v = v
        self.logL_edges = logL_edges
        self.chi_edges = chi_edges
        self.num_L = self.N_bins[0].shape[0]
        self.num_chi = self.N_bins[0].shape[1]

    def save(self, filepath, filter):
        '''
        saves the current Cell object to .csv and .npy files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)

        Keyword Input(s): None

        Output(s): None

        Note(s): filter should be the same size as the t array from NCell.
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.num_sat_types, self.num_rb_types, self.alt, self.dh, self.v, self.num_L, self.num_chi])
        csv_file.close()

        # write easy arrays
        Cl_array, Cnl_array = np.array(self.C_l)[filter], np.array(self.C_nl)[filter]
        np.save(filepath + "dataCl.npy", Cl_array)
        np.save(filepath + "dataCnl.npy", Cnl_array)
        np.save(filepath + "logL.npy", self.logL_edges)
        np.save(filepath + "chi.npy", self.chi_edges)

        # write N_bins values
        bin_path = filepath + "N_bins/"
        os.mkdir(bin_path)
        index = 0
        for i in range(len(self.N_bins)):
            if filter[i]:
                np.save(bin_path + str(index), self.N_bins[i])
                index += 1

        # write satellites and rockets
        for i in range(self.num_sat_types):
            sat_path = filepath + 'Satellite' + str(i) + '/'
            os.mkdir(sat_path)
            self.satellites[i].save(sat_path, filter)
        for i in range(self.num_rb_types):
            rb_path = filepath + 'RocketBody' + str(i) + '/'
            os.mkdir(rb_path)
            self.rockets[i].save(rb_path, filter)

    def load(filepath):
        '''
        builds a Cell object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        cell : Cell object build from loaded data
        '''

        cell = Cell.__new__(Cell) # create blank Cell

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, this extracts it
            cell.num_sat_types = int(row[0])
            cell.num_rb_types = int(row[1])
            cell.alt = float(row[2])
            cell.dh = float(row[3])
            cell.v = float(row[4])
            cell.num_L = int(row[5])
            cell.num_chi = int(row[6])
        csv_file.close()

        # load basic arrays
        cell.C_l = np.load(filepath + "dataCl.npy")
        cell.C_nl = np.load(filepath + "dataCnl.npy")
        cell.logL_edges = np.load(filepath + "logL.npy")
        cell.chi_edges = np.load(filepath + "chi.npy")

        # load N_bins values
        cell.N_bins = []
        bin_path = filepath + "N_bins"
        i = 0
        while True:
            try:
                N_bin = np.load(bin_path + str(i))
                cell.N_bins.append(N_bin)
            except OSError:
                break
            i += 1

        # load satellites and rockets
        cell.satellites = []
        cell.rockets = []
        for i in range(cell.num_sat_types):
            sat_path = filepath + 'Satellite' + str(i) + '/'
            cell.satellites.append(Satellite.load(sat_path))
        for i in range(cell.num_rb_types):
            rb_path = filepath + 'RocketBody' + str(i) + '/'
            cell.rockets.append(RocketBody.load(rb_path))

        return cell
