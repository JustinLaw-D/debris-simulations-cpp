# shell of Object classes in Python, for reading/writing purposes

import numpy as np
import csv

class Satellite:

    def __init__(self, S_i, S_di, D_i, m, sigma, lam, del_t, tau_do, target_alt, up_time, fail_t, 
                 alpha, P, AM, tau, C, expl_rate_L, expl_rate_D):
        '''
        constructor method for Satellite class

        Parameter(s):
        S_i : initial number of live satellites of this type
        S_di : initial number of de-orbiting satellites of this type
        D_i : initial number of derelict satellites of this type
        m : mass of each satellite (kg)
        sigma : collision cross-section of each satellite (m^2)
        lam : launch rate of the satellites (1/yr)
        del_t : mean satellite lifetime (yr)
        fail_t : ascending satellite failure lifetime (yr)
        tau_do : mean time for satellite to de-orbit from shell (yr)
        target_alt : target final altitude for the satellite type (km)
        up_time : amount of time it takes a satellite to ascend through the band (yr)
        alpha : tuple of (alphaS, alphaD, alphaN, alphaR), which are the fraction of collisions a
                live satellite fails to avoid with a live satellite, derelict, trackable debris,
                and rocket body respectively
        P : post-mission disposal probability
        AM : area-to-mass ratio of the satellite (m^2/kg)
        tau : atmospheric drag lifetime of a satellite (yr)
        C : fit constant for explosions
        expl_rate_L : number of explosions that occur in a 1yr period with a population 
                      of 100 live satellites
        expl_rate_D : number of explosions that occur in a 1yr period with a population 
                      of 100 derelict satellites
                      
        Keyword Parameter(s): None

        Output(s): Instance of Satellite class

        Note(s): preforms no validity checks on given values
        '''

        self.S = [S_i]
        self.S_d = [S_di]
        self.D = [D_i]
        self.m = m
        self.sigma = sigma
        self.lam = lam
        self.del_t = del_t
        self.fail_t = fail_t
        self.tau_do = tau_do
        self.target_alt = target_alt
        self.up_time = up_time
        self.alphaS, self.alphaD, self.alphaN, self.alphaR = alpha
        self.P = P
        self.AM = AM
        self.tau = tau
        self.C = C
        self.expl_rate_L = expl_rate_L
        self.expl_rate_D = expl_rate_D

    def save(self, filepath, filter):
        '''
        saves the current Satellite object to .csv and .npy files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)
        Keyword Input(s):
        compress : whether or not to save the data in a compressed format (default True)
        Output(s): None
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.m, self.sigma, self.lam, self.del_t, self.fail_t, self.tau_do, self.target_alt, 
                             self.up_time, self.alphaS, self.alphaD, self.alphaN, self.alphaR, self.P, 
                             self.AM, self.tau, self.C, self.expl_rate_L, self.expl_rate_D])
        csv_file.close()

        # save data
        S_array, Sd_array, D_array = np.array(self.S)[filter], np.array(self.S_d)[filter], np.array(self.D)[filter]
        np.save(filepath + "dataS.npy", S_array)
        np.save(filepath + "dataSd.npy", Sd_array)
        np.save(filepath + "dataD.npy", D_array)

    def load(filepath):
        '''
        builds a Satellite object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        sat : Satellite object build from loaded data
        '''

        sat = Satellite.__new__(Satellite) # make a blank satellite

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, but this extracts it
            sat.m, sat.sigma, sat.lam, sat.del_t = float(row[0]), float(row[1]), float(row[2]), float(row[3])
            sat.fail_t, sat.tau_do, sat.target_alt, sat.up_time = float(row[4]), float(row[5]), float(row[6]), float(row[7])
            sat.alphaS, sat.alphaD, sat.alphaN, sat.alphaR = float(row[8]), float(row[9]), float(row[10]), float(row[11])
            sat.P, sat.AM, sat.tau, sat.C = float(row[12]), float(row[13]), float(row[14]), float(row[15])
            sat.expl_rate_L, sat.expl_rate_D = float(row[16]), float(row[17])
        csv_file.close()

        # load data
        sat.S = np.load(filepath + 'dataS.npy').tolist()
        sat.D = np.load(filepath + 'dataD.npy').tolist()
        sat.S_d = np.load(filepath + 'dataSd.npy').tolist()

        return sat

class RocketBody:

    def __init__(self, num, m, sigma, lam, AM, tau, C, expl_rate):
        '''
        constructor method for RocketBody class

        Parameter(s):
        num : initial number of rocket bodies of this type
        m : mass of each rocket body (kg)
        sigma : collision cross-section of each rocket body (m^2)
        lam : launch rate of the rocket bodies (1/yr)
        AM : area-to-mass ratio of a rocket body (m^2/kg)
        tau : atmospheric drag lifetime of a rocket body (yr)
        C : fit constant for explosions
        expl_rate : number of explosions that occur in a 1yr period with a population 
                    of 100 rocket bodies

        Keyword Parameter(s): None

        Output(s): Instance of Rocket class

        Note(s): preforms no validity checks on given values
        '''

        self.num = [num]
        self.m = m
        self.sigma = sigma
        self.lam = lam
        self.AM = AM
        self.tau = tau
        self.C = C
        self.expl_rate = expl_rate

    def save(self, filepath, filter):
        '''
        saves the current Rocket object to .csv and .npy files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)

        Keyword Input(s): None

        Output(s): None
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.m, self.sigma, self.lam, self.AM, self.tau, self.C, self.expl_rate])
        csv_file.close()

        # save data
        num_array = np.array(self.num)[filter]
        np.save(filepath + "data.npy", num_array)

    def load(filepath):
        '''
        builds a RocketBody object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        rb : RocketBody object build from loaded data
        '''

        rb = RocketBody.__new__(RocketBody) # creates empty instance

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # only one row, but this extracts it
            rb.m, rb.sigma, rb.lam, rb.AM = float(row[0]), float(row[1]), float(row[2]), float(row[3])
            rb.tau, rb.C, rb.expl_rate = float(row[4]), float(row[5]), float(row[6])
        csv_file.close()

        # load data
        rb.num = np.load(filepath + "data.npy").tolist()

        return rb    