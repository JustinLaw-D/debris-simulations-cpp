# shell of NCell class, for file reading/writing

from CellShell import *
from ObjectsShell import *
from BreakupModel import *
import numpy as np
from copy import deepcopy
import os
import shutil
import csv

G = 6.67430e-11 # gravitational constant (N*m^2/kg^2)
Re = 6371 # radius of Earth (km)
Me = 5.97219e24 # mass of Earth (kg)

class NCell:

    def __init__(self, S, S_d, D, N_l, target_alts, alt_edges, lam, R_i=None, lam_rb=None, up_time=None, del_t=None, fail_t=None,
                 expl_rate_L=None, expl_rate_D=None, C_sat=None, sigma_sat=None, expl_rate_R=None, C_rb=None, sigma_rb=None, 
                 v=None, delta=None, alphaS=None, alphaD=None, alphaN=None, alphaR=None, P=None, m_s=None, m_rb=None, 
                 AM_sat=None, AM_rb=None, tau_do=None, L_min=1e-3, L_max=1, num_L=10, chi_min=-2, chi_max=1.5, num_chi=10,
                 update_period=1/12):
        '''
        Constructor for NCell class
    
        Parameter(s):
        S : list of initial number of live satellites in each shell of each type (list of arrays)
        S_d : list of initial number of deorbiting satellites in each shell of each type (list of arrays)
        D : list of initial number of derelict satellites in each shell of each type (list of arrays)
        N_l : initial number of catestrophically lethal debris in each shell (array)
        target_alts : list of target altitude of each satellite type (array, km)
        alt_edges : edges of the altitude bands to be used (array, km)
        lam : launch rate of satellites of each type (array, 1/yr)

        Keyword Parameter(s):
        events : the discrete events occuring in the system (list of Event objects, default no events)
        R_i : list of rocket bodies in each shell of each type (list of lists, default no rocket bodies)
        lam_rb : launch rate of rocket bodies of each type into the each shell (list of arrays, 1/yr, default all 0)
        up_time : ascention time of satellites of each type in each shell (list of arrays, yr, default all 1/10yr)
        del_t : mean satellite lifetime of each type in each shell (list of lists, yr, default 5yr)
        fail_t : ascending satellite failure lifetime (list of lists, yr, default 1000yr)
        expl_rate_L : number of explosions that occur in a 1yr period with a population of 100 live satellites for
                      each type of satellite (list of floats, default all 0)
        expl_rate_D : number of explosions that occur in a 1yr period with a population of 100 derelict satellites
                      for each type of satellite (list of floats, default all 0)
        C_sat : fit constant for explosions of each type of satellite (list of floats, default all 1)
        sigma_sat : satellite cross-section of each type (list, m^2, default 10m^2)
        expl_rate_R : number of explosions that occur in a 1yr period with a population of 100 rocket bodies for
                      each type of rocket body (list of floats, default all 0)
        C_rb : fit constant for explosions of each type of rocket body (list of floats, default all 1)
        sigma_rb : rocket cross-section of each type (list, m^2, default 10m^2)
        v : relative collision speed in each shell (list, km/s, default 10km/s)
        delta : initial ratio of the density of disabling to catestrophic debris in each shell (list, default 10)
        alphaS : fraction of collisions with another live satellite that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default 0)
        alphaD : fraction of collisions with another derelict that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default alphaN)
        alphaN : fraction of collisions with trackable debris that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default 0.2)
        alphaR : fraction of collisions with a rocket body that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default alphaN)
        P : post-mission disposal probability for satellites of each type in each shell (list of lists, default 0.95)
        m_s : mass of the satallites of each type (list, kg, default 250kg)
        m_s : mass of the rocket bodies of each type (list, kg, default 250kg)
        AM_sat : area-to-mass ratio of the satallites of each type (list, m^2/kg, default 1/(20*2.2)m^2/kg)
        AM_rb : area-to-mass ratio of the rocket bodies of each type (list, m^2/kg, default 1/(20*2.2)m^2/kg)
        tau_do : average deorbiting time for satellites of each type in each shell (list of lists, yr, default decay_time/10)
        L_min : minimum characteristic length to consider (m, default 1mm)
        L_max : maximum characteristic length to consider (m, default 1m)
        num_L : number of debris bins in characteristic length (default 10)
        chi_min : minimum log10(A/M) to consider (log10(m^2/kg), default -3)
        chi_max : maximum log10(A/M) to consider (log10(m^2/kg), default 3)
        num_chi : number of debris bins in log10(A/M) (default 10)
        update_period : how often to update the decay lifetimes (yr, default 1/12)

        Output(s):
        NCell instance

        Note: no size checks are done on the arrays, the program will crash if any of the arrays differ in size.
        shells are assumed to be given in order of ascending altitude. if you only want to pass values in the
        keyword argument for certain shells, put None in the list for all other shells. internally, cells have
        padded space in their arrays, use the getter functions to clean those up.
        '''

        # convert Nones to array of Nones
        if R_i is None:
            R_i = [[]]*len(S)
        if lam_rb is None:
            lam_rb = [None]*len(S)
        if up_time is None:
            up_time = [None]*len(S)
        if del_t is None:
            del_t = [None]*len(S)
        if fail_t is None:
            fail_t = [None]*len(S)
        if v is None:
            v = [None]*len(S)
        if delta is None:
            delta = [10]*len(S)
        if alphaS is None:
            alphaS = [None]*len(S)
        if alphaD is None:
            alphaD = [None]*len(S)
        if alphaN is None:
            alphaN = [None]*len(S)
        if alphaR is None:
            alphaR = [None]*len(S)
        if P is None:
            P = [None]*len(S)
        if tau_do is None:
            tau_do = [None]*len(S)

        self.alts = np.zeros(len(alt_edges)-1)
        self.dh = np.zeros(self.alts.shape)
        for i in range(len(alt_edges)-1):
            self.dh[i] = alt_edges[i+1]-alt_edges[i]
            self.alts[i] = (alt_edges[i]+alt_edges[i+1])/2
        self.num_L = num_L
        self.num_chi = num_chi
        self.time = 0 # index of current time step
        self.lupdate_time = 0 # index of last time drag lifetimes were updated
        self.t = [0] # list of times traversed
        self.cells = [] # start list of cells
        # generate bins for log10(L), chi
        self.logL_edges = np.linspace(np.log10(L_min), np.log10(L_max), num=num_L+1)
        self.chi_edges = np.linspace(chi_min, chi_max, num=num_chi+1)
        self.update_period = update_period

        for i in range(0, len(S)): # iterate through shells

            # convert Nones to array of Nones
            if lam_rb[i] is None:
                lam_rb[i] = [None]*len(R_i[i])
            if up_time[i] is None:
                up_time[i] = [None]*len(S[i])
            if del_t[i] is None:
                del_t[i] = [None]*len(S[i])
            if fail_t[i] is None:
                fail_t[i] = [None]*len(S[i])
            if expl_rate_L is None:
                expl_rate_L = [None]*len(S[i])
            if expl_rate_D is None:
                expl_rate_D = [None]*len(S[i])
            if C_sat is None:
                C_sat = [None]*len(S[i])
            if sigma_sat is None:
                sigma_sat = [None]*len(S[i])
            if expl_rate_R is None:
                expl_rate_R = [None]*len(R_i[i])
            if C_rb is None:
                C_rb = [None]*len(R_i[i])
            if sigma_rb is None:
                sigma_rb = [None]*len(R_i[i])
            if alphaS[i] is None:
                alphaS[i] = [None]*len(S[i])
            if alphaD[i] is None:
                alphaD[i] = [None]*len(S[i])
            if alphaN[i] is None:
                alphaN[i] = [None]*len(S[i])
            if alphaR[i] is None:
                alphaR[i] = [None]*len(S[i])
            if P[i] is None:
                P[i] = [None]*len(S[i])
            if m_s is None:
                m_s = [None]*len(S[i])
            if m_rb is None:
                m_rb = [None]*len(R_i[i])
            if AM_sat is None:
                AM_sat = [None]*len(S[i])
            if AM_rb is None:
                AM_rb = [None]*len(R_i[i])
            if tau_do[i] is None:
                tau_do[i] = [None]*len(S[i])

            sat_list = []

            for j in range(len(S[0])): # iterate through satellite types, and generate object for each
                
                # convert Nones to default values
                if up_time[i][j] is None:
                    up_time[i][j] = 1/10
                if del_t[i][j] is None:
                    del_t[i][j] = 5
                if fail_t[i][j] is None:
                    fail_t[i][j] = 5
                if expl_rate_L[j] is None:
                    expl_rate_L[j] = 0
                if expl_rate_D[j] is None:
                    expl_rate_D[j] = expl_rate_L[j]
                if C_sat[j] is None:
                    C_sat[j] = 1
                if sigma_sat[j] is None:
                    sigma_sat[j] = 10
                if alphaS[i][j] is None:
                    alphaS[i][j] = 0
                if alphaN[i][j] is None:
                    alphaN[i][j] = 0.2
                if alphaD[i][j] is None:
                    alphaD[i][j] = alphaN[i][j]
                if alphaR[i][j] is None:
                    alphaR[i][j] = alphaN[i][j]
                if P[i][j] is None:
                    P[i][j] = 0.95
                if m_s[j] is None:
                    m_s[j] = 250
                if AM_sat[j] is None:
                    AM_sat[j] = 1/(20*2.2)

                # compute atmospheric drag lifetime for satallites in the shell
                if tau_do[i][j] is None:
                    tau_do[i][j] = 0 # this value is used to communicate the value not being set
                sat = Satellite(S[i][j], S_d[i][j], D[i][j], m_s[j], sigma_sat[j], lam[j], del_t[i][j],
                                tau_do[i][j], target_alts[j], up_time[i][j], fail_t[i][j], (alphaS[i][j], alphaD[i][j],
                                alphaN[i][j], alphaR[i][j]), P[i][j], AM_sat[j], C_sat[j], expl_rate_L[j], expl_rate_D[j])
                sat_list.append(sat)

            rb_list = []

            for j in range(len(R_i[0])): # iterate through rocket types, and generate object for each
                
                # convert Nones to default values
                if lam_rb[i][j] is None:
                    lam_rb[i][j] = 0
                if expl_rate_R[j] is None:
                    expl_rate_R[j] = 0
                if C_rb[j] is None:
                    C_rb[j] = 1
                if sigma_rb[j] is None:
                    sigma_rb[j] = 10
                if m_rb[j] is None:
                    m_rb[j] = 250
                if AM_rb[j] is None:
                    AM_rb[j] = 1/(20*2.2)

                # compute atmospheric drag lifetime for rocket bodies in the shell
                rb = RocketBody(R_i[i][j], m_rb[j], sigma_rb[j], lam_rb[i][j], AM_rb[j], C_rb[j], expl_rate_R[j])
                rb_list.append(rb)

            # calculate decay paremeters for debris, initial debris values
            N_initial = np.zeros((num_L, num_chi))
            # generate initial distributions
            lethal_L = np.log10(randL(N_l[i], 1e-1, L_max, 'expl')) # explosions are the main source https://www.esa.int/esapub/bulletin/bullet109/chapter16_bul109.pdf
            nlethal_L = np.log10(randL(delta[i]*N_l[i], L_min, 1e-1, 'expl'))
            for j in range(num_L):
                bin_L = 0
                bin_bot_L, bin_top_L = self.logL_edges[j], self.logL_edges[j+1]
                bin_L += len(lethal_L[(bin_bot_L < lethal_L) & (lethal_L < bin_top_L)])
                bin_L += len(nlethal_L[(bin_bot_L < nlethal_L) & (nlethal_L < bin_top_L)])
                N_initial[j,0] = bin_L # put everything in the lowest A/M bin

            # initialize cell
            cell = Cell(sat_list, rb_list, N_initial, self.logL_edges, self.chi_edges, self.alts[i], self.dh[i], v=v[i])
            self.cells.append(cell)
            if i == len(S) - 1: self.upper_N = deepcopy(N_initial) # take the debris field above to be initial debris of top

        self.num_cells = len(self.cells)

    def save(self, filepath, name, gap=0):
        '''
        saves the current NCell object to .csv and .npy files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        name : name of the object, must be a valid unix folder name (string)

        Keyword Input(s):
        gap : largest acceptable time gap between saved data points (yr, default 0 i.e. save all data)

        Output(s): None
        Note(s): adherence to the "gap" value is approximate, and may behave strangely if the time step is 
                 close to the gap size. any previously saved data with the same name is overwritten.
        '''

        true_path = filepath + name + '/'
        try:
            os.mkdir(true_path) # make the folder representing the object
        except FileExistsError:
            shutil.rmtree(true_path)
            os.mkdir(true_path)

        # write parameters
        csv_file = open(true_path + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.num_L, self.num_chi, self.num_cells, self.update_period])
        csv_file.close()

        # write easy arrays
        t_arr = np.array(self.t)
        filter = np.full(t_arr.shape, False) # build filter based on time steps
        if t_arr.size > 0:
            prev_t = t_arr[0]
            filter[0] = True
            for i in range(1, t_arr.size):
                if t_arr[i] - prev_t >= gap:
                    prev_t = t_arr[i]
                    filter[i] = True
        np.save(true_path + "alts.npy", self.alts)
        np.save(true_path + "dhs.npy", self.dh)
        np.save(true_path + "t.npy", t_arr[filter])
        np.save(true_path + "logL.npy", self.logL_edges)
        np.save(true_path + "chi.npy", self.chi_edges)

        # save the Cells
        for i in range(self.num_cells):
            cell_path = true_path + "cell" + str(i) + "/"
            os.mkdir(cell_path)
            self.cells[i].save(cell_path, filter)

    def load(filepath):
        '''
        builds an NCell object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        atmos : NCell object build from loaded data
        '''

        atmos = NCell.__new__(NCell) # empty initialization

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, this extracts it
            atmos.num_L = int(row[0])
            atmos.num_chi = int(row[1])
            atmos.num_cells = int(row[2])
            atmos.update_period = float(row[3])
        csv_file.close()

        # load in simple numpy arrays
        atmos.alts = np.load(filepath + "alts.npy")
        atmos.dh = np.load(filepath + "dhs.npy")
        atmos.t = np.load(filepath + "t.npy").tolist()
        atmos.time = len(atmos.t) - 1 # set time to the end of the data
        atmos.lupdate_time = atmos.time
        atmos.logL_edges = np.load(filepath + "logL.npy")
        atmos.chi_edges = np.load(filepath + "chi.npy")

        # get Cells
        atmos.cells = []
        for i in range(atmos.num_cells):
            cell_path = filepath + "cell" + str(i) + "/"
            atmos.cells.append(Cell.load(cell_path))

        return atmos