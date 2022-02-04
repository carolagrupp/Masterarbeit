# Sources:
# ------------------------------------------------------------------------------
# Imports: 
from bisect import bisect_left, bisect_right
import csv
import bisect
from json import loads
import sys
import hdf5storage as h5s
import numpy as np
from pfh import building
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Abbreviations
# ------------------------------------------------------------------------------
# p...  load
# r...  response
# ms... model scale
# fs... full scale
# L...  lift
# D...  drag
# M...  moment
# F...  force
# H...  (at) height of building
# sp..  sample
# fq... frequency
# dn... direction
# ------------------------------------------------------------------------------
class DataProp():
    def __init__(self):
        #Eigenschaften der Daten
        self.Data = []
        self.scale_ms = []
        self.frequency_ms = [] #Hz, s. Tanaka et all
        
    def calcModelFloorForces_F(self, buildingProp, loads):
        """Calculating floor forces from wind tunnel model base forces in model scale
        with an assumed distribution after alpha_v
        :param buildingProp: full scale properties of the building
        :type buildProp: :class
        """
        alpha_v = loads.alpha_v
        

    def readData_F(self,fname, buildingProp):
        """Loading wind tunnel model characteristics for force measurement data
            """
        # Load matlab file
        self.fname = fname
        mat = h5s.loadmat(self.fname)

        # Force measurements
        self.cfd = mat['CFD']
        self.cfl = mat['CFL']
        self.cfx = mat['CFX']
        self.cfy = mat['CFY']
        self.cfz = mat['CFZ']
        self.cmd = mat['CMD']
        self.cml = mat['CML']
        self.cmx = mat['CMX']
        self.cmy = mat['CMY']
        self.cmz = mat['CMZ']
        self.vH = mat['Vh']                                                                             # =Velocity auf Höhe H über Zeit

        self.vH_ms = np.mean(self.vH)
        #oder self.vH_ms = v_fs /10 , da Velocity Scale in Datenblatt mit 1/10 angegeben??

        # Anzahl der genommenen Datenpunkten
        self.NumberOfData = float(mat['NumberOfData'][0])

        # PowerLawIndex = alpha
        self.PowerLawIndex = float(mat['PowerLawIndex'][0])

        # WindDirection of current data
        self.WindDirection = float(mat['WindDirection'][0])


        #self.Datenreihen = [self.cfd,self.cfl,self.cfx,self.cfy,self.cfz,self.cmd,self.cml,self.cmx,self.cmy,self.cmz,self.vh_ms]

        self.H_ms = 0.400 #in m
        self.B_ms = 0.100

        # Time scales
        self.fq_sp = 200 #in Hz
        self.dT = 1 / self.fq_sp
        #self.T = Sample Period
        #self.nT = int(self.T * self.fq_sp)

        # cf = F/(vBH)
        self.BFD_ms = []
        self.BFL_ms = []
        for i in range(0,len(self.cfd)):
            self.BFD_ms.append(self.cfd[i]*self.vH_ms*self.H_ms*self.B_ms)       #CFD ist Coefficient, jetzt * vH_ms, H und B
            self.BFL_ms.append(self.cfl[i]*self.vH_ms*self.H_ms*self.B_ms)

        #cm = M/(vBH**2)
        self.BMD_ms = []
        self.BML_ms = []
        for i in range(0,len(self.cmd)):
            self.BMD_ms.append(self.cmd[i]*self.vH_ms*self.H_ms**2*self.B_ms)
            self.BML_ms.append(self.cml[i]*self.vH_ms*self.H_ms**2*self.B_ms)

        #with open("Masterarbeit/WindData/01_Square_alpha_015_0.CSV") as csvdatei:
            #csv_reader_object = csv.reader(csvdatei, delimiter=';')
            #lines = csvdatei.readlines()

    def calcModelBaseForces_P(self):
        """Calculating base forces from wind tunnel model forces in model scale
        """
        # Sort in forces
        self.BFD_ms = np.zeros(self.nT)
        self.BFL_ms = np.zeros(self.nT)
        self.BMD_ms = np.zeros(self.nT)
        self.BML_ms = np.zeros(self.nT)

        for i, F_i in enumerate(self.F_p.T, 0):
            # Face 1 -> +DragF, +DragM
            if self.face[i] == 1: #and buildProp.dn == 'D':
                self.BFD_ms = self.BFD_ms + F_i
                self.BMD_ms = self.BMD_ms + F_i * self.z[i] 
            # Face 2 -> +LiftF, -LiftM
            elif self.face[i] == 2: # and buildProp.dn == 'L':
                self.BFL_ms = self.BFL_ms + F_i
                self.BML_ms = self.BML_ms - F_i * self.z[i] 
            # Face 3 -> -DragF, -DragM
            elif self.face[i] == 3: # and buildProp.dn == 'D':
                self.BFD_ms = self.BFD_ms - F_i
                self.BMD_ms = self.BMD_ms - F_i * self.z[i] 
            # Face 4 -> -LiftF, +LiftM
            elif self.face[i] == 4:# and buildProp.dn == 'L':
                self.BFL_ms = self.BFL_ms - F_i
                self.BML_ms = self.BML_ms + F_i * self.z[i]
    
    def calcModelFloorForces_P(self):
        """Calculating floor forces from wind tunnel model forces in model scale
        """
        # Sort in drag / lift forces for each floor
        self.LFD_ms = np.zeros((self.nz, self.nT))
        self.LFL_ms = np.zeros((self.nz, self.nT))

        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p.T, 0):
            # get index where to sort to
            j = np.where(self.z_lev == self.z[i])
            # Face 1 -> +DragF, +DragM
            if self.face[i] == 1:# and buildProp.dn == 'D':
                self.LFD_ms[j,:] = self.LFD_ms[j,:] + F_i
            # Face 2 -> +LiftF, -LiftM
            elif self.face[i] == 2:# and buildProp.dn == 'L':
                self.LFL_ms[j,:] = self.LFL_ms[j,:] + F_i
            # Face 3 -> -DragF, -DragM
            elif self.face[i] == 3:# and buildProp.dn == 'D':
                self.LFD_ms[j,:] = self.LFD_ms[j,:] - F_i
           # Face 4 -> -LiftF, +LiftM
            elif self.face[i] == 4:# and buildProp.dn == 'L':
                self.LFL_ms[j,:] = self.LFL_ms[j,:] - F_i


    def readData_P(self, fname):
        """Loading wind tunnel model characteristics for pressure coefficient data
            """
        # Load matlab file
        self.fname = fname
        mat = h5s.loadmat(self.fname)

        # Geometrie
        self.H_ms = float(mat['Building_height'][0])

        # Measurement locations
        self.Loc_ms = mat['Location_of_measured_points']

        # Coordinates [in m], number and face of measurement points
        self.x = self.Loc_ms[0]
        self.z = self.Loc_ms[1]     #Koordinate in die Höhe (von unten gemessen und oben beginnend)
        self.n = self.Loc_ms[2]
        self.face = self.Loc_ms[3]

        # Measurement area [in m2]
        self.A_i = 0.02 * 0.02  #Wieso 0.02?

        # Air density
        self.rho_air = 1.25

        # Get levels & sort them
        self.z_lev = np.unique(self.z)              # sortiert von klein nach groß
        self.z_lev = self.z_lev[::-1]               # dreht z_lev wieder um -> wieder von oben nach unten
        self.nz = len(self.z_lev)

        # Wind speeds
        self.vH_ms = float(mat['Uh_AverageWindSpeed'][0])

        # Time scales
        self.fq_sp = float(mat['Sample_frequency'][0])
        self.dT = 1 / self.fq_sp
        self.T = float(mat['Sample_period'][0])
        self.nT = int(self.T * self.fq_sp)
        
        #Wind Direction of current data
        self.WindDirection = float(mat['Wind_direction_angle'][0])

        #Pressure coeff. 
        self.Cp = mat['Wind_pressure_coefficients']

        #--------------------------------------------------------------
        # Transform Pressure Coeff. in base force
        # -------------------------------------------------------------
        # Calculate forces on area [in kN]
        self.F_p = self.A_i * self.Cp * 0.5 * self.rho_air  * self.vH_ms  ** 2 * 10 ** -3      #F = (Cp*0.5*Roh*V**2)*A

        # Calculating base forces
        self.calcModelBaseForces_P()

        # Calculating floor forces
        self.calcModelFloorForces_P()

    
    def scaleFactors(self, buildingProp, loads):
    
        # Skaliert die Modelldaten auf die reale Modellgröße
        self.lambda_l = buildingProp.h_total/self.H_ms
        self.lambda_v = loads.vH_fs/self.vH_ms
        self.lambda_fq = self.lambda_v/self.lambda_l
        self.lambda_t = 1/self.lambda_fq

        # Skaliert Fußpunkt Kräfte und Moment von Modell in reale Größe
        self.lambda_F = self.lambda_v ** 2 * self.lambda_l ** 2
        self.lambda_M = self.lambda_v ** 2 * self.lambda_l ** 3


    def scaleData(self, buildingProp, loads):
        self.scaleFactors(buildingProp, loads)

        self.BFD_fs = []
        self.BFL_fs = []
        for i in range(0,len(self.BFD_ms)):
            self.BFD_fs.append(self.BFD_ms[i] * self.lambda_F)
            self.BFL_fs.append(self.BFL_ms[i] * self.lambda_F)

        self.BMD_fs = []
        self.BML_fs = []
        for i in range(0,len(self.BMD_ms)):
            self.BMD_fs.append(self.BMD_ms[i] * self.lambda_M)
            self.BML_fs.append(self.BML_ms[i] * self.lambda_M)
                
        self.T_fs = self.T*self.lambda_t
        self.fq_fs = self.fq_sp/self.lambda_fq
        self.dT = 1/self.fq_fs
        self.nT_fs = int(self.T_fs * self.fq_fs)

    def scalebuild(self):
        """scale wind tunnel model data to full scale
        """
        # Coordinates [in m], number and face of measurement points
        self.x_fs = self.x * self.lambda_l
        self.z_fs = self.z * self.lambda_l

        # Measurement area [in m2]
        self.A_i_fs = self.A_i * self.lambda_l ** 2

        # Get levels & sort them
        self.z_lev_fs = self.z_lev * self.lambda_l
        self.nz_fs = len(self.z_lev)

        # Scale floor forces
        self.LFD_fs = self.LFD_ms* self.lambda_F
        self.LFL_fs = self.LFL_ms* self.lambda_F

    def getClosest(self, node, z_lev):
        """Returns Position of node in z_lev_fs and the closest number bigger and smaller
         """
        
        pos = bisect.bisect_right(z_lev, node)
        
        if pos == 0:
            pos = len(self.z_lev_fs_r_inv) - pos - 2
            return pos, z_lev[1], z_lev[0]
        elif pos == len(z_lev):
            pos = len(self.z_lev_fs_r_inv) - pos
            return pos, z_lev[-1], z_lev[-2]

        # Richtig rum?
        higher = z_lev[pos]
        lower = z_lev[pos-1]

        pos = len(self.z_lev_fs_r_inv) - pos - 1

        return pos, higher, lower



    def interpolateForces(self, buildingProp):
        """interpolates full scale Floor forces to Floor geometry of building
        """

        #Coordinates [in m] for building floors (from top to bottom, cooordinate from bottom)
        self.nodes_coords_geom = np.linspace(buildingProp.h_total, 0, buildingProp.n+1)

        #Coordinates [in m] full scale for floor forces (from top to bottom, cooordinate from bottom)
        self.z_lev_fs_r = []
        for i in range(0,len(self.z_lev_fs)):
            self.z_lev_fs_r.append(round(self.z_lev_fs[i], 2))

        idx = len(self.z_lev_fs_r)-1
        self.z_lev_fs_r_inv = []

        while (idx >=0):
            self.z_lev_fs_r_inv.append(self.z_lev_fs_r[idx])
            idx -= 1


        self.LFD = [0]
        self.LFL = [0]

        for i in range(1,len(self.nodes_coords_geom)-1):
            pos, higher, lower = self.getClosest(self.nodes_coords_geom[i], self.z_lev_fs_r_inv)
            # Interpolieren
            f_D = self.LFD_fs[pos]+(self.LFD_fs[pos]-self.LFD_fs[pos+1])/(higher-lower)*(self.nodes_coords_geom[i]-higher)
            f_L = self.LFL_fs[pos]+(self.LFL_fs[pos]-self.LFL_fs[pos+1])/(higher-lower)*(self.nodes_coords_geom[i]-higher)
            self.LFD.append(f_D)
            self.LFL.append(f_L)

        anzahl_messpunkte = len(self.LFD_fs[0])
        f_fusspunkt = np.array([0]*anzahl_messpunkte)
        self.LFD.append(f_fusspunkt)                                              # 0 für Kraft an Fußpunkt anhängen
        self.LFL.append(f_fusspunkt)
        self.LFD[0] = self.LFD[1]                                       # an erstem Punkt über dem ersten Messpunkt wird der Wert aus dem Geschoss darunter übernommen
        self.LFL[0] = self.LFL[1]

        






