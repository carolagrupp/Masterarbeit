# ------------------------------------------------------------------------------
# Description:  Import Winddaten
#
# ------------------------------------------------------------------------------
# Author:    Carola Grupp
# Created:      2022-01-07  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports: 
from re import A
import numpy as np
#from Masterarbeit.pfh.building import buildingProp
#from pfh import building
from pfh import fea
from WindData import data
from WindData import response
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

def grabData(str_, element1,element2=None,element3=None,element4=None):     #buildingProp, materialProp, loads
    Directions = [0,5,10,15,20,25,30,35,40,45]
    buildingProp = building()   #für Implementierung
    DataProp = data.DataProp()

    schlankheit = 4     #Für Anfang Implementierung
    
    #schlankheit = buildingProp.schlankheit
    
    #for i in Directions:
        #Durchführung aller Schritte
    directions = str(Directions[0])
    if schlankheit == 8:
        DataProp.Data = 'Force'
        DataProp.scale_ms = 1000
        if buildingProp.alpha_v == 0.16:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\8\\015\\time_series_1_{}.mat".format(directions)
        elif buildingProp.alpha_v == 0.30:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\8\\027\\time_series_1_{}.mat".format(directions)
    elif schlankheit == 4:
        DataProp.Data = 'Pressure'
        DataProp.scale_ms = 400
        if buildingProp.alpha_v == 0.16:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\4\\1_6\\time_series_of_point_wind_pressure_{}.mat".format(directions)
        elif buildingProp.alpha_v == 0.30:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\4\\1_4\\time_series_of_point_wind_pressure_{}.mat".format(directions)
    elif schlankheit == 5:
        DataProp.Data = 'Pressure'
        DataProp.scale_ms = 400
        if buildingProp.alpha_v == 0.16:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\5\\1_6\\time_series_of_point_wind_pressure_{}.mat".format(directions)
        elif buildingProp.alpha_v == 0.30:
            fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\5\\1_4\\time_series_of_point_wind_pressure_{}.mat".format(directions)
    else:
        raise('Wähle eine Gebäudeschlankheit von 4,5 oder 8')
        #+Abbruch
    
    # Calculate wind velocity full scale
    buildingProp.calcvH_fs()

    # Load wind tunnel model properties, TPU Database files
    if DataProp.Data == 'Force':
        DataProp.readData_F(fname, buildingProp)
    if DataProp.Data == 'Pressure':
        DataProp.readData_P(fname)

    # Scale Data
    DataProp.scaleData(buildingProp)

    # Scale wind tunnel model
    if DataProp.Data == 'Pressure':
        DataProp.scalebuild()


    w_max = buildingProp.w_max  # nachher loads.w_max
    a_max = buildingProp.a_max

    # hier Massenberechnung mit aktuellen ts des Tragwerks, vorher schon Bemessung auf Vertikallasten über Calcelementloads?
    #buildingProp.calcBuildMass()    # Brauche ich das?
    # nur Stiffnessberechnung des aktuellen Systems für Eigenfrequenz
    str_.buildingStiffness(buildingProp, materialProp, element1, element2, element3, element4)

    # Create analysis model
    feModel = fea.feModel(buildingProp)

    # Compute eigenfrequencies
    feModel.getEigenfrequency()

    # Calc generalized quantities (K_gen, M_gen)
    feModel.calcGeneralizedQuantitites()

    # Calc response forces
    responseForces_D = response.responseForces(DataProp.BMD_fs, DataProp.dT, feModel.fq_e, buildingProp.D, feModel.fq_e, 360)
    #Wieso zweimal fq_e als Input? response nimmt zweites mal als nue, weiso 360?
    responseForces_L = response.responseForces(DataProp.BML_fs, DataProp.dT, feModel.fq_e, buildingProp.D, feModel.fq_e, 360)


    # Iterativ: Anfang bis Response + folgendes
    # Calc response deflections
    # LFL/LFD = FLoor Forces: Wie Berechnung für HFFB Modelldaten?
    responseDeflections_D = response.responseDeflection(feModel, responseForces_D, DataProp.LFD_fs)
    w_tip_D = responseDeflections_D.delta_tip_r_max

    if w_tip_D > w_max:
        t = t + a

    responseDeflections_L = response.responseDeflection(feModel, responseForces_L, DataProp.LFL_fs)
    w_tip_L = responseDeflections_L.delta_tip_r_max
    
    if w_tip_L > w_max:
        t = t + a

    #Erhöhen von t, wenn w_tip größer als w_max



    # Iterativ: Anfang bis Response + folgendes
    # Calc response accelerations
    responseAccelerations_D = response.responseAccelerations(feModel, DataProp.BMD_fs, DataProp.dT, feModel.fq_e, buildingProp.D, feModel.fq_e, 360)
    a_rms_D = responseAccelerations_D.a_r_rms

    responseAccelerations_L = response.responseAccelerations(feModel, DataProp.BML_fs, DataProp.dT, feModel.fq_e, buildingProp.D, feModel.fq_e, 360)
    a_rms_L = responseAccelerations_L.a_r_rms

    # Erhöhen von t, wenn a größer als a_max




class building():
    def __init__(self):
        #Schlankheit 1:8
        #self.h_total = 256 # Schlankheit von 8
        self.h_total = 128 # Schlankheit von 4
        self.b_total = 32
        self.n_abschnitt = 8
        self.h_geschoss = 4
        self.n = 64
        self.D = 0 #Damping
        self.v_bk = 22.5                            #m/s        nachher aus gui
        self.GK = 2 #=loads.GK
        self.w_max = self.h_total/600 # = loads.w_max
        self.a_max = 10 * 0.003 # ISO 10137 1-yr return limit for offices
        if self.GK == 2:
            self.alpha_v = 0.16
        elif self.GK == 4:
            self.alpha_v = 0.30

    def calcBuildMass(self):
        #building masses (in t) Dummywerte
        self.M_DL_Floor = 25600         # ständige Lasten aus Decken
        self.M_DL_Col = 401             # ständige Lasten aus Stützen
        self.M_DL_IWall = 3200          # ständige Lasten aus Kern?

        self.M_DL_CWall = 0

        self.M_SDL_Floor = 1.0 * 0.100 * (self.b_total**2) * self.n        # veränderliche Lasten aus Decken
        self.M_LL_Floor = 0.2 * 0.250 * (self.b_total**2) * self.n         # Nutzlasten 

        self.M_DL_Tot = self.M_DL_Floor + self.M_DL_Col + self.M_DL_IWall + \
            self.M_DL_CWall
        self.M_SLS_Tot = self.M_DL_Floor + self.M_DL_Col + self.M_DL_IWall + \
            self.M_DL_CWall + self.M_SDL_Floor + self.M_LL_Floor

    def recalcBuildMass(self):
        self.M_DL_Tot = self.M_DL_Floor + self.M_DL_Col + self.M_DL_IWall + \
            self.M_DL_CWall
        self.M_SLS_Tot = self.M_DL_Floor + self.M_DL_Col + self.M_DL_IWall + \
            self.M_DL_CWall + self.M_SDL_Floor + self.M_LL_Floor

    def calcvH_fs(self):
        #mittlere Windgeschwindigkeit vb auf 10m Höhe mit Wiederkehrwahrscheinlichkeit von 50 Jahren (DIN 3.4)

        h = self.h_total
        if self.GK == 2:
            v_fs = 1.00*self.v_bk*(h/10)**0.16
        elif self.GK == 4:
            v_fs = 0.56*self.v_bk*(h/10)**0.30
        
        # Get return period / probability of excedence
        R = 50
        p = 1 / R

        # According to DIN EN 1991-1-4 / NA Abs. 4.2, Anmerkung 5
        K = 0.1
        n = 1

        # Get probability factor (DIN EN 1991-1-4 gl 4.2)
        cprob   = ((1 - K * np.log(-np.log(1-p))) / (1 - K * np.log(-np.log(0.98)))) ** n

        # Calculate wind speeds at different return periods
        self.vH_fs = cprob * v_fs #vielleicht in building damit buildingProp.vH_fs


    











