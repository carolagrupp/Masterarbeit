# ------------------------------------------------------------------------------
# Description:  Objekte in denen alle Eigenschaften des Gebäudes, der Elemente, 
#               des Materials, die Lasten und Grenzwerte gespeichert werden
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-09      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke

# Co-Author:    Carola Grupp
# Created:      2021-11-02  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:
import numpy as np
# ------------------------------------------------------------------------------


class buildingProp():
    def __init__(self):
        # Geometric properties
        #Eingaben
        
        self.n=0 #[-]
        self.n_abschnitt=0
        self.h_total=0
        self.h_geschoss=1
        self.posOut=[]
        self.z=0        #Anzahl Iterationen Tragfähigkeitsnachweis

        
        self.I=[]
        self.GA=[]

        self.tragwerk='-'
        self.G_außenStützen=[0]
        self.G_innenStützen=[0]
        self.G_aussteifung=[0]
        self.G_total=[0]
        self.G_decken=0
        self.K = []             #Liste Federsteifigkeiten Outrigger

        self.D = 0

        self.multi_n=[]
        self.multi_h_total=[]
        self.multi_schlankheit=[]
        
        self.multi_G_aussteifung=[]
        self.multi_G_außenStützen=[]
        self.multi_G_innenStützen=[]
        self.multi_G_total=[]
        self.multi_G_totalOhnePFH=[]
        self.multi_G_decken=[]

        self.multiPar_G_aussteifung=[]
        self.multiPar_G_außenStützen=[]
        self.multiPar_G_innenStützen=[]
        self.multiPar_G_total=[]
        self.multiPar_G_decken=[]
        self.multiPar_G_totalOhnePFH=[]
        self.multi_p=[]
        self.parameter=''

        self.multi_G_aussteifung_opt=[]
        self.multi_G_außenStützen_opt=[]
        self.multi_G_innenStützen_opt=[]
        self.multi_G_total_opt=[]
        self.multi_G_decken_opt=[]
        self.multi_G_totalOhnePFH_opt=[]
        self.multi_p_opt=[]

        self.w_EI=[0]
        self.w_GA=[0]

        self.Teta_i_EI=[0]
        self.Teta_i_GA=[0]

        self.multiPar_w_EI = []
        self.multiPar_w_GA = []
        self.multiPar_w = []

        self.multiPar_t_soll_1 = []
        self.multiPar_t_rand_1 = []
        self.multiPar_t_eck_1 = []     
        self.multiPar_t_soll_2 = []    
        self.multiPar_t_rand_2 = []
        self.multiPar_t_eck_2 = []    
        self.multiPar_t_soll_3 = []    
        self.multiPar_t_rand_3 = []
        self.multiPar_t_eck_3 = []    
        self.multiPar_t_soll_4 = []    
        self.multiPar_t_rand_4 = []
        self.multiPar_t_eck_4 = []
        
        self.multi_eigenFrequenz=[]

        self.i_wechselNW='none'

        self.alpha_outriggerFulfilled = True
        # 2D Parameteroptimierung
        self.varPosOut = False
        self.varPosOut_2D = False

        self.multiPar_2D_G_aussteifung = []
        self.multiPar_2D_G_außenStützen = []
        self.multiPar_2D_G_innenStützen = []
        self.multiPar_2D_G_total = []
        self.multiPar_2D_G_decken = []
        self.multiPar_2D_G_totalOhnePFH = []
        self.multiPar_2D_w_EI = []
        self.multiPar_2D_w_GA = []
        self.multiPar_2D_w = []

        # Grundriss:
        # ES   RS   RS   RS   ES
        # RS   IS   IS   IS   RS
        # RS   IS   IS   IS   RS
        # RS   IS   IS   IS   RS
        # ES   RS   RS   RS   ES

    
    def calcFloorWeight(self,materialProp):

        self.G_decken=(self.b_raster*4)**2*self.n*materialProp.G_decke/1000



class loads():
    def __init__(self):
        self.M=[0]
        self.V=[0]
        self.windData = False
        self.dynamicAnalysis = False

        '''__init__(self,gd,qd,gd_fassade,GK,qb,Psi_q,Psi_w):
        self.gd=gd #[kN/m²]
        self.qd=qd #[kN/m²]
        self.gd_fassade=gd_fassade #[kN/m]

        self.GK=GK
        self.qb=qb #[kN/m²]

        self.Psi_q=Psi_q #[-]
        self.Psi_w=Psi_w #[-]'''

    def calculateInternalForces(self,buildingProp):
        b_total=buildingProp.b_total
        n=buildingProp.n
        h_geschoss=buildingProp.h_geschoss
        h=buildingProp.h_total
        
        GK=self.GK
        wk=self.wk #[kN/m²]
        
        self.M=[] #[kNm]
        self.V=[] #[kN]     
        self.F_p=[] #[kN]

        # Abweichende Verläufe/Parameter je nach Geländekategorie
        #EC1-1-4 NA Tab. NA B.2 z_min, a, b, c, d: aus Integration von qp(z)*btotal, so ergibt sich V aus einer Integration und M aus der erneuten Integration von V
        if GK==1:
            a=0.644143/2.6          # Der Faktor 2.6 ist bereits in wk eingerechnet s. app.py Z.408 w_k
            b=1.410673/2.6
            c=0.766530/2.6
            exp=1.19
            z_min=2
            d=1.9/2.6

        if GK==2:
            a=0.435060/2.1
            b=0.974535/2.1
            c=0.539475/2.1
            exp=1.24
            z_min=4
            d=1.7/2.1

        if GK==3:
            a=0.258962/1.6
            b=0.598203/1.6
            c=0.339241/1.6
            exp=1.31
            z_min=8
            d=1.5/1.6

        if GK==4:
            a=0.130333/1.1
            b=0.312798/1.1
            c=0.182466/1.1
            exp=1.40
            z_min=16
            d=1.3/1.1

        if GK==0:      # Gleichlast, als Berechnungstest
            a=0.5
            b=1
            c=0.5
            exp=1
            z_min=0
            d=1

        # Berücksichtigung von Belastung z<zmin? Eher rausstreichen, da nicht relevant
        M_zmin=wk*b_total*a*z_min**(exp+1)-wk*b_total*b*h**exp*z_min+wk*b_total*c*h**(exp+1)
        V_zmin=wk*b_total*b*z_min**exp-wk*b_total*b*h**exp

        # Momente auf Höhe der Geschosse
        for i in range(0,n+1):
            z=h-i*h_geschoss #Höhe von unten
            
            if z < z_min:
                M_i=1/2*wk*b_total*d*(z-z_min)**2+V_zmin*(z-z_min)+M_zmin
            else:
                M_i=wk*b_total*a*z**(exp+1)-wk*b_total*b*h**exp*z+wk*b_total*c*h**(exp+1)
            

            self.M.append(M_i)
            #self.M.append(0)

        # Querkräfte auf Höhe der Geschosse
        for i in range(0, n+1):
            z=h-i*h_geschoss

            if z < z_min:
                V_i=wk*b_total*d*(z-z_min)+V_zmin
            else:
                V_i=wk*b_total*b*z**exp-wk*b_total*b*h**exp

            self.V.append(-V_i) # Vorzeichen wird umgekehrt damit nachher positive Ergebnisse

        # Knotenkräfte
        self.V.append(0) # Damit Knotenkraft am obersten Knoten berechnet werden kann

        for i in range (0,n):
            F=1/3*(self.V[i]-self.V[i-1])+2/3*(self.V[i+1]-self.V[i]) 
            # Die Differenz der Querkräfte entspricht der Windlast die in dem Abschnitt angreift
            # Die Windlast jedes Abschnitts wird auf der sicheren Seite liegen zu 2/3 auf den oberen Knoten und zu 1/3 auf den unteren Knoten aufgeteilt
            self.F_p.append(F)
        
        self.F_p.append(0)  # am untersten Knoten soll die Kraft Null sein (geht direkt ins Lager-keine Verformung)

        del self.V[-1] # löscht die 0 wieder aus V

    def calcvH_fs(self, buildingProp):
        #mittlere Windgeschwindigkeit vb auf 10m Höhe mit Wiederkehrwahrscheinlichkeit von 50 Jahren (DIN 3.4)

        h = buildingProp.h_total
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
        self.vH_fs = cprob * v_fs

    def grabData(self, buildingProp, DataProp, direction, schlankheit):

        # Schnittkräfte aus Normberechnung löschen
        self.M=[0]
        self.V=[0]

        if schlankheit == 8:
            DataProp.Data = 'Force'
            DataProp.scale_ms = 1000
            if self.alpha_v == 0.16:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\8\\015\\time_series_1_{}.mat".format(direction)
            elif self.alpha_v == 0.30:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\8\\027\\time_series_1_{}.mat".format(direction)
        elif schlankheit == 4:
            DataProp.Data = 'Pressure'
            DataProp.scale_ms = 400
            if self.alpha_v == 0.16:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\4\\1_6\\time_series_of_point_wind_pressure_{}.mat".format(direction)
            elif self.alpha_v == 0.30:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\4\\1_4\\time_series_of_point_wind_pressure_{}.mat".format(direction)
        elif schlankheit == 5:
            DataProp.Data = 'Pressure'
            DataProp.scale_ms = 400
            if self.alpha_v == 0.16:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\5\\1_6\\time_series_of_point_wind_pressure_{}.mat".format(direction)
            elif self.alpha_v == 0.30:
                fname = "C:\\Users\\carol\\github\\Masterarbeit\\WindData\\5\\1_4\\time_series_of_point_wind_pressure_{}.mat".format(direction)
        else:
            raise('Wähle eine Gebäudeschlankheit von 4,5 oder 8')

        # Calculate wind velocity full scale
        self.calcvH_fs(buildingProp)

        # Load wind tunnel model properties, TPU Database files
        if DataProp.Data == 'Force':
            DataProp.readData_F(fname, buildingProp)
        if DataProp.Data == 'Pressure':
            DataProp.readData_P(fname)

        # Scale Data
        DataProp.scaleData(buildingProp, self)

        # Scale wind tunnel model
        if DataProp.Data == 'Pressure':
            DataProp.scalebuild()
            DataProp.interpolateForces(buildingProp)





class materialProp():
    def __init__(self):
        
        self.E=0 #[kN/m²]
        self.GWP_decke=0
        self.PET_decke=0
        self.GWP_tragwerk=0
        self.PET_tragwerk=0
        #self.G=0 #[kN/m²]
        #self.gamma=0 #[kN/m³]

class elements():
    def __init__(self,A_einzug,l_fassade,typ,profil, typGenau):
        # Werte werden nur zur Übersicht hier angegeben
        self.A_einzug=A_einzug #[m²]        #Eingabe in str_.design
        self.l_fassade=l_fassade #[m]
        self.typ=typ
        self.typGenau = typGenau
        self.profil=profil
