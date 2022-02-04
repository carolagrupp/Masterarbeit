# ------------------------------------------------------------------------------
# Description:  Berechnung Kerntragwerk
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-22      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke

# Co-Author:    Carola Grupp
# Created:      2021-11-02  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
from pfh import building
from pfh import calculations
from pfh import str_core
from pfh import fea
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):
    #in calculations.calcElementWidth
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    x=buildingProp.x

    Gd=loads.gd
    Qd=loads.qd
    gd_fassade=loads.gd_fassade
    Psi_w=loads.Psi_w
    Psi_q=loads.Psi_q
    M=loads.M       #M jeweils auf Geschosshöhe
    gamma_g=loads.gamma_g
    gamma_w=loads.gamma_w

    gamma=materialProp.gamma

    A_einzug=element.A_einzug
    l_fassade=element.l_fassade

    #M_max aus Wind für aktuelles Geschoss
    if element.typ=='Kern':
        M_max=gamma_w*M[s]          #s=Geschosszahl von oben

    else:
        M_max=0                     #Tragwerkstypen außer Kern tragen keinen Wind

    M_kombi=Psi_w*M_max

    if  s==n and x*n_abschnitt < n: #wenn n/n_abschnitt nicht aufgeht, unten noch mehrere Geschosse mit ni<n_abschnitt
        Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g

    else:
        Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g

    N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s
    N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s

    return N_max,N_kombi,M_max,M_kombi


def buildingStiffness(buildingProp,materialProp,kern,element2,element3,element4):
    #wird in calculations.buildingDeflection und interstoryDrift aufgerufen
    # I und ggf. A berechnen für EI und GA
    I=[]
    GA=[]

    for i in range (0,len(kern.t)):     #Länge = Anzahl Abschnitte, kern.t = benötigte Querschnitsdicke nach calcElementWidth
        calculations.calcProfileProp(kern, buildingProp, materialProp, kern.t[i])
            
        I.append(kern.I)    #nur ein Wert
        GA.append(kern.A_min/2*5/6*materialProp.G)      # Halbe Kernfläche, da nur die Wände in Kraftrichtung ("Stegwände") die Schubkraft aufnehmen. Abmidnerung mit Kappa (s. Engelke S.14)

    buildingProp.I=I        #Länge Anzahl Abschnitte
    buildingProp.GA=GA
   

def design(buildingProp,loads,materialProp,DataProp):
    #in app.py.MainWondow.MainCalculation
    # Elemente:
    b_raster=buildingProp.b_raster
    A=b_raster**2
    innenStütze=building.elements(A,0,'Stütze','Vollprofil')       # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStütze=building.elements(1/2*A,b_raster,'Stütze','Vollprofil')
    eckStütze=building.elements(1/4*A,b_raster,'Stütze','Vollprofil')
    kernOhnePFH=building.elements(8*A,0,'Kern ohne PFH','Kern')
    kern=building.elements(8*A,0,'Kern','Kern')

    str_=str_core

    # Tragfähigkeitsnachweise:
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kernOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kern,buildingProp,loads,materialProp,str_)
    
    t0=kern.t[-1]       #Kerndicke ganz unten (letztes Element der Liste)

    # Dynamischer Nachweis
    if loads.dynamicAnalysis == True:
        calculations.calcDynamicElementWidth(buildingProp, loads, materialProp, DataProp, str_, kern)

    # Gebäudenachweise:
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,kern)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,kern)
    buildingProp.t_kern=kern.t

    #Vergleich Tragfähigkeit und Gebäudenachweise:
    if kern.t[-1] > t0:     
        buildingProp.NW_maßgebend='Verformung'
        print(kern.t[-1],t0)
    
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'

    # Massen:
    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,randStütze)
    calculations.calcWeight(buildingProp,materialProp,eckStütze)
    calculations.calcWeight(buildingProp,materialProp,kernOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,kern)

    buildingProp.G_innenStützen=innenStütze.G
    buildingProp.G_randStützen=[element*12 for element in randStütze.G]
    buildingProp.G_eckStützen=[element*4 for element in eckStütze.G]
    buildingProp.G_aussteifung=kern.G
    buildingProp.G_kernOhnePFH=kernOhnePFH.G

    buildingProp.G_außenStützen=[]
    buildingProp.G_total=[]
    buildingProp.G_totalOhnePFH=[]

    for i in range (0,len(innenStütze.G)):
        G_außenStützen=randStütze.G[i]*12+eckStütze.G[i]*4
        buildingProp.G_außenStützen.append(G_außenStützen)
        buildingProp.G_total.append(G_außenStützen+innenStütze.G[i]+kern.G[i])
        buildingProp.G_totalOhnePFH.append(G_außenStützen+innenStütze.G[i]+kernOhnePFH.G[i])
    
    
    # Eigenfrequenz
    calculations.calcEigenfrequency(buildingProp,loads,materialProp)
    
    buildingProp.t_innenStützen=innenStütze.t
    buildingProp.t_randStützen=randStütze.t
    buildingProp.t_eckStützen=eckStütze.t
    buildingProp.t_kern=kern.t
    buildingProp.t_diagonale=[]
    buildingProp.t_querstrebe=[]
    buildingProp.t_riegel=[]
    buildingProp.t_outrigger = []
    buildingProp.t_belt = []

    #feModel = fea.feModel(buildingProp, loads, materialProp)
    #moment = fea.feModel.calcBendingMoment(feModel)
 