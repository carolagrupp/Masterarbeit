# ------------------------------------------------------------------------------
# Description:  Berechnung Kerntragwerk
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-22      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
from pfh import building
from pfh import calculations
from pfh import str_core
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    x=buildingProp.x

    Gd=loads.gd
    Qd=loads.qd
    gd_fassade=loads.gd_fassade
    Psi_w=loads.Psi_w
    Psi_q=loads.Psi_q
    M=loads.M
    gamma_g=loads.gamma_g
    gamma_w=loads.gamma_w

    gamma=materialProp.gamma

    A_einzug=element.A_einzug
    l_fassade=element.l_fassade

    if element.typ=='Kern':
        M_max=gamma_w*M[s]

    else:
        M_max=0    

    M_kombi=Psi_w*M_max

    if  s==n and x*n_abschnitt < n:
        Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g

    else:
        Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g

    N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s
    N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s

    return N_max,N_kombi,M_max,M_kombi


def buildingStiffness(buildingProp,materialProp,kern,element2,element3,element4):
    # I und ggf. A berechnen für EI und GA
    I=[]
    GA=[]

    for i in range (0,len(kern.t)):
        calculations.calcProfileProp(kern, buildingProp, materialProp, kern.t[i])
            
        I.append(kern.I)
        GA.append(kern.A_min/2*5/6*materialProp.G)      # Halbe Kernfläche, da nur die Wände in Kraftrichtung ("Stegwände") die Schubkraft aufnehmen. Abmidnerung mit Kappa (oder näherung nach Land? nur stegfläche)

    buildingProp.I=I
    buildingProp.GA=GA
   

def design(buildingProp,loads,materialProp):
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
    
    t0=kern.t[-1]

    # Gebäudenachweise:
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,kern)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,kern)
    buildingProp.t_kern=kern.t

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
 