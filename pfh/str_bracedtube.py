# ------------------------------------------------------------------------------
# Description:  Berechnung Rahmentragwerk
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-05-13      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
#from premiumForHeight.pfh.building import elements
from pfh import building
from pfh import calculations
from pfh import str_bracedtube
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    x=buildingProp.x
    b_total=buildingProp.b_total

    Gd=loads.gd
    Qd=loads.qd
    gd_fassade=loads.gd_fassade
    M=loads.M
    V=loads.V
    gamma_g=loads.gamma_g
    gamma_w=loads.gamma_w

    gamma=materialProp.gamma

    A_einzug=element.A_einzug
    l_fassade=element.l_fassade

    if element.typ=='Diagonale':
        M_max=0
        M_kombi=0
        
        a=min(int(s/8+7/8)*8-4,n)    # Dadurch wird immer bis zu nächsten 8er Stufe aufgerundet, aber maximal bis zur Gesamtgeschossanzahl
        h=h_geschoss*8
        b=b_total
        d=(h**2+b**2)**0.5

        N_max=gamma_w*V[a]*d/b/4        # 4 Diagonalen tragen Last ab
        N_kombi=0


    if element.typ=='Mittelstütze':
        M_max=0
        M_kombi=0
        
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g

        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s
        N_kombi= 0


    if element.typ=='Außenstütze':
        M_max=0
        M_kombi=0
        
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            i=-1

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)

        a=min(int(s/4+3/4)*4,n)    # Dadurch wird immer bis zu nächsten 4er Stufe aufgerundet, aber maximal bis zur Gesamtgeschossanzahl
        F=gamma_w*M[a]/(b_total*5)          # Kraft aus Kräftepaar aus Moment

        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+F+2/5*buildingProp.G_diagonale[i]*gamma_g+3/8*buildingProp.G_querstrebe[i]*gamma_g
        N_kombi= 0


    if element.typ=='Außenstütze ohne PFH':
        M_max=0
        M_kombi=0
        
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            i=-1

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)

        a=min(int(s/4+3/4)*4,n)    # Dadurch wird immer bis zu nächsten 4er Stufe aufgerundet, aber maximal bis zur Gesamtgeschossanzahl
        F=0                        # Kraft aus Kräftepaar aus Moment

        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+F
        N_kombi= 0

    if element.typ=='Querstrebe':
        M_max=0
        M_kombi=0
        a=min(int(s/4+3/4)*4,n)

        N_max=gamma_w*(M[a]-M[a-4])/(b_total*5)*b_total/(h_geschoss*8)
        N_kombi=0

    return N_max,N_kombi,M_max,M_kombi


def buildingStiffness(buildingProp,materialProp,diagonale,außenStützen,element3,element4):
    # I und ggf. A berechnen für EI und GA
    I=[]
    GA=[]

    for i in range (0,len(außenStützen.t)):
        calculations.calcProfileProp(außenStützen,buildingProp,materialProp,außenStützen.t[i])
        calculations.calcProfileProp(diagonale,buildingProp,materialProp,diagonale.t[i])
            
        I_i=10*außenStützen.A*(buildingProp.b_raster*2)**2+4*außenStützen.A*(buildingProp.b_raster)**2#+16*außenStützen.I     # Steiner Anteil der Randstützen
        I.append(I_i)

        h=buildingProp.h_geschoss*8
        b=buildingProp.b_total
        d=(h**2+b**2)**0.5
        GA_i=4*materialProp.E*h*b**2*diagonale.A/d**3       # Schubsteifigkeit 2er Fachwerke nach Kaviani et al.
        GA.append(GA_i)
    
    buildingProp.I=I
    buildingProp.GA=GA


def shearStiffnessModification(buildingProp,diagonale,außenStützen,element3,element4,delta_t):
    
    diagonale.t = [element+delta_t for element in diagonale.t]


def bendingStiffnessModification(buildingProp,diagonale,außenStützen,element3,element4,delta_t):
    
    außenStützen.t = [element+delta_t for element in außenStützen.t]


def design(buildingProp,loads,materialProp):
    # Elemente:
    b_raster=buildingProp.b_raster
    A=b_raster**2
    innenStütze=building.elements(A,0,'Mittelstütze','Vollprofil')       # Einzugsfläche, Fassadenlänge, Typ, Profil
    diagonale=building.elements(0,0,'Diagonale','Quadratisches Hohlprofil')
    querstrebe=building.elements(0,0,'Querstrebe','Quadratisches Hohlprofil')
    außenStütze=building.elements(A/2,b_raster,'Außenstütze','Quadratisches Hohlprofil')
    außenStützeOhnePFH=building.elements(A/2,b_raster,'Außenstütze ohne PFH','Quadratisches Hohlprofil')

    str_= str_bracedtube

    # Tragfähigkeitsnachweise:
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(außenStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(diagonale, buildingProp, loads, materialProp, str_)
    calculations.calcElementWidth(querstrebe,buildingProp,loads,materialProp,str_)

    calculations.calcWeight(buildingProp,materialProp,diagonale)
    h=buildingProp.h_geschoss*8
    b=buildingProp.b_total
    d=(h**2+b**2)**0.5
    diagonale.G=[element*d/b for element in diagonale.G]

    buildingProp.G_diagonale=diagonale.G

    calculations.calcWeight(buildingProp,materialProp,querstrebe)
    x=2*b_raster/(4*buildingProp.h_geschoss)                                    # Alle 4 Geschosse gibt es eine Querstrebe mit l=2*b_raster
    querstrebe.G=[element*x for element in querstrebe.G]
    
    buildingProp.G_querstrebe=querstrebe.G

    calculations.calcElementWidth(außenStütze,buildingProp,loads,materialProp,str_)
    
    ts0=außenStütze.t[-1]
    td0=diagonale.t[-1]

    # Gebäudenachweise:
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,diagonale,außenStütze)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,diagonale,außenStütze)

    if außenStütze.t[-1] > ts0 or diagonale.t[-1] > td0:
        buildingProp.NW_maßgebend='Verformung'
    
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'
    
    # Massen:
    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,außenStützeOhnePFH)
    
    calculations.calcWeight(buildingProp,materialProp,diagonale)
    diagonale.G=[element*d/b for element in diagonale.G]

    calculations.calcWeight(buildingProp,materialProp,außenStütze)
    
    buildingProp.G_innenStützen=[element*9 for element in innenStütze.G]
    buildingProp.G_außenStützen=[element*16 for element in außenStütze.G]
    buildingProp.G_diagonalen=[element *8 for element in diagonale.G]
    buildingProp.G_querstreben=[element*4 for element in querstrebe.G]
    
    buildingProp.G_aussteifung=[]
    buildingProp.G_total=[]
    buildingProp.G_totalOhnePFH=[]

    for i in range (0,len(innenStütze.G)):
        buildingProp.G_aussteifung.append(buildingProp.G_diagonalen[i]+buildingProp.G_querstreben[i])
        buildingProp.G_total.append(buildingProp.G_innenStützen[i]+buildingProp.G_diagonalen[i]+buildingProp.G_außenStützen[i]+buildingProp.G_querstreben[i])
        buildingProp.G_totalOhnePFH.append(buildingProp.G_innenStützen[i]+außenStützeOhnePFH.G[i]*16)

    # Eigenfrequenz

    calculations.calcEigenfrequency(buildingProp,loads,materialProp)

    buildingProp.t_innenStützen=innenStütze.t
    buildingProp.t_randStützen=außenStütze.t
    buildingProp.t_eckStützen=außenStütze.t
    buildingProp.t_kern=[]
    buildingProp.t_diagonale=diagonale.t
    buildingProp.t_querstrebe=querstrebe.t
    buildingProp.t_riegel=[]
