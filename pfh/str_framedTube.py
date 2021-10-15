# ------------------------------------------------------------------------------
# Description:  Berechnung Framed Tube
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-07-23      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
#from premiumForHeight.pfh.building import buildingProp, materialProp
from numpy.lib import stride_tricks
from pfh import building
from pfh import calculations
from pfh import str_framedTube
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    b_total=buildingProp.b_total
    n=buildingProp.n
    n_stiele=buildingProp.n_stiele
    x=buildingProp.x

    Gd=loads.gd
    Qd=loads.qd
    gd_fassade=loads.gd_fassade
    Psi_w=loads.Psi_w
    Psi_q=loads.Psi_q
    M=loads.M
    V=loads.V
    gamma_w=loads.gamma_w
    gamma_g=loads.gamma_g

    gamma=materialProp.gamma

    A_einzug=element.A_einzug
    l_fassade=element.l_fassade


    if element.typ=='Riegel':
        M_max=1/2*1/(4*n_stiele)*h_geschoss/2*gamma_w*V[s]     
        M_kombi=0
        N_max=0
        N_kombi=0


    if element.typ=='Stiel':
        M_max=gamma_w*V[s]/(16*n_stiele)*h_geschoss/2      
        M_kombi=Psi_w*gamma_w*V[s]/(16*n_stiele)*h_geschoss/2
        
        I,G_w,t_w=calcInteriaMoment(buildingProp,materialProp,s,element,t)

        sigma_M=gamma_w*M[s]/I*b_total/2         # Anteil des Moments der nicht über ein Moment in den Stützen sondern über ein Kräftepaar abgetragen wird
        F=sigma_M*element.A                          # Kraft in den Außenstützen aus Kräftepaar, wenn Außenstützen dreimal so viel Kraft aufnehmen wie Innenstützen

        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            calculations.calcWeight(buildingProp, materialProp, buildingProp.riegel)
            Ng_riegel=buildingProp.riegel.G[-1]*gamma_g*10

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)
            calculations.calcWeight(buildingProp, materialProp, buildingProp.riegel)
            Ng_riegel=buildingProp.riegel.G[i]*gamma_g*10
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+Ng_riegel+Psi_w*F
        N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s+Ng_riegel+F


    if element.typ=='Stütze':
        M_max=0      
        M_kombi=0
       
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
        
        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s
        N_kombi= 0
    

    return N_max,N_kombi,M_max,M_kombi


def calcInteriaMoment(buildingProp,materialProp,s,stiel,t_stiel):
    # Berechnung von G_w:
    if s < buildingProp.n:
        i=int(s/buildingProp.n_abschnitt-1)
    
    else:
        i=-1

    riegel=buildingProp.riegel

    calculations.calcProfileProp(riegel,buildingProp,materialProp,riegel.t[i])
    calculations.calcProfileProp(stiel,buildingProp,materialProp,t_stiel)
    
    h=buildingProp.h_geschoss
    d_b=0# riegel.t[i]/100 #(riegel.t[i]*materialProp.verhältnis_td*2)/100
    d_c=0#t_stiel/100 #(t_stiel*materialProp.verhältnis_td*2)/100
    E=materialProp.E
    G=materialProp.G
    sp=buildingProp.b_raster/buildingProp.n_stiele
    
    d1=(h-d_b)**3/(12*E*stiel.I)
    d2=(h/sp)**2*(sp-d_c)**3/(12*E*riegel.I)
    
    d3=(h-d_b)/(G*stiel.A/2*5/6)            # Schubfläche = Halbe Querschnittsfläche * Schubkorrektiurfaktor 5/6
    d4=(h/sp)**2*(sp-d_c)/(G*riegel.A*2/3*5/6)    # Stegfläche ist 2/3 der Gesamtfläche, nicht 1/2 da Hohlprofil

    G_w=(h/stiel.A)/(d1+d2+d3+d4)

    h_total=buildingProp.h_total
    a=buildingProp.b_total/2
    
    m_w=(G_w*h_total**2)/(E*a**2)

    alpha_1=(2.57*m_w+1.12)/(m_w**2+2.94*m_w+0.64)      # Gleichmäßig verteilte Last nach Kwon
    alpha_2=(0.03*m_w+1.12)/(m_w**2+2.94*m_w+0.64)

    beta_1=(7.72*m_w+14.15)/(m_w**2+12.35*m_w+11.32)
    beta_2=(0.08*m_w+14.15)/(m_w**2+12.35*m_w+11.32)

    z=h_total-s*h

    alpha=alpha_1*(1-z/h_total)**2+alpha_2*(2*z/h_total-(z/h_total)**2)
    beta=beta_1*(1-z/h_total)**2+beta_2*(2*z/h_total-(z/h_total)**2)

    k=1/3*(1-2/5*alpha)+(1-2/3*beta)    # Korrekturfaktor
    #k2=(1-2/3*beta)+1/3                 # Korrekturfaktor ohne alpha

    tube_wirksam=k/(4/3)*100        # Wirksamkeit der Tube in Prozent
    #print('alpha=',alpha,'beta=',beta)
    t= stiel.A/sp
    
    I=4*t*a**3*k

    buildingProp.tube_d1=d1
    buildingProp.tube_d2=d2
    buildingProp.tube_d3=d3
    buildingProp.tube_d4=d4
    buildingProp.tube_wirksam=tube_wirksam

    return I,G_w,t



def buildingStiffness(buildingProp,materialProp,riegel,stiel,element3,element4):
    # I und ggf. A berechnen für EI und GA
    I=[]
    GA=[]

    for i in range (0,len(riegel.t)):
        
        if i==(len(riegel.t)-1):
            s=buildingProp.n       # s für Stelle an der berechnet wird
        
        else: 
            s=buildingProp.n_abschnitt*(i+1)

        #s=40

        buildingProp.riegel=riegel
        I_i,G_w,t=calcInteriaMoment(buildingProp,materialProp,s,stiel,stiel.t[i])

        calculations.calcProfileProp(riegel,buildingProp,materialProp,riegel.t[i])
        calculations.calcProfileProp(stiel,buildingProp,materialProp,stiel.t[i])

        I_i=16*buildingProp.n_stiele*stiel.I+stiel.A*((buildingProp.n_stiele*8+2)*(buildingProp.b_raster*2)**2+(8*buildingProp.n_stiele-4)*(buildingProp.b_raster)**2)
        #I_i=4*t*(buildingProp.b_total/2)**3*4/3
        I.append(I_i) # I bei Steifigkeit ohne Shear-Lag

        # d_b=0#riegel.t[i]/100
        # d_c=0#stiel.t[i]/100

        # h=buildingProp.h_geschoss
        # sp=buildingProp.b_raster/buildingProp.n_stiele

        # GA_riegel=2*4*buildingProp.n_stiele*12*materialProp.E*riegel.I/((sp-d_c)**3*h/sp**2)
        # GA_stiel=(buildingProp.n_stiele*8+2)*12*materialProp.E*stiel.I/((h-d_b)**3/h)    # GA_stiel=Summe(12*EI/h^2)=12*E/h^2*Summe(I)
        # GA_riegel2=2*4*buildingProp.n_stiele*(materialProp.G*riegel.A*2/3*5/6)/(buildingProp.h_geschoss/buildingProp.b_raster/buildingProp.n_stiele)
        # GA_stiel2=(buildingProp.n_stiele*8+2)*stiel.A*(materialProp.G*0.5*5/6)*h/(h-d_b)

        # GA_i=(GA_riegel**-1+GA_stiel**-1+GA_riegel2**-1+GA_stiel2**-1)**-1

        GA_i=G_w*buildingProp.b_total*t*2       # G*Stegfläche
        GA.append(GA_i)
    
    buildingProp.I=I
    buildingProp.GA=GA


def shearStiffnessModification (buildingProp,riegel,stiel,element3,element4,delta_t):
    # Die Schubsteifigkeit wird immer vom schwächsten Element bestimmt - siehe Formel

    d1=buildingProp.tube_d1
    d2=buildingProp.tube_d2
    d3=buildingProp.tube_d3
    d4=buildingProp.tube_d4
    
    if d4 == max(d1,d2,d3,d4) or d2 == max(d1,d2,d3,d4):  # Limitierendes Profil wird verstärkt
        riegel.t= [element+delta_t for element in riegel.t]

    else:
        stiel.t= [element+delta_t for element in stiel.t]


def bendingStiffnessModification (buildingProp,riegel,stiel,element3,element4,delta_t):

    stiel.t= [element+delta_t for element in stiel.t]


def interiaMomentModification (buildingProp,t_riegel,t_stiel,delta_t):

    d1=buildingProp.tube_d1
    d2=buildingProp.tube_d2
    d3=buildingProp.tube_d3
    d4=buildingProp.tube_d4
    tube_wirksam=buildingProp.tube_wirksam
    
    if tube_wirksam < buildingProp.tube_wirksam_min:      # G_w erhöhen damit Shear-Lag rduziert
        if d4 == max(d1,d2,d3,d4) or d2 == max(d1,d2,d3,d4):  # Limitierendes Profil wird verstärkt
            t_riegel=t_riegel+delta_t

        else:
            t_stiel=t_stiel+delta_t

    else:
        t_stiel=t_stiel+delta_t

    return t_riegel, t_stiel


def design(buildingProp,loads,materialProp):
    # Elemente:
    b_raster=buildingProp.b_raster
    A=b_raster**2

    A_stiele=buildingProp.b_total**2-(buildingProp.b_total-b_raster)**2
    n_stiele=buildingProp.n_stiele
    
    riegel=building.elements(A,0,'Riegel','Rechteckiges Hohlprofil')
    stiel=building.elements(A_stiele/(16*n_stiele),b_raster/n_stiele,'Stiel','Quadratisches Hohlprofil')
    innenStütze=building.elements(A,0,'Stütze','Quadratisches Hohlprofil')

    stielOhnePFH=building.elements(A_stiele/(16*n_stiele),b_raster/n_stiele,'Stütze','Quadratisches Hohlprofil')

    str_= str_framedTube

    # Tragfähigkeitsnachweise:
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(stielOhnePFH,buildingProp,loads,materialProp,str_)

    calculations.calcElementWidth(riegel,buildingProp,loads,materialProp,str_)

    # calculations.calcWeight(buildingProp,materialProp,riegel)
    # riegel.G=[element*buildingProp.b_raster/n_stiele/buildingProp.h_geschoss for element in riegel.G]
    # buildingProp.G_riegel=riegel.G

    buildingProp.riegel=riegel

    calculations.calcElementWidth(stiel,buildingProp,loads,materialProp,str_)

    riegel.t=buildingProp.riegel.t

    tr0=riegel.t[-1]
    ts0=stiel.t[-1]
   
    # Gebäudenachweise:
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,riegel,stiel)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,riegel,stiel)

    if riegel.t[-1] > tr0 or stiel.t[-1] > ts0:
        buildingProp.NW_maßgebend='Verformung'
    
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'
    
    # Massen:
    calculations.calcWeight(buildingProp,materialProp,riegel)
    riegel.G=[element*buildingProp.b_raster/n_stiele/buildingProp.h_geschoss for element in riegel.G]

    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,stielOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,stiel)

    buildingProp.G_innenStützen=[element*9 for element in innenStütze.G]
    buildingProp.G_aussteifung=[element*16*n_stiele for element in riegel.G]

    buildingProp.G_außenStützen=[element*16*n_stiele for element in stiel.G]
    buildingProp.G_total=[]
    buildingProp.G_totalOhnePFH=[]
  
    for i in range (0,len(riegel.G)):
        buildingProp.G_total.append(9*innenStütze.G[i]+16*n_stiele*riegel.G[i]+16*n_stiele*stiel.G[i])
        buildingProp.G_totalOhnePFH.append(9*innenStütze.G[i]+16*n_stiele*stielOhnePFH.G[i])

    # Eigenfrequenz
    calculations.calcEigenfrequency(buildingProp,loads,materialProp)

    buildingProp.t_innenStützen=innenStütze.t
    buildingProp.t_randStützen=stiel.t
    buildingProp.t_eckStützen=stiel.t
    buildingProp.t_kern=[]
    buildingProp.t_diagonale=[]
    buildingProp.t_querstrebe=[]
    buildingProp.t_riegel=riegel.t