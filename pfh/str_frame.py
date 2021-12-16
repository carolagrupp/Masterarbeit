# ------------------------------------------------------------------------------
# Description:  Berechnung Rahmentragwerk
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-16      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
from pfh import building
from pfh import calculations
from pfh import str_frame
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    b_raster=buildingProp.b_raster
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

    #M_max und N_max für aktuelles Geschoss 
    if element.typ=='Riegel':   #Wind und Vertikallasten
        M_max=((Gd+Qd)*b_raster+element.A*gamma)*b_raster**2*0.107 + 1/5*1/4*h_geschoss/2*gamma_w*V[s]     # q*l*0.107 plus Stützenmoment, welches sich aus der Horizontalkraft der Innenstützen (1/4 der Querkraft pro Rahmen) und der halben Geschosshöhe zusammensetzt
        M_kombi=0
        N_max=0
        N_kombi=0


    if element.typ=='Stiel_innen':
        M_max= 1/4*h_geschoss/2*gamma_w*V[s]/5      
        M_kombi=Psi_w*1/4*h_geschoss/2*gamma_w*V[s]/5

        if  s==n and x*n_abschnitt < n: #unvollständiger Abschnitt
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            Ng_riegel=buildingProp.G_riegel[-1]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)  # Faktor am Ende sorgt dafür dass für alle Stützentypen Riegellänge passt

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g   #EG aus aktuellem Abschnitt (noch nicht in Ng_darüberliegend)
            i=int(s/n_abschnitt-1)      #-1, da Stellen ab 0 existieren
            Ng_riegel=buildingProp.G_riegel[i]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)   #EG Riegel
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+Ng_riegel     # Gewicht von Riegeln lastet je Geschoss zusätzlich auf Stützen 
        N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s+Ng_riegel   # Bei N_max wird q als 1. veränderl. Einwirkung angesehen, Wind wird abgemindert


    if element.typ=='Stiel_außen':
        M_max=((Gd+Psi_q*Qd)*A_einzug*2/b_raster+element.A*gamma)*b_raster**2/12 + 1/8*h_geschoss/2*gamma_w*V[s]/5      #aus Streckenlast + aus H-Last, Kombinationsbeiwert q
        M_kombi=((Gd+Qd)*A_einzug*2/b_raster+element.A*gamma)*b_raster**2/12 + Psi_w*1/8*h_geschoss/2*gamma_w*V[s]/5    #Kombinationsbeiwert Wind
            
        M_kräftepaar=M[s]/5-h_geschoss/2*V[s]/5         # Anteil des Moments der nicht über ein Moment in den Stützen sondern über ein Kräftepaar abgetragen wird
        F=gamma_w*M_kräftepaar/(b_raster*4)                          # Kraft in den Außenstützen aus Kräftepaar, wenn Außenstützen dreimal so viel Kraft aufnehmen wie Innenstützen

        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            Ng_riegel=buildingProp.G_riegel[-1]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)
            Ng_riegel=buildingProp.G_riegel[i]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+Ng_riegel+Psi_w*F
        N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s+Ng_riegel+F


    if element.typ=='Riegel ohne PFH':  #nur Vertikallasten
        M_max=((Gd+Qd)*b_raster+element.A*gamma)*b_raster**2*0.107    # q*l*0.107 plus Stützenmoment, welches sich aus der Horizontalkraft der Innenstützen (1/4 der Querkraft pro Rahmen) und der halben Geschosshöhe zusammensetzt
        M_kombi=0
        N_max=0
        N_kombi=0


    if element.typ=='Stiel_innen ohne PFH':
        M_max=0
        M_kombi=0
            
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            Ng_riegel=buildingProp.G_riegel[-1]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)
            Ng_riegel=buildingProp.G_riegel[i]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+Ng_riegel
        N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s+Ng_riegel
    

    if element.typ=='Stiel_außen ohne PFH':
        M_max=((Gd+Psi_q*Qd)*A_einzug*2/b_raster+element.A*gamma)*b_raster**2/12
        M_kombi=((Gd+Qd)*A_einzug*2/b_raster+element.A*gamma)*b_raster**2/12
            
        if  s==n and x*n_abschnitt < n:
            Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g
            Ng_riegel=buildingProp.G_riegel[-1]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)

        else:
            Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g
            i=int(s/n_abschnitt-1)
            Ng_riegel=buildingProp.G_riegel[i]*gamma_g*10*(2*A_einzug/b_raster**2+0.5*l_fassade/b_raster)
        
        N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s+Ng_riegel
        N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s+Ng_riegel   # Bei N_max wird q als 1. veränderl. Einwirkung angesehen, Wind wird abgemindert


    return N_max,N_kombi,M_max,M_kombi


def buildingStiffness(buildingProp,materialProp,riegel,innenStütze,randStütze,eckStütze):
    # I und ggf. A berechnen für EI und GA
    I=[]
    GA=[]

    for i in range (0,len(riegel.t)):
        calculations.calcProfileProp(riegel,buildingProp,materialProp,riegel.t[i])
        calculations.calcProfileProp(innenStütze,buildingProp,materialProp,innenStütze.t[i])
        calculations.calcProfileProp(randStütze,buildingProp,materialProp,randStütze.t[i])
        calculations.calcProfileProp(eckStütze,buildingProp,materialProp,eckStütze.t[i])

        I_i= (9*innenStütze.I+12*randStütze.I+4*eckStütze.I                                                 # Trägheitsmomente aller Stützen
            +3*(2*innenStütze.A*buildingProp.b_raster**2+2*randStütze.A*(buildingProp.b_raster*2)**2)       # Steiner Anteil innere Rahmen
            +2*(2*randStütze.A*buildingProp.b_raster**2+2*eckStütze.A*(buildingProp.b_raster*2)**2))        # Steiner Anteil äußere Rahmen
        I.append(I_i)
        
        GA_riegelEI=5*4*12*materialProp.E*riegel.I/(buildingProp.b_raster*buildingProp.h_geschoss)
        GA_stielEI=12*materialProp.E/(buildingProp.h_geschoss**2)*(9*innenStütze.I+12*randStütze.I+4*eckStütze.I)     # GA_stiel=Summe(12*EI/h^2)=12*E/h^2*Summe(I)
        GA_riegelGA=4*5*(materialProp.G*riegel.A*2/3*5/6)/(buildingProp.h_geschoss/buildingProp.b_raster)
        GA_stielGA=(9*innenStütze.A+12*randStütze.A+4*eckStütze.A)*(materialProp.G*0.5*5/6)

        GA_i=(GA_riegelEI**-1+GA_stielEI**-1+GA_riegelGA**-1+GA_stielGA**-1)**-1
        GA.append(GA_i)
    
    buildingProp.I=I
    buildingProp.GA=GA
    buildingProp.GA_riegelEI=GA_riegelEI
    buildingProp.GA_stielEI=GA_stielEI
    buildingProp.GA_riegelGA=GA_riegelGA
    buildingProp.GA_stielGA=GA_stielGA

def shearStiffnessModification (buildingProp,riegel,innenStütze,randStütze,eckStütze,delta_t):
    # Die Schubsteifigkeit wird immer vom schwächsten Element bestimmt - siehe Formel

    #wenn Stiel kleinere GAs verursacht
    if buildingProp.GA_stielEI <= min(buildingProp.GA_riegelEI,buildingProp.GA_stielGA,buildingProp.GA_riegelGA) or buildingProp.GA_stielGA <= min(buildingProp.GA_riegelEI,buildingProp.GA_stielEI,buildingProp.GA_riegelGA):                      
        innenStütze.t= [element+delta_t for element in innenStütze.t]
        randStütze.t= [element+delta_t for element in randStütze.t]
        eckStütze.t= [element+delta_t for element in eckStütze.t]

    #wenn Riegel vergrößert werden muss
    else:
        riegel.t= [element+delta_t for element in riegel.t]


def bendingStiffnessModification (buildingProp,riegel,innenStütze,randStütze,eckStütze,delta_t):
    
    innenStütze.t= [element+delta_t for element in innenStütze.t]
    randStütze.t= [element+delta_t for element in randStütze.t]
    eckStütze.t= [element+delta_t for element in eckStütze.t]


def design(buildingProp,loads,materialProp):
    # Elemente:
    b_raster=buildingProp.b_raster
    A=b_raster**2
    
    riegel=building.elements(A,0,'Riegel','Rechteckiges Hohlprofil')
    #stiel=building.elements(A,0,'Stiel_innen','Quadratisches Hohlprofil')
    #stiel_außen=building.elements(A/2,b_raster,'Stiel_außen','Quadratisches Hohlprofil')
    innenStütze=building.elements(A,0,'Stiel_innen','Quadratisches Hohlprofil')       # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStütze=building.elements(1/2*A,b_raster,'Stiel_innen','Quadratisches Hohlprofil')#Randstütze in innerem Riegel
    randStütze_außen=building.elements(1/2*A,b_raster,'Stiel_außen','Quadratisches Hohlprofil')#Randstütze in äußerem Riegel
    eckStütze=building.elements(1/4*A,b_raster,'Stiel_außen','Quadratisches Hohlprofil')

    riegelOhnePFH=building.elements(A,0,'Riegel ohne PFH','Rechteckiges Hohlprofil')
    innenStützeOhnePFH=building.elements(A,0,'Stiel_innen ohne PFH','Quadratisches Hohlprofil')       # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStützeOhnePFH=building.elements(1/2*A,b_raster,'Stiel_außen ohne PFH','Quadratisches Hohlprofil')
    eckStützeOhnePFH=building.elements(1/4*A,b_raster,'Stiel_außen ohne PFH','Quadratisches Hohlprofil')

    str_= str_frame

    # Tragfähigkeitsnachweise:
    calculations.calcElementWidth(riegel,buildingProp,loads,materialProp,str_)

    calculations.calcWeight(buildingProp,materialProp,riegel)   #erstellt riegel.G in calcWeight, Abschnittsweise Liste = EG des Querschnitts aller Geschosse darüber
    riegel.G=[element*buildingProp.b_raster/buildingProp.h_geschoss for element in riegel.G]    #Umrechnen da Länge b_raster und nicht h_geschoss
    buildingProp.G_riegel=riegel.G

    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStütze_außen,buildingProp,loads,materialProp,str_)

    for i in range (0,len(randStütze.t)):
        randStütze.t[i]=max(randStütze.t[i],randStütze_außen.t[i])

    calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)

    #calculations.calcElementWidth(stiel,buildingProp,loads,materialProp,str_)
    #calculations.calcElementWidth(stiel_außen,buildingProp,loads,materialProp,str_)

    #for i in range (0,len(stiel.t)):
    #    stiel.t[i]=max(stiel.t[i],stiel_außen.t[i])

    calculations.calcElementWidth(riegelOhnePFH,buildingProp,loads,materialProp,str_)

    calculations.calcWeight(buildingProp,materialProp,riegelOhnePFH)
    riegelOhnePFH.G=[element*buildingProp.b_raster/buildingProp.h_geschoss for element in riegelOhnePFH.G]
    buildingProp.G_riegelOhnePFH=riegelOhnePFH.G

    calculations.calcElementWidth(innenStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStützeOhnePFH,buildingProp,loads,materialProp,str_)

    #calculations.calcElementWidth(stielOhnePFH,buildingProp,loads,materialProp,str_)
    #calculations.calcElementWidth(stielOhnePFH_außen,buildingProp,loads,materialProp,str_)
    
    ts0=randStütze.t[-1]    #QS der Randstütze im untersten Abschnitt
    tr0=riegel.t[-1]        #QS des Riegels im untersten Abschnitt

    #for i in range (0,len(stielOhnePFH.t)):
    #    stielOhnePFH.t[i]=max(stielOhnePFH.t[i],stielOhnePFH_außen.t[i])    
    
    # Gebäudenachweise:
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,riegel,innenStütze,randStütze,eckStütze)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,riegel,innenStütze,randStütze,eckStütze)
    
    if randStütze.t[-1] > ts0 or riegel.t[-1] > tr0:    #Vergößerung des QSs in Verformungsnachweisen
        buildingProp.NW_maßgebend='Verformung'
    
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'

    # Massen:
    calculations.calcWeight(buildingProp,materialProp,riegel)
    riegel.G=[element*buildingProp.b_raster/buildingProp.h_geschoss for element in riegel.G]
    #calculations.calcWeight(buildingProp,materialProp,stiel)
    #calculations.calcWeight(buildingProp,materialProp,stielOhnePFH)

    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,randStütze)
    calculations.calcWeight(buildingProp,materialProp,eckStütze)

    calculations.calcWeight(buildingProp,materialProp,innenStützeOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,randStützeOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,eckStützeOhnePFH)

    buildingProp.G_innenStützen=[element*9 for element in innenStütze.G]    #Gesamtgewicht je Abschnitt
    buildingProp.G_aussteifung=[element*4*5*2 for element in riegel.G]      #Gesamtgewicht der Riegel je Abschnitt

    buildingProp.G_außenStützen=[]
    buildingProp.G_total=[]
    buildingProp.G_totalOhnePFH=[]

    #for i in range (0,len(stiel.G)):
     #   buildingProp.G_total.append(25*stiel.G[i]+40*riegel.G[i])
      #  buildingProp.G_totalOhnePFH.append(25*stielOhnePFH.G[i]+40*riegelOhnePFH.G[i])
    
    for i in range (0,len(riegel.G)):
        buildingProp.G_außenStützen.append(12*randStütze.G[i]+4*eckStütze.G[i])
        buildingProp.G_total.append(9*innenStütze.G[i]+40*riegel.G[i]+12*randStütze.G[i]+4*eckStütze.G[i])
        buildingProp.G_totalOhnePFH.append(9*innenStützeOhnePFH.G[i]+40*riegelOhnePFH.G[i]+12*randStützeOhnePFH.G[i]+4*eckStützeOhnePFH.G[i])
    

    # Eigenfrequenz
    calculations.calcEigenfrequency(buildingProp,loads,materialProp)

    buildingProp.t_innenStützen=innenStütze.t
    buildingProp.t_randStützen=randStütze.t
    buildingProp.t_eckStützen=eckStütze.t
    buildingProp.t_kern=[]
    buildingProp.t_diagonale=[]
    buildingProp.t_querstrebe=[]
    buildingProp.t_riegel=riegel.t