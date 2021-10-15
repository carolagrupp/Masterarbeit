# ------------------------------------------------------------------------------
# Description:  Manuelle Berechnungen der Elementeigenschaften, des Nachweis der 
#               Tragfähigkeit, der Schubverformung
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-09      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
from pfh import fea

# ------------------------------------------------------------------------------

def calcProfileProp (element,buildingProp,materialProp,t):
    b_kern=buildingProp.b_raster*2
    minderung_A=materialProp.minderung_A
    minderung_I=materialProp.minderung_I
    verhältnis_td=materialProp.verhältnis_td

    if element.profil=='Vollprofil':
        element.A=(t/100)*(t/100)
        element.A_min=element.A
        element.I=(t/100)**4/12
        element.W=element.I*2/(t/100)
    
    if element.profil=='Kern':
        element.A=b_kern*t/100*4
        element.A_min=minderung_A*element.A
        element.I=((b_kern+t/100)**4/12 - (b_kern-t/100)**4/12)*minderung_I
        element.W=element.I*2/(b_kern+t/100)

    if element.profil=='Rechteckiges Hohlprofil':
        h=t/100         # h und b=t/200 sind Außenmaße nicht Achsmaße wie beim Kern
        b=t/200
        d=h*verhältnis_td
        element.A=h*b-(h-d*2)*(b-d*2)
        element.A_min=element.A
        element.I=h**3*b/12-(h-d*2)**3*(b-d*2)/12
        element.W=element.I*2/h
                            
    if element.profil=='Quadratisches Hohlprofil':
        h=t/100         # h und b=h sind Außenmaße
        d=t/100*verhältnis_td
        element.A=h**2-(h-d*2)**2
        element.A_min=element.A
        element.I=h**4/12-(h-d*2)**4/12
        element.W=element.I*2/h  


def calcElementWidth(element,buildingProp,loads,materialProp,str_): 
    'Nachweis der Tragfähigkeit der Elemente'
    # Extrahieren der benötigten Eigenschaften aus den Objekten:   
    n=buildingProp.n
    n_abschnitt=buildingProp.n_abschnitt
    x=buildingProp.x
    h_geschoss=buildingProp.h_geschoss
        
    f=materialProp.f
    gamma=materialProp.gamma

    if loads.alpha_a=="Nein" or element.A_einzug == 0:
        alpha_a=1
    
    else:   
        alpha_a=0.5+10/element.A_einzug

    #Startwerte:
    t=materialProp.t_min
    a=materialProp.delta_t
    b=[]
    Ng_darüberliegend=0

    if x*n_abschnitt < n:
        z=x+2
    
    else:
        z=x+1

    #Äußere Schleife: Wiederholen für alle Abschnitte
    for i in range(1,z):

        if i==x+1:
            s=n        # s für Stelle an der berechnet wird
        
        else: 
            s=n_abschnitt*i
        
        sigma=2*f

        if loads.alpha_n=="Nein":
            alpha_n=1

        else:
            alpha_n=0.7+0.6/(n_abschnitt*i)
        
        alpha=min(alpha_a,alpha_n)

        #Beginn der inneren Schleife: Ermitteln von t in dem jeweiligen Abschnitt
        while sigma > f:
            calcProfileProp(element,buildingProp,materialProp,t)
            N_max,N_kombi,M_max,M_kombi=str_.calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t)
                
            sigma_Nmax=M_kombi/element.W+N_max/element.A_min     # Abdeckung beider Kombinationen aus den veränderlichen Lasten durch Wind und Verkehr
            sigma_Mmax=M_max/element.W+N_kombi/element.A_min
            
            sigma=max(sigma_Nmax,sigma_Mmax)
         
            if sigma>f:
                if buildingProp.tragwerk=='Framed Tube' and element.typ=='Stiel':
                    buildingProp.riegel.t[i-1],t = str_.interiaMomentModification(buildingProp,buildingProp.riegel.t[i-1],t,a)

                else:
                    t=t+a
        
        if i==x+1:
            Ng_darüberliegend=Ng_darüberliegend+element.A*gamma*(n-x*n_abschnitt)*h_geschoss*1.35

        else:
            Ng_darüberliegend=Ng_darüberliegend+element.A*gamma*n_abschnitt*h_geschoss*1.35
       
        b.append(t)             # t in Liste übergeben
      
    element.t=b


def calcShearDeformation(buildingProp,loads):
    'Schubverfromung berechnen'
    # Berechnung der Verformung erfolgt nicht wie sonst von oben nach unten sondern von unten nach oben
    M=loads.M
    x=buildingProp.x
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    
    w_GA_i=0
    w_GA=[0]
    l=1
    a=1
    
    if x*n_abschnitt<n:
        GA=buildingProp.GA[-a]
        
        for j in range (0,n-x*n_abschnitt):
            w_GA_i=w_GA_i+(M[-l]-M[-(l+1)])/GA
            w_GA.insert(0,w_GA_i)
            l=l+1
                    
        a=a+1
       
    for i in range (0,x):
        GA=buildingProp.GA[-a]
        
        for j in range (0,n_abschnitt):
            w_GA_i=w_GA_i+(M[-l]-M[-(l+1)])/GA
            w_GA.insert(0,w_GA_i)
            l=l+1
        
        a=a+1
       
    return w_GA


def buildingDeflection(buildingProp,loads,materialProp,str_,element1,element2=None,element3=None,element4=None):
    'Nachweis der Verformung'
    
    w = 2*loads.w_max
    delta_t=materialProp.delta_t
    buildingProp.mue=[]

    for i in range (0,len(element1.t)):
        buildingProp.mue.append(0)

    while w > loads.w_max:
        str_.buildingStiffness(buildingProp,materialProp,element1,element2,element3,element4)        
        
        # Maximale Verformung berechnen
        feModel=fea.feModel(buildingProp,loads,materialProp)
        w_EI = fea.feModel.calcStaticWindloadDeflection(feModel)
        
        w_GA=calcShearDeformation(buildingProp, loads)
               
        w = w_EI[0] + w_GA[0]
        
        if element2 == None and w > loads.w_max:
            element1.t= [element+delta_t for element in element1.t]

        if element2 != None and w > loads.w_max and w_EI[0] > loads.w_verhältnis*w_GA[0]:      # Biegeverformung größer als erwünscht
            str_.bendingStiffnessModification(buildingProp,element1,element2,element3,element4,delta_t)
        
        if element2 != None and w > loads.w_max and w_EI[0] <= loads.w_verhältnis*w_GA[0]:       # Schubverfromung größer als erwünscht
            str_.shearStiffnessModification(buildingProp,element1,element2,element3,element4,delta_t)
           
    buildingProp.w_EI=w_EI # Wert an Spitze
    buildingProp.w_GA=w_GA # Liste mit Werten an jedem Geschoss

    buildingProp.w=w


def interstoryDrift(buildingProp,loads,materialProp,str_,element1,element2=None,element3=None,element4=None):
    'Nachweis interstory Drift'
    
    Teta_i=2*loads.maxTeta_i
    delta_t=materialProp.delta_t
    buildingProp.mue=[]

    for i in range (0,len(element1.t)):
        buildingProp.mue.append(0)

    while Teta_i > loads.maxTeta_i:
        str_.buildingStiffness(buildingProp,materialProp,element1,element2,element3,element4)
        
        # Knotenverformungen und Verformungsdifferenz berechnen
        feModel=fea.feModel(buildingProp,loads,materialProp)                    # Modell zu Ermitllung der Biegeverformung jedes Knotens
        buildingProp.w_GA=calcShearDeformation(buildingProp, loads)             # Funktion zur Ermittlung des Schubverformung jedes Knotens
        Teta_i, w_EI_max = fea.feModel.calcInterstoryDrift(feModel,buildingProp)          # Funktion zur Ermitllung der Verfromungsdifferenz der Knoten aus Schub und Biegung
        
        #if Teta_i > loads.maxTeta_i:
         #   element1.t = [element+2.5 if element == 27.5 else element+5 for element in element1.t]
          #  if element2 != None:
           #     element2.t= [element+2.5 if element == 27.5 else element+5 for element in element2.t]

        if element2 == None and Teta_i > loads.maxTeta_i:
            element1.t= [element+delta_t for element in element1.t]

        if element2 != None and Teta_i > loads.maxTeta_i and w_EI_max > loads.w_verhältnis*buildingProp.w_GA[0]:      # Biegeverformung an Gebäudespitze größer als erwünscht
            str_.bendingStiffnessModification(buildingProp,element1,element2,element3,element4,delta_t)
        
        if element2 != None and Teta_i > loads.maxTeta_i and w_EI_max <= loads.w_verhältnis*buildingProp.w_GA[0]:     # Schubverfromung an Gebäudespitze größer als erwünscht
            str_.shearStiffnessModification(buildingProp,element1,element2,element3,element4,delta_t)


def calcWeight(buildingProp, materialProp, element):
    element.G=[]
    G=0

    for i in range (0,len(element.t)):
        calcProfileProp(element,buildingProp,materialProp,element.t[i])
        G=G+element.A*materialProp.gamma/10*buildingProp.n_abschnitt*buildingProp.h_geschoss
        element.G.append (G)

    if buildingProp.n_abschnitt*buildingProp.x < buildingProp.n:
        element.G[-1]=element.G[-1]-element.A*materialProp.gamma/10*((buildingProp.x+1)*buildingProp.n_abschnitt-buildingProp.n)*buildingProp.h_geschoss


def calcEigenfrequency(buildingProp,loads,materialProp):

    buildingProp.mue=[]

    mue_konstant=((loads.gamma_gdyn*(loads.gk1+loads.gk2)+loads.gamma_qdyn*loads.qk)*buildingProp.b_total**2/buildingProp.h_geschoss
                +loads.gk_fassade*loads.gamma_gdyn*buildingProp.b_total*4/buildingProp.h_geschoss)/10     # Decken und Nutzlasten auf die Höhe aufgeteilt "/10" für Umrechnung von kN auf t

    for i in range (0,len(buildingProp.G_total)):  
        
        if i == 0:
            mue_tragwerk=buildingProp.G_total[i]*loads.gamma_gdyn/(buildingProp.n_abschnitt*buildingProp.h_geschoss)

        elif buildingProp.n_abschnitt*i > buildingProp.n:
            mue_tragwerk=(buildingProp.G_total[i]-buildingProp.G_total[i-1])*loads.gamma_gdyn/((buildingProp.n-buildingProp.n_abschnitt*(i-1))*buildingProp.h_geschoss)

        else:  
            mue_tragwerk=(buildingProp.G_total[i]-buildingProp.G_total[i-1])*loads.gamma_gdyn/(buildingProp.n_abschnitt*buildingProp.h_geschoss)

        buildingProp.mue.append(mue_konstant+mue_tragwerk)

    feModel=fea.feModel(buildingProp,loads,materialProp)
    fq_e=fea.feModel.getEigenfrequency(feModel)

    buildingProp.eigenFrequenz_EI=fq_e

    if buildingProp.w_EI[0]==0:
        buildingProp.eigenFrequenz=fq_e

    else:
        buildingProp.eigenFrequenz=fq_e*(1/(1+buildingProp.w_GA[0]/buildingProp.w_EI[0]))**0.5

    print('Eigenfrequenz EI =',fq_e)
    print('Eigenfrequenz genähert=',buildingProp.eigenFrequenz)
