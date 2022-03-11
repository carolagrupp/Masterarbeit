# ------------------------------------------------------------------------------
# Description:  Manuelle Berechnungen der Elementeigenschaften, des Nachweis der 
#               Tragfähigkeit, der Schubverformung
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
from pfh import fea
#from WindData.data import DataProp
#from pfh.building import loads
from WindData import response

# ------------------------------------------------------------------------------

def calcProfileProp (element,buildingProp,materialProp,t):
    #für calcElementWidth und calcWeight
    b_kern=buildingProp.b_raster*2          #in m
    minderung_A=materialProp.minderung_A    #in app.py submit_materialProp ist Eingabe in GUI Profil
    minderung_I=materialProp.minderung_I
    verhältnis_td=materialProp.verhältnis_td
    h_geschoss = buildingProp.h_geschoss
    

    if element.profil=='Vollprofil':
        element.A=(t/100)*(t/100)
        element.A_min=element.A
        element.I=(t/100)**4/12
        element.W=element.I*2/(t/100)
    
    elif element.profil=='Kern': #in str_core
        element.A=(b_kern*t/100)*4        #in m, durch 100 da t in cm
        element.A_min=minderung_A*element.A
        element.I=((b_kern+t/100)**4/12 - (b_kern-t/100)**4/12)*minderung_I #Volles Quadrat - inneres Quadrat
        element.W=element.I*2/(b_kern+t/100)    #I/Hebelarm

    elif element.profil=='Rechteckiges Hohlprofil':
        h=t/100         # h und b=t/200 sind Außenmaße nicht Achsmaße wie beim Kern
        b=t/(2*100)     #Breite = halbe Höhe t
        d=h*verhältnis_td
        element.A=h*b-(h-d*2)*(b-d*2)
        element.A_min=element.A
        element.I=h**3*b/12-(h-d*2)**3*(b-d*2)/12
        element.W=element.I*2/h
                            
    elif element.profil=='Quadratisches Hohlprofil':
        h=t/100         # h und b=h sind Außenmaße
        d=t/100*verhältnis_td
        element.A=h**2-(h-d*2)**2
        element.A_min=element.A
        element.I=h**4/12-(h-d*2)**4/12
        element.W=element.I*2/h  

    elif element.profil == 'Rechteckiges Vollprofil':
        h = h_geschoss
        b = t/100
        element.A = h*b
        element.A_min = minderung_A*element.A
        element.I = b*h**3/12
        element.W = element.I*2/h
        



def calcElementWidth(element,buildingProp,loads,materialProp,str_): 
    'Nachweis der Tragfähigkeit der Elemente'
    # in str_.design für alle Tragwerkstypen (Innenstütze, ...)
    # Extrahieren der benötigten Eigenschaften aus den Objekten:   
    n=buildingProp.n
    n_abschnitt=buildingProp.n_abschnitt
    x=buildingProp.x        #def. in app.py submit_buildingProp:int(n/n_abschnitt)
    h_geschoss=buildingProp.h_geschoss
        
    f=materialProp.f        #Festigkeit Material
    gamma=materialProp.gamma
    gamma_g=loads.gamma_g

    #Abminderungsbeiwert Einzugsfläche (EC1-1-1/NA, 6.3.1.2(10))
    if loads.alpha_a=="Nein" or element.A_einzug == 0:
        alpha_a=1
    
    else:   
        alpha_a=0.5+10/element.A_einzug     #A_einzug angegeben in str_.design
    buildingProp.alpha_a = alpha_a

    #Startwerte:
    t=materialProp.t_min
    a=materialProp.delta_t
    b=[]
    Ng_darüberliegend=0

    if x*n_abschnitt < n:   #ein weiterer Iterationsschritt nötig, da x für Ganzzahl abgerundet wurde
        z=x+2
        
    
    else:               #wenn x*n_abschnitt = n
        z=x+1

    buildingProp.z = z
    buildingProp.j = 0        #Iterationsparameter für Outrigger

    if element.typ == 'Kern':
        buildingProp.sigma_kern = []
    elif element.typGenau == 'innenStütze':
        buildingProp.sigma_innenStütze = []
    elif element.typGenau == 'randStütze':
        buildingProp.sigma_randStütze = []
    elif element.typGenau == 'eckStütze':
        buildingProp.sigma_eckStütze = []
     
    
    
    #Äußere Schleife: Wiederholen für alle Abschnitte
    for i in range(1,z):

        if i==x+1:      # nur wenn z=x+2
            s=n        # s für Stelle an der berechnet wird
        
        else: # s Geschosse über aktuellem Abschnitt
            s=n_abschnitt*i
        
        sigma=2*f   # Anfangswert, damit sigma für Schleife größer als f ist

        if loads.alpha_n=="Nein":
            alpha_n=1

        else:
            alpha_n=0.7+0.6/(n_abschnitt*i)
        
        alpha=min(alpha_a,alpha_n)      # Nach Norm darf nur einer der beiden Werte angewendet werden
        
        buildingProp.i_aktuell = i
        buildingProp.ErsteBerechnungSigma = True
        
        if buildingProp.tragwerk == 'Outrigger' and buildingProp.Iteration == True and element.typ == 'Stütze': # Dann t bereits vorhanden
            t = element.t[i-1]

        # Beginn der inneren Schleife: Ermitteln von t in dem jeweiligen Abschnitt
        while sigma > f:
            calcProfileProp(element,buildingProp,materialProp,t)        #A, I und W für aktuelles t
            if buildingProp.tragwerk == 'Outrigger' and buildingProp.Iteration == True:
                buildingProp.outriggerAbschnitt = False
                buildingProp.anzahlOutInAbschnitt = 0
                for k in buildingProp.posOut_abschnitt:
                    if i == k+1:
                        buildingProp.anzahlOutInAbschnitt += 1
                        buildingProp.outriggerAbschnitt = True
            N_max,N_kombi,M_max,M_kombi=str_.calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t)
            
            
            sigma_Nmax=M_kombi/element.W+N_max/element.A_min     # Abdeckung beider Kombinationen aus den veränderlichen Lasten durch Wind und Verkehr
            sigma_Mmax=M_max/element.W+N_kombi/element.A_min
            
            sigma=max(sigma_Nmax,sigma_Mmax)
         
            if sigma>f:     #Querschnittsvergrößerung
                if buildingProp.tragwerk=='Framed Tube' and element.typ=='Stiel':
                    buildingProp.riegel.t[i-1],t = str_.interiaMomentModification(buildingProp,buildingProp.riegel.t[i-1],t,a)
                    
                else:
                    t=t+a   #Erhöhen um delta_t
                buildingProp.ErsteBerechnungSigma = False
        # Belastung der Elements speichern
        if buildingProp.tragwerk == 'Outrigger':
            if element.typ == 'Kern':
                buildingProp.sigma_kern.append(sigma)
            elif element.typGenau == 'innenStütze':
                buildingProp.sigma_innenStütze.append(sigma)
            elif element.typGenau == 'randStütze':
                buildingProp.sigma_randStütze.append(sigma)
            elif element.typGenau == 'eckStütze':
                buildingProp.sigma_eckStütze.append(sigma)

        # Ng auf darunterliegendes Geschoss aus EG des Elements
        if i==x+1:
            Ng_darüberliegend=Ng_darüberliegend+element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g

        else:
            Ng_darüberliegend=Ng_darüberliegend+element.A*gamma*n_abschnitt*h_geschoss*gamma_g
       
        b.append(t)             # t in Liste übergeben
      
    element.t=b         #Länge = Anzahl Abschnitte


def calcDynamicElementWidth(buildingProp, loads, materialProp, DataProp, str_, element1, element2=None, element3=None, element4=None):

    loads.windData = True

    delta_t=materialProp.delta_t

    if loads.GK == 2:
            loads.alpha_v = 0.16
    elif loads.GK == 4:
            loads.alpha_v = 0.30

    Directions_belastung = [0,5,10,15,20,25,30,35,40,45]
    #for direction in Directions_belastung:
        #Durchführung aller Schritte
    schlankheit = buildingProp.schlankheit
    direction = str(Directions_belastung[0])

    loads.grabData(buildingProp, DataProp, direction, schlankheit)

    buildingProp.mue=[]

    for i in range (0,len(element1.t)): #element1.t = kern.t (siehe str_core.design)
        buildingProp.mue.append(0)                                          # oder anders, s. fea auf tik server

    achsen_response = [[DataProp.BMD_fs, DataProp.LFD, 'D'], [DataProp.BML_fs, DataProp.LFL, 'L']]      # LFD und LFL sind Kräfte auf Geschosshöhen

    for achse in [0,1]:

        # Übergabe der aktuellen Lasten
        loads.F_p_j = achsen_response[achse][1]

        w_max = loads.w_max
        
        w_tip = 2*w_max

        # hier Massenberechnung mit aktuellen ts des Tragwerks, vorher schon Bemessung auf Vertikallasten über Calcelementloads?
        # kein buildingProp.calcBuildMass()
        while w_tip > w_max:

            # Stiffnessberechnung des aktuellen Systems für Eigenfrequenz
            str_.buildingStiffness(buildingProp, materialProp, element1, element2, element3, element4)

            # Create analysis model
            feModelDyn = fea.feModel(buildingProp, loads, materialProp)

            # Compute eigenfrequencies
            feModelDyn.getEigenfrequency()

            # Calc generalized quantities (K_gen, M_gen)
            feModelDyn.calcGeneralizedQuantitites()

            # Calc response forces
            responseForces = response.responseForces(achsen_response[achse][0], DataProp.dT, feModelDyn.fq_e, buildingProp.D, feModelDyn.fq_e, 360)
            #Wieso zweimal fq_e als Input? response nimmt zweites mal als nue, weiso 360?

            # Calc response deflections
            # LFL/LFD = FLoor Forces: Wie Berechnung für HFFB Modelldaten?
            responseDeflections = response.responseDeflection(feModelDyn, responseForces, achsen_response[achse][1])
            w_tip = responseDeflections.delta_tip_r_max
            #w_GA dafür noch nicht berücksichtigt, hinzufügen?

            if w_tip > w_max:
                if element2 == None:        #zB core hat nur ein element, bei dem der QS vergrößert werden kann
                    element1.t= [element+delta_t for element in element1.t]
                #else:
                    #siehe wie buildingDeflection

        a_max = loads.a_max     #muss noch in gui hinzugefügt werden 

        a_rms = 2*a_max

        while a_rms > a_max:

            #buildingProp.calcBuildMass()    # Brauche ich das?
            # nur Stiffnessberechnung des aktuellen Systems für Eigenfrequenz
            str_.buildingStiffness(buildingProp, materialProp, element1, element2, element3, element4)

            # Create analysis model
            feModelDyn = fea.feModel(buildingProp, loads, materialProp)

            # Compute eigenfrequencies
            feModelDyn.getEigenfrequency()

            # Calc generalized quantities (K_gen, M_gen)
            feModelDyn.calcGeneralizedQuantitites()
            # Iterativ: Anfang bis Response + folgendes
            # Calc response accelerations
            responseAccelerations = response.responseAccelerations(feModelDyn, achsen_response[achse][0], DataProp.dT, feModelDyn.fq_e, buildingProp.D, feModelDyn.fq_e, 360)
            a_rms = responseAccelerations.a_r_rms

            # Erhöhen von t, wenn a größer als a_max
            if a_rms > a_max:
                if element2 == None:        #zB core hat nur ein element, bei dem der QS vergrößert werden kann
                    element1.t= [element+delta_t for element in element1.t]
                #else:
                    #siehe wie buildingDeflection

        if achsen_response[achse][2] == 'D':
            buildingProp.dyn_w_tip_D = w_tip
            buildingProp.dyn_a_rms_D = a_rms
            buildingProp.dyn_a_max_D = responseAccelerations.a_r_max

        elif achsen_response[achse][2] == 'L':
            buildingProp.dyn_w_tip_L = w_tip
            buildingProp.dyn_a_rms_L = a_rms
            buildingProp.dyn_a_max_L = responseAccelerations.a_r_max

    print(buildingProp.dyn_w_tip_D)






def calcShearDeformation(buildingProp,loads):
    'Schubverfromung berechnen'
    # Berechnung der Verformung erfolgt nicht wie sonst von oben nach unten sondern von unten nach oben
    if buildingProp.tragwerk == 'Outrigger':
        M = buildingProp.moment
    else:
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
        
        for j in range(0,n-x*n_abschnitt):
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
    #in str_.design
    w = 2*loads.w_max       #loads.w_max=h_total/Eingabe w_max, in m, 2x, damit der Startwert für die while Schleife größer als w_max ist
    delta_t=materialProp.delta_t    #Schrittweite Querschnittserhöhung
    buildingProp.mue=[]

    for i in range (0,len(element1.t)): #element1.t = kern.t (siehe str_core.design)
        buildingProp.mue.append(0)

    while w > loads.w_max:
        str_.buildingStiffness(buildingProp,materialProp,element1,element2,element3,element4)
        
        # Maximale Verformung berechnen
        feModel = fea.feModel(buildingProp,loads,materialProp)
        w_EI = fea.feModel.calcStaticWindloadDeflection(feModel)    #Angabe Liste jedes Geschoss
        
        if buildingProp.Nachweis_GZT == False:
            # Auslesen Biegemomente über die Höhe
            moment = fea.feModel.calcBendingMoment(feModel)
            
            buildingProp.moment_ges = moment                                                # kNm                 Array mit 2 Momenten je Element bzw Geschoss
            
            buildingProp.moment = [0]                                                                           # Array mit maximalem Moment je Geschoss
            for i in range(0,buildingProp.n):
                bm_max = max(abs(moment[2*i]), abs(moment[2*i+1]))
                buildingProp.moment.append(bm_max)

        w_GA = calcShearDeformation(buildingProp, loads)
               
        w = w_EI[0] + w_GA[0]
        
        if element2 == None and w > loads.w_max:        #zB core hat nur ein element, bei dem der QS vergrößert werden kann
            element1.t= [element+delta_t for element in element1.t]

        if element2 != None and w > loads.w_max and w_EI[0] > loads.w_verhältnis*w_GA[0]:      # Biegeverformung größer als erwünscht
            str_.bendingStiffnessModification(buildingProp,materialProp,element1,element2,element3,element4,delta_t)
        
        if element2 != None and w > loads.w_max and w_EI[0] <= loads.w_verhältnis*w_GA[0]:       # Schubverfromung größer als erwünscht
            str_.shearStiffnessModification(buildingProp,materialProp,element1,element2,element3,element4,delta_t)
           
    buildingProp.w_EI=w_EI # Liste mit Werten an jedem Geschoss
    buildingProp.w_GA=w_GA # Liste mit Werten an jedem Geschoss

    buildingProp.w=w


def interstoryDrift(buildingProp,loads,materialProp,str_,element1,element2=None,element3=None,element4=None):
    'Nachweis interstory Drift'
    
    Teta_i=2*loads.maxTeta_i    #2x, damit Startwert für Teta_i> als maxTeta ist für Schleife
    delta_t=materialProp.delta_t
    buildingProp.mue=[]

    for i in range (0,len(element1.t)):
        buildingProp.mue.append(0)

    while Teta_i > loads.maxTeta_i:
        str_.buildingStiffness(buildingProp,materialProp,element1,element2,element3,element4)
        
        # Knotenverformungen und Verformungsdifferenz berechnen
        feModel=fea.feModel(buildingProp,loads,materialProp)                    # Modell zu Ermitllung der Biegeverformung jedes Knotens
        if buildingProp.Nachweis_GZT == False:
            # Auslesen Biegemomente über die Höhe
            moment = fea.feModel.calcBendingMoment(feModel)
            
            buildingProp.moment_ges = moment                                                # kNm                 Array mit 2 Momenten je Element bzw Geschoss
            
            buildingProp.moment = [0]                                                                           # Array mit maximalem Moment je Geschoss
            for i in range(0,buildingProp.n):
                bm_max = max(abs(moment[2*i]), abs(moment[2*i+1]))
                buildingProp.moment.append(bm_max)

        buildingProp.w_GA=calcShearDeformation(buildingProp, loads)             # Funktion zur Ermittlung des Schubverformung jedes Knotens
        Teta_i, w_EI_max = fea.feModel.calcInterstoryDrift(feModel,buildingProp)          # Funktion zur Ermitllung der Verfromungsdifferenz der Knoten aus Schub und Biegung
        
        #if Teta_i > loads.maxTeta_i:
         #   element1.t = [element+2.5 if element == 27.5 else element+5 for element in element1.t]
          #  if element2 != None:
           #     element2.t= [element+2.5 if element == 27.5 else element+5 for element in element2.t]

        if element2 == None and Teta_i > loads.maxTeta_i:
            element1.t= [element+delta_t for element in element1.t]

        if element2 != None and Teta_i > loads.maxTeta_i and w_EI_max > loads.w_verhältnis*buildingProp.w_GA[0]:      # Biegeverformung an Gebäudespitze größer als erwünscht
            str_.bendingStiffnessModification(buildingProp,materialProp,element1,element2,element3,element4,delta_t)
        
        if element2 != None and Teta_i > loads.maxTeta_i and w_EI_max <= loads.w_verhältnis*buildingProp.w_GA[0]:     # Schubverfromung an Gebäudespitze größer als erwünscht
            str_.shearStiffnessModification(buildingProp,materialProp,element1,element2,element3,element4,delta_t)


def calcWeight(buildingProp, materialProp, element):
    element.G=[]
    G=0

    for i in range (0,len(element.t)):
        calcProfileProp(element,buildingProp,materialProp,element.t[i])
        if element.typ == 'Outrigger':
            G=G+element.A*materialProp.gamma/10*buildingProp.b_raster
        elif element.typ == 'Belt Truss':
            G=G+4*element.A*materialProp.gamma/10*buildingProp.b_total
        else:
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
