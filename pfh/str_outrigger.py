# ------------------------------------------------------------------------------
# Description:  Berechnung Kerntragwerk mit Outrigger
#
# ------------------------------------------------------------------------------
# Author:    Carola Grupp
# Created:      2021-11-25  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
from Masterarbeit.pfh.building import materialProp
from Masterarbeit.pfh.calculations import calcProfileProp
from pfh import building
from pfh import calculations
from pfh import str_outrigger
from pfh import fea
import numpy
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
    
    gamma_g=loads.gamma_g
    gamma=materialProp.gamma
    gamma_w=loads.gamma_w

    A_einzug=element.A_einzug
    l_fassade=element.l_fassade

    if  s==n and x*n_abschnitt < n: #wenn n/n_abschnitt nicht aufgeht, unten noch mehrere Geschosse mit ni<n_abschnitt
        Ng=element.A*gamma*(n-x*n_abschnitt)*h_geschoss*gamma_g

    else:
        Ng=element.A*gamma*n_abschnitt*h_geschoss*gamma_g

    N_max= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha)+l_fassade*gd_fassade)*s
    N_kombi= Ng_darüberliegend+Ng+(A_einzug*(Gd+Qd*alpha*Psi_q)+l_fassade*gd_fassade)*s

    if buildingProp.Iteration == False:
        M_max=0
    
    elif buildingProp.Iteration == True:
        if buildingProp.outriggerAbschnitt == True:  #Überprüfung ob in Outrigger Abschnitt
            
            if element.typ == "Stütze":
                j = buildingProp.j
                #Ändern da nun Moment und kein sigma mehr
                deltaN = gamma_w*buildingProp.sigmaStütze_outrigger[j]*buildingProp.A_min_aktuell
                buildingProp.N_max_delta += Psi_w*deltaN
                buildingProp.N_kombi_delta += deltaN
                buildingProp.j += 1

        if element.typ == "Kern":
            moment_core = buildingProp.moment[buildingProp.i_aktuell]
            M_max = moment_core*buildingProp.W_aktuell
                
        if element.typ == 'Stütze':
            M_max = 0
            N_max += buildingProp.N_max_delta
            N_kombi += buildingProp.N_kombi_delta
    
    M_kombi=Psi_w*M_max

    return N_max,N_kombi,M_max,M_kombi


def arrangementOutrigger(buildingProp):

    n=buildingProp.n
    n_outrigger = buildingProp.n_outrigger

    #Bestimmen der Outrigger Geometrie (nach Park, Lee et all 2016)
    if n_outrigger == 0:
        buildingProp.posOut = []

    elif n_outrigger == 1:
        buildingProp.posOut = [int(0.61*n)]

    elif n_outrigger == 2:
        buildingProp.posOut = [int(0.73*n), int(0.41*n) ]

    elif n_outrigger == 3:
        buildingProp.posOut = [int(0.80*n), int(0.55*n), int(0.33*n) ]

    elif n_outrigger == 4:
        buildingProp.posOut = [int(0.84*n), int(0.64*n), int(0.47*n), int(0.29*n) ]
    
    #Umwandlung in Abschnitte
    buildingProp.posOut_abschnitt = []
    for i in range(0, len(buildingProp.posOut)):
        buildingProp.posOut_abschnitt[i] = buildingProp.posOut[i]//buildingProp.n_abschnitt #Liste mit Liste der Anschnitte in denen Outrigger sind


def outriggerStiffness(I_outriggerExakt, buildingProp, i):
    #Berechnung der effektiven Biegesteifigkeit
    I_eff = I_outriggerExakt[i]*(1+buildingProp.b_raster/buildingProp.b_raster)**3

    return I_eff

def outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp):
    #Kopplung Outrigger und Stützen
    b_total = buildingProp.b_total
    I_effOutrigger = []
    K = []
    buildingProp.H_outrigger = []
    A_randStütze = []
    A_eckStütze = []
    t_randStützeOutrigger = []
    t_eckStützeOutrigger = []
    I_outriggerExakt = []

    #Aktualisierung Querschnitsswerte
    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        t_randStützeOutrigger.append(buildingProp.randStütze.t[n])
        t_eckStützeOutrigger.append(buildingProp.eckStütze.t[n])
        t_outrigger = buildingProp.outrigger.t[i]
        calcProfileProp(outrigger, buildingProp, materialProp, t_outrigger)
        I_outriggerExakt.append(outrigger.I)
    
    #Berechnung Steifigkeiten
    for i in range(0,len(buildingProp.posOut)):
        I_effOutrigger[i] = outriggerStiffness(I_outriggerExakt, buildingProp, materialProp,i)
        calcProfileProp(randStütze, buildingProp, materialProp, t_randStützeOutrigger[i])
        A_randStütze.append(randStütze.A_min)
        calcProfileProp(eckStütze, buildingProp, materialProp, t_eckStützeOutrigger[i])
        A_eckStütze.append(eckStütze.A_min)
    
    for i in buildingProp.posOut:
        H = i*buildingProp.h_geschoss    #Höhe des Outriggers
        buildingProp.H_outrigger.append(H)

    for i in range(0, len(buildingProp.posOut)):
        I_outrigger = 4*I_effOutrigger[i]    #4 Outrigger in Windrichtung
        A_stützen = 6*A_randStütze[i] + 4*A_eckStütze[i]
        K[i] = 12*materialProp.E*I_outrigger*A_stützen*b_total**2/(A_stützen*b_total**3+24*buildingProp.H_outrigger[i]*I_outrigger)

    return K

def calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n):
    #Outrigger
    b_total = buildingProp.b_total
    h = buildingProp.h_total
    t=kern.t[n]
    calculations.calcProfileProp(kern, buildingProp, materialProp, t)
    I_core = kern.I_min
    t = outrigger.t[i]
    calculations.calcProfileProp(outrigger, buildingProp, materialProp, t)
    I_out = outrigger.I_min
    beta_ist = materialProp.E/materialProp.E*I_core*b_total/(4*I_out*h)
    t = randStütze.t[n]
    calculations.calcProfileProp(randStütze, buildingProp, materialProp, t)
    A_rand = randStütze.A_min
    t = eckStütze.t[n]
    calculations.calcProfileProp(eckStütze, buildingProp, materialProp, t)
    A_eck = eckStütze.A_min
    alpha_ist = 2*materialProp.E*I_core/(b_total**2*materialProp.E*(6*A_rand+4*A_eck))

    return beta_ist, alpha_ist

def controllStiffnesRatios(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze):
    a=materialProp.delta_t
    Verhältnis1 = []
    Verhältnis2 = []
    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger
    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        #Berechnung Steifigkeitsverhältnis
        beta_ist, alpha_ist = calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n)
        #Runden für Optimierung (4,4 etc ist gut genug)
        beta_ist_round = round(beta_ist)
        alpha_ist_round = round(alpha_ist)
        if beta_ist_round < beta_outrigger:    #Erhöhung Kern
            kern.t[n] += a 
            Verhältnis1[n] = False
        elif beta_ist_round > beta_outrigger:    #Erhöhung Outrigger
            outrigger.t[i] += a
            Verhältnis1[n] = False
        elif alpha_ist_round < alpha_outrigger:
            kern.t[n] += a 
            Verhältnis2[n] = False
        elif alpha_ist_round > alpha_outrigger:
            eckStütze.t[n] += a 
            randStütze.t[n] += a
            Verhältnis2[n] = False

        elif beta_ist_round == beta_outrigger:
            Verhältnis1[n] = True
        elif alpha_ist_round == beta_outrigger:
            Verhältnis2[n] = True

    Test = True

    for i in buildingProp.posOut_abschnitt:
        if Verhältnis1[i] == False:
            Test = False
        if Verhältnis2[i] == False:
            Test = False

    return Test



def buildingStiffness(buildingProp,materialProp,kern,element2,element3,element4):
    #wird in calculations.buildingDeflection und interstoryDrift aufgerufen
    # I und ggf. A berechnen für EI und GA
    # anpassen? (Für fea gleich, oder?)
    I=[]
    GA=[]

    for i in range (0,len(kern.t)):     #Länge = Anzahl Abschnitte, kern.t = benötigte Querschnitsdicke nach calcElementWidth
        calculations.calcProfileProp(kern, buildingProp, materialProp, kern.t[i])
            
        I.append(kern.I)    #nur ein Wert
        GA.append(kern.A_min/2*5/6*materialProp.G)      # Halbe Kernfläche, da nur die Wände in Kraftrichtung ("Stegwände") die Schubkraft aufnehmen. Abmidnerung mit Kappa (s. Engelke S.14)

    buildingProp.I_kern=I        #Länge Anzahl Abschnitte
    buildingProp.GA_core=GA
   

def design(buildingProp,loads,materialProp):
    #in app.py.MainWondow.MainCalculation
    arrangementOutrigger(buildingProp)

    # Elemente:
    b_raster=buildingProp.b_raster
    b_total =buildingProp.b_total
    h_geschoss=buildingProp.h_geschoss
    h=buildingProp.h_total
    A=b_raster**2
    innenStütze=building.elements(A,0,'Stütze','Vollprofil')       # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStütze=building.elements(1/2*A,b_raster,'Stütze','Vollprofil')
    randStützeOhnePFH=building.elements(1/2*A,b_raster,'Stütze ohen PFH','Vollprofil')
    eckStütze=building.elements(1/4*A,b_raster,'Stütze','Vollprofil')
    eckStützeOhnePFH=building.elements(1/4*A,b_raster,'Stütze ohne PFH','Vollprofil')
    kernOhnePFH=building.elements(8*A,0,'Kern ohne PFH','Kern')
    kern=building.elements(8*A,0,'Kern','Kern')
    outrigger=building.elements(0,0,'Outrigger','Rechteckiges Vollprofil')

    str_=str_outrigger
    buildingProp.Iteration = False

    # Tragfähigkeitsnachweise Vertikallasten (Für Angangswert t):
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kernOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kern,buildingProp,loads,materialProp,str_)
    
    t0=kern.t[-1]       #Kerndicke ganz unten (letztes Element der Liste)

    #Startwert für Steifigkeiten Outriggern
    #Steifigkeitsverhältnis Kern zu Outrigger
    buildingStiffness(buildingProp,materialProp,kern)   #I des Kerns
    #buildingProp.beta_outrigger = 4  #Vielleicht über Eingabe in gui und Optimierung
    #buildingProp.alpha_outrigger = 3
    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger
    a=materialProp.delta_t
    #buildingProp.randStützeOutrigger.t =[]
    #buildingProp.eckStützeOutrigger.t =[]

    for i in buildingProp.posOut:
        #Bestimmen A_outrigger (4 Outrigger je Richtung)
        I_outrigger_i = materialProp.E/materialProp.E*buildingProp.I_kern[i//buildingProp.n_abschnitt]*b_total/(beta_outrigger*h*4)
        buildingProp.outrigger.t[i]=I_outrigger_i*12/(h_geschoss)**3*100  #in cm
        #calcProfileProp(outrigger,buildingProp,materialProp, buildingProp.outrigger.t[i])
        #buildingProp.outrigger.W.append(outrigger.W)    #Abspeichern W je Outrigger
        #buildingProp.outrigger.I.append(outrigger.I)
    for i in buildingProp.posOut_abschnitt:
        #Bestimmen A-Stützen (10 Stützen je Richtung)
        A_stütze = 2*materialProp.E/materialProp.E*buildingProp.I_kern[i]/(alpha_outrigger*b_total**2*10)
        t_stütze = numpy.sqrt(A_stütze)*100         #in cm
        while t_stütze > buildingProp.randStütze.t[i]: #Vergrößern von t wenn nötig
            buildingProp.randStütze.t[i] += a     #Ersetzen von t für ganzen Abschnitt
        #buildingProp.randStützeOutrigger.t.append(buildingProp.randStütze.t[i])
        while t_stütze > buildingProp.eckStütze.t[i]: #Vergrößern von t wenn nötig
            buildingProp.eckStütze.t[i] += a     #Ersetzen von t für ganzen Abschnitt
        #buildingProp.eckStützeOutrigger.t.append(buildingProp.eckStütze.t[i])
    
    
    Tragfähigkeitsnachweis = False
     

    while Tragfähigkeitsnachweis == False:
        buildingProp.Iteration = True

        
        #Tragfähigkeitsnachweise Aussteifung:
    
        #Federsteifigkeiten der Outriggergeschossen:
        buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)

        #Steifigkeit Kern
        buildingStiffness(buildingProp,materialProp,kern)   #nur core, passt das?

        #fea: calcBendingMoment
        feModel = fea.feModel(buildingProp, loads, materialProp)
        moment = fea.feModel.calcBendingMoment(feModel)
        buildingProp.moment = moment    #Anpassen sigma in calcelementLoad
        print(moment)

        #Auflagerkräfte, Länge müsste = 3+n_outrigger
        reactions = fea.feModel.calcreactions(feModel)
        print(reactions)

        print(reactions[1])
        print(moment[buildingProp.posOut[1]])
    
        #Auslesen Auflagerkraft an Feder
        buildingProp.moment_spring = []
        I_outrigger  = []
        for n, i in enumerate(buildingProp.posOut):
            calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
            I_outrigger[n] = outrigger.I_min

        for i in range(0, len(buildingProp.posOut)):
            moment = reactions[i]     #drei für Einspannung unten an Ende der Auflistung?
            buildingProp.moment_spring.append(moment)
            verdrehung = moment/buildingProp.K[i]
            #Kraft/Fläche (kN/m^2)
            spannung_stütze = (3*verdrehung-moment*2*materialProp.E*4*I_outrigger[i]/buildingProp.b_total)\
                *materialProp.E*buildingProp.b_total/(6*buildingProp.H_outrigger[i])

            buildingProp.sigmaStütze_outrigger.append(spannung_stütze)

        buildingProp.N_max_delta = 0 #Belastung aus Outriggern darüber
        buildingProp.N_kombi_delta = 0

        #Beanspruchung Stützen (Auflagerkraft Feder + Vertikallasten inkl aus Federn darüber)
        calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(eckStütze, buildingProp,loads,materialProp,str_)

        #Beanspruchung Kern (Beanspruchung fea + Vertikallasten)
        calculations.calcElementWidth(kern, buildingProp,loads,materialProp,str_)
        
        #Beanspruchung Outrigger
        #Moment aus Feder in fea
        f = materialProp.f
        
        for n,i in enumerate(buildingProp.posOut):
            M_max_outrigger = loads.gamma_w*reactions[n]
            while sigma_outrigger > f:
                M_g = outrigger.t[i]*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2/8
                calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
                sigma_outrigger = (M_max_outrigger+M_g)/outrigger.W
                if sigma_outrigger > f:
                    outrigger.t[i] += a
        

        #Kontrolle beta_ist/alpha_ist = beta_outrigger/alpha_outrigger
        # inkl. Erhöhen t mit Berücksichtigung Steifigkeitsverhältnis
        Test = controllStiffnesRatios(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)
        #Test = True

        if Test == True:
            Tragfähigkeitsnachweis = True


       
    #Nachweis belt
    #reicht 10 cm dicke Wand um gesamten grundriss aus?



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
    calculations.calcWeight(buildingProp,materialProp,outrigger)



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
 