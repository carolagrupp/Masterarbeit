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
    b_raster = buildingProp.b_raster 

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
        if buildingProp.j == 0:
            buildingProp.N_max_delta = 0 #Belastung aus Outriggern darüber
            buildingProp.N_kombi_delta = 0

        if buildingProp.outriggerAbschnitt == True and element.typ == "Stütze" and buildingProp.ErsteBerechnungSigma == True:  #Überprüfung ob in Outrigger Abschnitt
            j = buildingProp.j
            #deltaN = N aus Outriggerwirkung + N_EG_Outrigger + N_EG_belt
            deltaN = gamma_w*buildingProp.sigmaStütze_outrigger[j]*element.A\
                +3*buildingProp.t_outrigger_now[buildingProp.posOut[j]]/100*h_geschoss*gamma*gamma_g*buildingProp.b_raster/8*2/5\
                +h_geschoss*b_raster*buildingProp.t_belt_now[buildingProp.posOut[j]]/100*gamma*gamma_g
            buildingProp.N_max_delta += Psi_w*deltaN
            buildingProp.N_kombi_delta += deltaN
            buildingProp.j += 1

        if element.typ == "Kern":
            moment_core = buildingProp.moment[buildingProp.i_aktuell*n_abschnitt-1]
            M_max = moment_core     #nicht: *buildingProp.W_aktuell, da keine Spannung
                
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
        buildingProp.posOut = [n-int(0.61*n)]

    elif n_outrigger == 2:
        buildingProp.posOut = [n-int(0.73*n), n-int(0.41*n) ]

    elif n_outrigger == 3:
        buildingProp.posOut = [n-int(0.80*n), n-int(0.55*n), n-int(0.33*n) ]

    elif n_outrigger == 4:
        buildingProp.posOut = [n-int(0.84*n), n-int(0.64*n), n-int(0.47*n), n-int(0.29*n) ]
    
    else: 
        print(str('Die gewählte Outrigger Anzahl ist nicht verfügbar. Bitte Zahlen zwischen 0 und 4 angeben!'))
        quit()

    #Umwandlung in Abschnitte
    buildingProp.posOut_abschnitt = []
    for i in range(0, len(buildingProp.posOut)):
        buildingProp.posOut_abschnitt.append(buildingProp.posOut[i]//buildingProp.n_abschnitt)
        #buildingProp.posOut_abschnitt.append(buildingProp.x-buildingProp.posOut[i]//buildingProp.n_abschnitt) #Liste mit Liste der Anschnitte in denen Outrigger sind


def outriggerStiffness(I_outriggerExakt, buildingProp):
    #Berechnung der effektiven Biegesteifigkeit
    I_eff = I_outriggerExakt*(1+buildingProp.b_raster/buildingProp.b_raster)**3

    return I_eff

def outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp):
    #Kopplung Outrigger und Stützen
    b_total = buildingProp.b_total
    I_effOutrigger = []
    K = []
    buildingProp.H_outrigger = []
    #Summen von A und I jeweils von einer Hälfte (2 Outrigger, 2 Eckstützen, 3 Randstützen)
    A_randStütze = []
    A_eckStütze = []
    t_randStützeOutrigger = []
    t_eckStützeOutrigger = []
    I_outriggerExakt = []

    #Aktualisierung Querschnitsswerte
    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        t_randStützeOutrigger.append(randStütze.t[n])
        t_eckStützeOutrigger.append(eckStütze.t[n])
        t_outrigger = outrigger.t[i]
        calcProfileProp(outrigger, buildingProp, materialProp, t_outrigger)
        I_aktuell = outrigger.I
        I_outriggerExakt.append(I_aktuell)
    
    #Berechnung Steifigkeiten
    for i in range(0,len(buildingProp.posOut)):
        I_effOutrigger_i = outriggerStiffness(I_outriggerExakt[i], buildingProp)
        I_effOutrigger.append(I_effOutrigger_i)
        buildingProp.I_effOutrigger = I_effOutrigger
        calcProfileProp(randStütze, buildingProp, materialProp, t_randStützeOutrigger[i])
        A_randStütze.append(randStütze.A_min)
        calcProfileProp(eckStütze, buildingProp, materialProp, t_eckStützeOutrigger[i])
        A_eckStütze.append(eckStütze.A_min)
    
    for i in buildingProp.posOut:
        H = (buildingProp.n-i)*buildingProp.h_geschoss    #Höhe des Outriggers
        buildingProp.H_outrigger.append(H)

    for i in range(0, len(buildingProp.posOut)):
        I_outriggerGesamt = 2*I_effOutrigger[i]    #4 Outrigger in Windrichtung
        A_stützen = 3*A_randStütze[i] + 2*A_eckStütze[i]
        K_i = 12*materialProp.E*I_outriggerGesamt*A_stützen*b_total**2/(A_stützen*b_total**3+24*buildingProp.H_outrigger[i]*I_outriggerGesamt)
        K.append(K_i)

    return K

def calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n):
    #Outrigger
    b_total = buildingProp.b_total
    h = buildingProp.h_total
    t=kern.t[n]
    calculations.calcProfileProp(kern, buildingProp, materialProp, t)
    I_core = kern.I
    t = outrigger.t[i]
    calculations.calcProfileProp(outrigger, buildingProp, materialProp, t)
    I_out = outrigger.I
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
    Verhältnis1 = [None]*buildingProp.x
    Verhältnis2 = [None]*buildingProp.x
    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger
    fehler = 0.4
    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        #Berechnung Steifigkeitsverhältnis
        beta_ist, alpha_ist = calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n)
        #Runden für Optimierung (4,4 etc ist gut genug)
        #beta_ist_round = round(beta_ist)
        #alpha_ist_round = round(alpha_ist)
        delta_beta = beta_outrigger-beta_ist
        delta_alpha = alpha_outrigger-alpha_ist

        if abs(delta_beta) <= fehler:
            Verhältnis1[n] = True
        elif delta_beta > fehler:    #Erhöhung Kern
            kern.t[n] += a 
            Verhältnis1[n] = False
        elif delta_beta < -fehler:    #Erhöhung Outrigger
            outrigger.t[i] += a
            Verhältnis1[n] = False

            
        if abs(delta_alpha) <= fehler:
            Verhältnis2[n] = True   
        elif delta_alpha > fehler:
            kern.t[n] += a 
            Verhältnis2[n] = False
        elif delta_alpha < -fehler:
            eckStütze.t[n] += a 
            randStütze.t[n] += a
            Verhältnis2[n] = False

        

    Test = True

    for i in buildingProp.posOut_abschnitt:
        if Verhältnis1[i] == False:
            Test = False
        if Verhältnis2[i] == False:
            Test = False

    return Test

def checkChange(buildingProp, t_check_kern, t_check_eckstütze, t_check_randstütze, t_check_outrigger, kern, eckStütze, randStütze, outrigger):
    length = len(kern.t)
    Check_kern = [False]*length
    Check_eckstütze = [False]*length
    Check_randstütze = [False]*length
    Check_outrigger = []
    fehler = 5 #in cm, ist das sinnvoll?
    #Vergleich t_alt zu t_neu
    for i in range(0, len(kern.t)):
        delta_kern = kern.t[i] - t_check_kern[i]
        if delta_kern < fehler:
            Check_kern[i] = True
        delta_eckstütze = eckStütze.t[i] - t_check_eckstütze[i]
        if delta_eckstütze < fehler:
            Check_eckstütze[i] = True
        delta_randstütze = randStütze.t[i] - t_check_randstütze[i]
        if delta_randstütze < fehler:
            Check_randstütze[i] = True

    for i in buildingProp.posOut:
        delta_outrigger = outrigger.t[i] - t_check_outrigger[i]
        if delta_outrigger < fehler:
            Check_outrigger.append(True)
        else:
            Check_outrigger.append(False)

    Test = True
    for i in range(0, length):
        if Check_kern[i] == False:
            Test = False
        if Check_eckstütze[i] == False:
            Test = False
        if Check_randstütze[i] == False:
            Test = False
    for i in range(0, len(Check_outrigger)):
        if Check_outrigger[i] == False:
            Test = False
    
    return Test

def checkChangeGebäudenachweis(t_check_kern, kern):
    length = len(kern.t)
    Check_kern = [False]*length
    fehler = 5 #in cm, ist das sinnvoll?
    #Vergleich t_alt zu t_neu
    for i in range(0, len(kern.t)):
        delta_kern = kern.t[i] - t_check_kern[i]
        if delta_kern < fehler:
            Check_kern[i] = True
    

    Test = True
    for i in range(0, length):
        if Check_kern[i] == False:
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

    buildingProp.I=I        #Länge Anzahl Abschnitte
    buildingProp.GA=GA
   

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
    belt = building.elements(0,0,'Belt Truss','Rechteckiges Vollprofil')

    str_=str_outrigger
    buildingProp.Iteration = False

    # Tragfähigkeitsnachweise Vertikallasten (Für Elemente ohne Aussteifung):
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kernOhnePFH,buildingProp,loads,materialProp,str_)
        
    calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kern,buildingProp,loads,materialProp,str_)
    
    #Startwert für Steifigkeiten Outriggern
    #Steifigkeitsverhältnis Kern zu Outrigger
    #buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)   #I des Kerns
    #beta_outrigger = buildingProp.beta_outrigger
    #alpha_outrigger = buildingProp.alpha_outrigger
    a = materialProp.delta_t
    f = materialProp.f
    #b = [0] * buildingProp.n
    outrigger.t = [0] * buildingProp.n

    for i in buildingProp.posOut:
        outrigger.t[i] = materialProp.t_min

        #Bestimmen A_outrigger (4 Outrigger je Richtung)
        #I_outrigger_i = materialProp.E/materialProp.E*buildingProp.I[i//buildingProp.n_abschnitt]*b_total/(beta_outrigger*h*4)
        #b[i]=I_outrigger_i*12/(h_geschoss)**3*100  #in cm
        #outrigger.t[i] = round(b[i])
    #for i in buildingProp.posOut_abschnitt:
        #Bestimmen A-Stützen (10 Stützen je Richtung)
        #A_stütze = 2*materialProp.E/materialProp.E*buildingProp.I[i]/(alpha_outrigger*b_total**2*10)
        #t_stütze = numpy.sqrt(A_stütze)*100         #in cm
        #while t_stütze > randStütze.t[i]: #Vergrößern von t wenn nötig
            #randStütze.t[i] += a     #Ersetzen von t für ganzen Abschnitt
        #while t_stütze > eckStütze.t[i]: #Vergrößern von t wenn nötig
            #eckStütze.t[i] += a     #Ersetzen von t für ganzen Abschnitt
    
    buildingProp.Iteration_counter = 0
    Tragfähigkeitsnachweis = False
     

    while Tragfähigkeitsnachweis == False:
        buildingProp.Iteration = True
        buildingProp.Iteration_counter += 1
        #Tragfähigkeitsnachweise Aussteifung:
    
        #Federsteifigkeiten der Outriggergeschossen:
        buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)

        #Steifigkeit Kern
        buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)   #nur core, passt das?

        #fea: calcBendingMoment
        buildingProp.mue=[]
        for i in range (0,len(kern.t)): 
            buildingProp.mue.append(0)

        feModel = fea.feModel(buildingProp, loads, materialProp)
        #Platzhalter
        moment = []
        #for i in range(0,buildingProp.n):
            #moment.append(460000*(i/(buildingProp.n-1)))
        
        moment = fea.feModel.calcBendingMoment(feModel)
        buildingProp.moment = moment    #Anpassen sigma in calcelementLoad
        print(moment)
        print(len(moment))

        #Auflagerkräfte, Länge müsste = 3+n_outrigger
        reactions = fea.feModel.calcreactions(feModel)
        print(reactions)

        print(reactions[2])
        #print(moment[buildingProp.posOut[1]])
    
        #Auslesen Auflagerkraft an Feder
        buildingProp.moment_spring = []
        buildingProp.sigmaStütze_outrigger = []
        I_outrigger  = [None]*len(buildingProp.posOut)
        for n, i in enumerate(buildingProp.posOut):
            calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
            I_outrigger[n] = outriggerStiffness(outrigger.I, buildingProp)

        for i in range(0, len(buildingProp.posOut)):
            #reactions = [H_Fußpunkt, V_Fußpunkt, M_Fußpunkt, M_Feder_1 (unterste), ...] ?
            moment_spring = reactions[i+3]
            buildingProp.moment_spring.append(moment_spring)
            #verdrehung = moment_spring/buildingProp.K[i]
            #Kraft/Fläche (kN/m^2)
            #spannung_stütze = (3*verdrehung-moment_spring*buildingProp.b_total/(2*materialProp.E*4*I_outrigger[i]))\
                #*materialProp.E*buildingProp.b_total/(6*buildingProp.H_outrigger[i])
            Kraft_stütze = moment_spring*2/buildingProp.b_total
            n = buildingProp.posOut_abschnitt[i]
            spannung_stütze = Kraft_stütze/(6*(randStütze.t[n]/100)**2+4*(eckStütze.t[n]/100**2))
            #vielleicht spannung erst in calcelementload berechnen?

            buildingProp.sigmaStütze_outrigger.append(spannung_stütze)

        #Beanspruchung Outrigger
        #Moment aus Feder in fea
        for n,i in enumerate(buildingProp.posOut):
            M_max_outrigger = loads.gamma_w*reactions[n+3]
            sigma_outrigger = 2*f
            while sigma_outrigger > f:
                M_g = outrigger.t[i]/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2/8
                calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
                sigma_outrigger = (M_max_outrigger/4+M_g)/outrigger.W   #da 4 Outrigger je Richtung
                if sigma_outrigger > f:
                    outrigger.t[i] += a
        
        #print(outrigger.t)
        
        #Nachweis belt
        #beansprucht durch Moment aus Auflagerkräften Stützen und Outrigger
        #Mmax in 1./2. Feld
        belt.t = [0] * buildingProp.n
        A_Out = []
    
        #in reactions: Fußpunkt, unterste Feder, ...
        for n,i in enumerate(buildingProp.posOut):
            verdrehung = reactions[2+buildingProp.n_outrigger-n]//buildingProp.K[n]
            A_Out.append(3/8*outrigger.t[i]/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster+\
                3*loads.gamma_w*materialProp.E*buildingProp.I_effOutrigger[n]*verdrehung/b_raster**2)

        for n,i in enumerate(buildingProp.posOut):
            t = materialProp.t_min
            M_belt_Out1 = 1/5*b_raster*A_Out[n] #M_max Feld 1
            M_belt_Out2 = 3/10*b_raster*A_Out[n] #M_max Feld 2
            sigma = 2*f
            while sigma > f:
                M_belt_EG1 = 0.077*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2
                M_belt_EG2 = 0.036*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2
                M_max_belt = max(M_belt_EG1+M_belt_Out1, M_belt_EG2+M_belt_Out2)
                calcProfileProp(belt, buildingProp, materialProp, t)
                sigma = M_max_belt/belt.W
                t += a
            belt.t[i] = t
        
        #Übergabe QS-Werte für calcElementLoads
        buildingProp.t_belt_now = belt.t
        buildingProp.t_outrigger_now = outrigger.t

        #Beanspruchung Stützen (Auflagerkraft Feder + Vertikallasten inkl aus Federn darüber)
        calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(eckStütze, buildingProp,loads,materialProp,str_)

        #Beanspruchung Kern (Beanspruchung fea + Vertikallasten)
        calculations.calcElementWidth(kern, buildingProp,loads,materialProp,str_)
        
        #Kontrolle beta_ist = beta_outrigger/alpha_ist = alpha_outrigger
        # inkl. Erhöhen t mit Berücksichtigung Steifigkeitsverhältnis
        #Test = controllStiffnesRatios(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)
        if buildingProp.Iteration_counter == 1:
            Test = False
        else:
            Test = checkChange(buildingProp, t_check_kern, t_check_eckstütze, t_check_randstütze, t_check_outrigger, kern, eckStütze, randStütze, outrigger)    

        #Speicherung t-Werte für Check in nächster Iteration
        t_check_kern = kern.t
        t_check_eckstütze = eckStütze.t
        t_check_randstütze = randStütze.t
        t_check_outrigger = outrigger.t


        if Test == True:
            Tragfähigkeitsnachweis = True
        
        
    print('Iterationsanzahl = ', buildingProp.Iteration_counter)
        
    
    # Gebäudenachweise:
    t0=kern.t[-1]       #Kerndicke ganz unten (letztes Element der Liste)
    Tragfähigkeitsnachweis = False
    Iteration_counter = 0
    while Tragfähigkeitsnachweis == False:
        Iteration_counter += 1
        calculations.buildingDeflection(buildingProp,loads,materialProp,str_,kern)
        calculations.interstoryDrift(buildingProp,loads,materialProp,str_,kern)
        buildingProp.t_kern=kern.t

        if Iteration_counter == 1:
            Test = False
        else:
            Test = checkChangeGebäudenachweis(t_check_kern, kern)    

        #Speicherung t-Werte für Check in nächster Iteration
        t_check_kern = kern.t

        if Test == True:
            Tragfähigkeitsnachweis = True


    #Vergleich Tragfähigkeit und Gebäudenachweise:
    if kern.t[-1] > t0:     
        buildingProp.NW_maßgebend='Verformung'
        print(kern.t[-1],t0)
    
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'

    # Massen:
    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,randStützeOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,randStütze)
    calculations.calcWeight(buildingProp,materialProp,eckStützeOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,eckStütze)
    calculations.calcWeight(buildingProp,materialProp,kernOhnePFH)
    calculations.calcWeight(buildingProp,materialProp,kern)
    calculations.calcWeight(buildingProp,materialProp,outrigger)
    calculations.calcWeight(buildingProp,materialProp,belt)


    buildingProp.G_innenStützen=innenStütze.G
    buildingProp.G_randStützenOhnePFH=[element*12 for element in randStützeOhnePFH.G]
    buildingProp.G_randStützen=[element*12 for element in randStütze.G]
    buildingProp.G_eckStützenOhnePFH=[element*4 for element in eckStützeOhnePFH.G]
    buildingProp.G_eckStützen=[element*4 for element in eckStütze.G]
    buildingProp.G_kernOhnePFH=kernOhnePFH.G
    
    #Outrigger in Abschnitte umwandeln
    buildingProp.G_outrigger = []
    buildingProp.G_belt = []
    x = buildingProp.x
    n = buildingProp.n
    for i in range(1, x+1):
        buildingProp.G_outrigger.append(4*outrigger.G[i*buildingProp.n_abschnitt-1])
        buildingProp.G_belt.append(belt.G[i*buildingProp.n_abschnitt-1])
    if x*buildingProp.n_abschnitt < n:
        buildingProp.G_outrigger.append(4*outrigger.G[-1])
        buildingProp.G_belt.append(belt.G[-1])
    
    buildingProp.G_aussteifung = [0]*len(kern.G)
    for i in range(0, len(kern.G)):
        buildingProp.G_aussteifung[i] = kern.G[i]+buildingProp.G_outrigger[i]+buildingProp.G_belt[i] 
    
    buildingProp.G_außenStützen=[]
    buildingProp.G_total=[]
    buildingProp.G_totalOhnePFH=[]

    for i in range (0,len(innenStütze.G)):
        G_außenStützen=randStütze.G[i]*12+eckStütze.G[i]*4
        buildingProp.G_außenStützen.append(G_außenStützen)
        buildingProp.G_total.append(G_außenStützen+innenStütze.G[i]+kern.G[i]+buildingProp.G_outrigger[i]+buildingProp.G_belt[i])
        buildingProp.G_totalOhnePFH.append(randStützeOhnePFH.G[i]*12+eckStützeOhnePFH.G[i]*4+innenStütze.G[i]+kernOhnePFH.G[i])
    
    # Eigenfrequenz
    calculations.calcEigenfrequency(buildingProp,loads,materialProp)
    
    buildingProp.t_innenStützen=innenStütze.t
    buildingProp.t_randStützen=randStütze.t
    buildingProp.t_eckStützen=eckStütze.t
    buildingProp.t_kern=kern.t
    buildingProp.t_outrigger = []
    buildingProp.t_belt = []
    for i in buildingProp.posOut:
        buildingProp.t_outrigger.append(outrigger.t[i])
        buildingProp.t_belt.append(belt.t[i])
    buildingProp.t_diagonale=[]
    buildingProp.t_querstrebe=[]
    buildingProp.t_riegel=[]
 