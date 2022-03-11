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
from Masterarbeit.pfh.building import buildingProp, materialProp
from Masterarbeit.pfh.calculations import calcProfileProp
from pfh import building
from pfh import calculations
from pfh import str_outrigger
from pfh import fea
import numpy
# ------------------------------------------------------------------------------

def calcElementLoads(buildingProp,loads,materialProp,element,s,alpha,Ng_darüberliegend,t):                      #in calculations.calcElementWidth
    """Berechnet die Lasten des übergebenen Elements in dem aktuellen Abschnitt
        Lasten aus Geschossen darüber (EG, q) und Lasten aus der Aussteifung
        :param s: Anzahl Geschosse über aktuellem Abschnitt
        :type coords: int
        :param alpha: Abminderungsbeiwert Einzugsfläche (EC1-1-1/NA, 6.3.1.2(10))
        :type alpha: float
        :param Ng_darüberliegend: Normalkraft aus Element darüber [kN]
        :type Ng_darüberliegend: float
        :param t: Dimension des aktuellen Elements [cm]
        :type t: int
        """

    # Übergaben Geometrie
    h_geschoss=buildingProp.h_geschoss
    n_abschnitt=buildingProp.n_abschnitt
    n=buildingProp.n
    x=buildingProp.x
    b_raster = buildingProp.b_raster 

    # Lasten
    Gd=loads.gd
    Qd=loads.qd
    gd_fassade=loads.gd_fassade
    Psi_w=loads.Psi_w
    Psi_q=loads.Psi_q
    
    # Sicherheitsbeiwert und Wichte
    gamma_g=loads.gamma_g
    gamma=materialProp.gamma
    gamma_w=loads.gamma_w

    # Einzugslänge und -fläche des aktuellen Elements
    A_einzug=element.A_einzug
    l_fassade=element.l_fassade

    # Eigenlasten des Elements in diesem Abschnitt
    if  s == n and x * n_abschnitt < n:                                                                                 # wenn n/n_abschnitt nicht aufgeht, unten noch mehrere Geschosse mit ni<n_abschnitt
        Ng = element.A * gamma * (n - x * n_abschnitt) * h_geschoss * gamma_g
    else:
        Ng = element.A * gamma * n_abschnitt * h_geschoss * gamma_g

    N_max   = Ng_darüberliegend + Ng + (A_einzug*(Gd+Qd*alpha)       + l_fassade*gd_fassade)*s
    N_kombi = Ng_darüberliegend + Ng + (A_einzug*(Gd+Qd*alpha*Psi_q) + l_fassade*gd_fassade)*s

    if buildingProp.Iteration == False:
        M_max = 0
    
    elif buildingProp.Iteration == True:
        if buildingProp.j == 0:                                                                                     # vor erstem Outrigger
            buildingProp.N_max_delta = 0                                                                            # Belastung aus Outriggern darüber
            buildingProp.N_kombi_delta = 0

        if buildingProp.outriggerAbschnitt == True and element.typ == "Stütze" and buildingProp.ErsteBerechnungSigma == True:  # Überprüfung ob in Outrigger Abschnitt
            j = buildingProp.j
            deltaNOut = 0

            for outrigger_i in range(0,buildingProp.anzahlOutInAbschnitt):
                deltaNOut += gamma_w*buildingProp.kraftStütze_outrigger[ j + outrigger_i]

            deltaN = deltaNOut\
                +3*buildingProp.t_outrigger_now[buildingProp.posOut[j]]/100*h_geschoss*gamma*gamma_g*buildingProp.b_raster/8*2/5\
                +h_geschoss*b_raster*buildingProp.t_belt_now[buildingProp.posOut[j]]/100*gamma*gamma_g                          # deltaN = N aus Outriggerwirkung + N_EG_Outrigger + N_EG_belt
            
            buildingProp.N_max_delta += Psi_w*deltaN
            buildingProp.N_kombi_delta += deltaN

            buildingProp.j += buildingProp.anzahlOutInAbschnitt                                                                 # Erhöhen der berücksichtigten Outrigger

        if element.typ == "Kern":                                                                                               # Zusatzbelastung Kern aus Aussteifung
            bm = buildingProp.moment
            i_aktuell = buildingProp.i_aktuell                                                                                  # Nummer des Abschnitts
            bm_max_abschnitt = 0
            if x * n_abschnitt < n:
                n_max = n - x*n_abschnitt
                for delta_i in range(0,n_max):
                    bm_i = bm[(i_aktuell-1)*n_abschnitt-delta_i+1]
                    if bm_i > bm_max_abschnitt:
                        bm_max_abschnitt = bm_i
            else:
                for delta_i in range(0,n_abschnitt):
                    bm_i = bm[i_aktuell*n_abschnitt-delta_i]
                    if bm_i > bm_max_abschnitt:
                        bm_max_abschnitt = bm_i
            M_max = gamma_w * bm_max_abschnitt
                
        if element.typ == 'Stütze':                                                                                             # Zusatzbelastung Stütze aus Aussteifung
            M_max = 0
            N_max += buildingProp.N_max_delta
            N_kombi += buildingProp.N_kombi_delta
    
    M_kombi=Psi_w*M_max

    return N_max,N_kombi,M_max,M_kombi


def arrangementOutrigger(buildingProp):
    """Bestimmt die Anordnung der Outrigger über die Höhe nach optimierten Werten der Literatur
        """

    n=buildingProp.n
    n_outrigger = buildingProp.n_outrigger

    
    
    if buildingProp.varPosOut == True:       #Bei Optimierung über Höhe
        buildingProp.posOut = [n-int(buildingProp.xi*n)]

    elif buildingProp.varPosOut_2D == True:       #Bei 2D Optimierung über Höhe
        if buildingProp.xi_1 > buildingProp.xi_2:
            buildingProp.posOut = [n-int(buildingProp.xi_1*n), n-int(buildingProp.xi_2*n)]
        if buildingProp.xi_2 > buildingProp.xi_1:
            buildingProp.posOut = [n-int(buildingProp.xi_2*n), n-int(buildingProp.xi_1*n)]
        #Was wenn sie gleich sind?
    else:
        # Bestimmen der Outrigger Geometrie (nach Park, Lee et all 2016)
        if n_outrigger == 0:
            buildingProp.posOut = []

        elif n_outrigger == 1:
            buildingProp.posOut = [n-int(0.61*n)]#[0] # ganz oben für Untersuchung[n-int(0.61*n)]

        elif n_outrigger == 2:
            buildingProp.posOut = [n-int(0.73*n), n-int(0.41*n) ]

        elif n_outrigger == 3:
            buildingProp.posOut = [n-int(0.80*n), n-int(0.55*n), n-int(0.33*n) ]

        elif n_outrigger == 4:
            buildingProp.posOut = [n-int(0.84*n), n-int(0.64*n), n-int(0.47*n), n-int(0.29*n) ]
        
        else: 
            print(str('Die gewählte Outrigger Anzahl ist nicht verfügbar. Bitte Zahlen zwischen 0 und 4 angeben!'))
            quit()

    # Umwandlung in Abschnitte
    buildingProp.posOut_abschnitt = []
    for posOut in buildingProp.posOut:
        buildingProp.posOut_abschnitt.append(posOut//buildingProp.n_abschnitt)              # Liste mit Liste der Abschnitte in denen Outrigger sind

def changeColumn(buildingProp, materialProp, randStütze, eckStütze):
    # Nicht in Verwendung
    """Berechnet die Dimension der Stützen über das Steifigkeitsverhältnis zum Kern, 
        dabei wird zwischen Eck- und Randstütze unterschieden
        """
    I = buildingProp.I
    E = materialProp.E

    alpha_outrigger = buildingProp.alpha_outrigger
    b_total = buildingProp.b_total
    a = materialProp.delta_t
    #A_stütze_erf_ges = []
    #A_vorh_ges = []

    for abschnittOut in buildingProp.posOut_abschnitt:
        A_stütze_erf = 2*E/E*I[abschnittOut]/(alpha_outrigger*b_total**2)                # m**2
        A_stütze_erf_ges = 5*A_stütze_erf
        A_eck_vorh = (eckStütze.t[abschnittOut]/100)**2
        A_rand_vorh = (randStütze.t[abschnittOut]/100)**2
        A_stütze_vorh_ges = 2*A_eck_vorh + 3*A_rand_vorh

        eta = A_stütze_erf_ges/A_stütze_vorh_ges

        A_eck_neu = eta*A_eck_vorh
        t_eck_neu = numpy.sqrt(A_eck_neu)*100
        A_rand_neu = eta*A_rand_vorh
        t_rand_neu = numpy.sqrt(A_rand_neu)*100
        buildingProp.t_stütze.append(t_rand_neu)        #auch noch für Eck bzw gesamt
        while t_rand_neu > randStütze.t[abschnittOut]:                                           # Vergrößern von t wenn nötig
            randStütze.t[abschnittOut] += a                                                    # Ersetzen von t für ganzen Abschnitt
        while t_eck_neu > eckStütze.t[abschnittOut]:                                            # Vergrößern von t wenn nötig
            eckStütze.t[abschnittOut] += a                                                     # Ersetzen von t für ganzen Abschnitt

        A_neu_vorh = 2*(eckStütze.t[abschnittOut]/100)**2 + 3*(randStütze.t[abschnittOut]/100)**2

def ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze):
    """Berechnet die Dimension des Outriggers über das Steifigkeitsverhältnis zum Kern 
    und überprüft das Steifigkeitsverhältnis zwischen Stützen und Kern, bei Bedarf wird der Stützen QS erhöht
        """

    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger
    b_total = buildingProp.b_total
    b_raster = buildingProp.b_raster
    h_geschoss = buildingProp.h_geschoss
    h = buildingProp.h_total
    a = materialProp.delta_t

    # Berechnung Steifigkeit des Kerns
    buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)

    I = buildingProp.I
    E = materialProp.E

    for posOut in buildingProp.posOut:
        I_outEff_soll_i = E/E*I[posOut//buildingProp.n_abschnitt]*b_total/(beta_outrigger*h)          # m**4
        I_outrigger_soll_i = I_outEff_soll_i/(1+b_raster/b_raster)**3
        b_erf=I_outrigger_soll_i*12/(h_geschoss)**3*100                                  # cm
        b_vorh = outrigger.t[posOut]
        while b_vorh < b_erf:
            b_vorh += a
        outrigger.t[posOut] = b_vorh


    buildingProp.t_stütze = []

    # Erhöhen der Stützenquerschnitte, wenn kleiner als erforderliches Steifigkeitsverhältnisses
    for abschnittOut in buildingProp.posOut_abschnitt:
        A_stütze = 2*E/E*I[abschnittOut]/(alpha_outrigger*b_total**2)                # m**2
        t_stütze = numpy.sqrt(A_stütze)*100                                         # cm
        buildingProp.t_stütze.append(t_stütze)

        while t_stütze > randStütze.t[abschnittOut]:                                           # Vergrößern von t wenn nötig
            randStütze.t[abschnittOut] += a                                                    # Ersetzen von t für ganzen Abschnitt
        while t_stütze > eckStütze.t[abschnittOut]:                                            # Vergrößern von t wenn nötig
            eckStütze.t[abschnittOut] += a                                                     # Ersetzen von t für ganzen Abschnitt
    
        # Stützenabschnitte unter Outrigger sollten auch mindestens so groß wie Outriggerstütze sein
        for abschnitte in range(abschnittOut+1, len(randStütze.t)):
            while randStütze.t[abschnittOut] > randStütze.t[abschnitte]:                                           # Vergrößern von t wenn nötig
                randStütze.t[abschnitte] += a                                                    # Ersetzen von t für ganzen Abschnitt
            while eckStütze.t[abschnittOut] > eckStütze.t[abschnitte]:                                            # Vergrößern von t wenn nötig
                eckStütze.t[abschnitte] += a                                                     # Ersetzen von t für ganzen Abschnitt
    #Alternative Unterscheidung Eck- und Randstütze
    #changeColumn(buildingProp, materialProp, randStütze, eckStütze)


def outriggerStiffness(I_outrigger, buildingProp):
    """Berechnet die effektive Biegesteifigkeit des aktuellen Outriggers
        :param I_outrigger: Flächenträgheitsmoment des Outriggerquerschnitts [m**4]
        :type I_outrigger: float
        """

    b_raster = buildingProp.b_raster

    I_outEff = I_outrigger*(1+b_raster/b_raster)**3

    return I_outEff

def outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp):
    """Berechnet die Federsteifigkeiten der Outrigger
        """
    # Übergabe der Gebäudebreite
    b_total = buildingProp.b_total
    n_abschnitt = buildingProp.n_abschnitt
    n = buildingProp.n
    x=buildingProp.x
    h_geschoss = buildingProp.h_geschoss

    I_outEff = []
    K = []
    buildingProp.H_outrigger = []

    # Summen von A und I jeweils von einer Hälfte (2 Outrigger, 2 Eckstützen, 3 Randstützen)
    A_randStütze = []                                                                                   # m**2
    A_eckStütze = []
    t_randStützeOutrigger = []                                                                          # cm
    t_eckStützeOutrigger = []
    I_outrigger = []                                                                               # m**4
    A_stützen = []
    K_stütze = []

    # Berechnung der Federkonstante Stützen in Reihe
    l_abschnitte = [n_abschnitt * h_geschoss] * x                                      # Liste mit Längen der Abschnitte
    if x*n_abschnitt < n:
        l_abschnitte.append(h_geschoss * (n - n_abschnitt * x))

    for abschnitt in range(len(randStütze.t)):                                                             # Liste der Summe der Stützen je Abschnitt (von einer Seite: 3*Rand+2*Eck)
        A_i_rand = (randStütze.t[abschnitt]/100)**2
        A_i_eck = (eckStütze.t[abschnitt]/100)**2
        A_i_ges = 3*A_i_rand + 2*A_i_eck
        A_stützen.append(A_i_ges)

    for posOut in buildingProp.posOut:
        H = (buildingProp.n-posOut)*buildingProp.h_geschoss                                                 # Höhe des Outriggers
        buildingProp.H_outrigger.append(H)

    j = 0
    
    for abschnittOut in buildingProp.posOut_abschnitt:
        K_stütze_inv = 0
        l_berechnet = 0
        anzahl_abschnitte = len(randStütze.t)-abschnittOut                      # Wie viele Querschnittsabschnitte unter dem aktuellen Outrigger
        for abschnitte in range(anzahl_abschnitte-1):                           # Berechnung für jeden Querschnitt, beginnend von unten
            abschnitt = len(l_abschnitte) -1 -abschnitte                        # Umrechnung auf Abschnitte von unten
            K_stütze_inv += l_abschnitte[abschnitt] /(materialProp.E*A_stützen[abschnitt])
            l_berechnet += l_abschnitte[abschnitt]
        abschnitt = len(l_abschnitte) -(anzahl_abschnitte-1)
        if abschnittOut < abschnitt:
            l_übrig = buildingProp.H_outrigger[j] - l_berechnet
            abschnitt -= 1
            K_stütze_inv += l_übrig /(materialProp.E*A_stützen[abschnitt])
        K_stütze.append(1/K_stütze_inv)
        j += 1


    # Alles folgende könnte in eine Schleife

    # Aktualisierung Querschnittswerte
    for posOut in buildingProp.posOut:
        abschnittOut = posOut//buildingProp.n_abschnitt
        #t_randStützeOutrigger.append(randStütze.t[abschnittOut])
        #t_eckStützeOutrigger.append(eckStütze.t[abschnittOut])
        t_outrigger = outrigger.t[posOut]
        calcProfileProp(outrigger, buildingProp, materialProp, t_outrigger)
        I_outrigger.append(outrigger.I)
    
    # Berechnung einzelne Querschnitssflächen und Steifigkeiten
    for i in range(0,len(buildingProp.posOut)):
        I_outEff_i = outriggerStiffness(I_outrigger[i], buildingProp)
        I_outEff.append(I_outEff_i)
        buildingProp.I_effOutrigger = I_outEff
        #calcProfileProp(randStütze, buildingProp, materialProp, t_randStützeOutrigger[i])
        #A_randStütze.append(randStütze.A_min)
        #calcProfileProp(eckStütze, buildingProp, materialProp, t_eckStützeOutrigger[i])
        #A_eckStütze.append(eckStütze.A_min)
    
    

    for i in range(0, len(buildingProp.posOut)):
        I_outriggerGesamt = 2*I_outEff[i]                                                         # 4 Outrigger in Windrichtung
        #A_stützen = 3*A_randStütze[i] + 2*A_eckStütze[i]
        #K_i = 12*materialProp.E*I_outriggerGesamt*A_stützen*b_total**2/(A_stützen*b_total**3+24*buildingProp.H_outrigger[i]*I_outriggerGesamt)
        K_i = 12*materialProp.E*I_outriggerGesamt*K_stütze[i]*b_total**2/(K_stütze[i]*b_total**3+24*materialProp.E*I_outriggerGesamt)
        K.append(K_i)

    return K                                                                                            # kNm


def checkChange(buildingProp, materialProp, t_check_kern, t_check_eckstütze, t_check_randstütze, t_check_outrigger, t_check_belt, kern, eckStütze, randStütze, outrigger, belt):
    """Überprüft ob sich Dimensionen der Bauteile seit der letzten Iteration verändert haben
        Wenn nicht, Abbruch der Schleife
        """
    length = len(kern.t)
    Check_kern = [False]*length
    Check_eckstütze = [False]*length
    Check_randstütze = [False]*length
    Check_outrigger = []
    Check_belt = []
    fehler = 0.9*materialProp.delta_t                                                              # cm, ist das sinnvoll?

    # Vergleich t_alt zu t_neu
    for i in range(0, len(kern.t)):
        delta_kern = kern.t[i] - t_check_kern[i]
        if abs(delta_kern) < fehler:
            Check_kern[i] = True
        delta_eckstütze = eckStütze.t[i] - t_check_eckstütze[i]
        if abs(delta_eckstütze) < fehler:
            Check_eckstütze[i] = True
        delta_randstütze = randStütze.t[i] - t_check_randstütze[i]
        if abs(delta_randstütze) < fehler:
            Check_randstütze[i] = True

    for posOut in buildingProp.posOut:
        delta_outrigger = outrigger.t[posOut] - t_check_outrigger[posOut]
        if abs(delta_outrigger) < fehler:
            Check_outrigger.append(True)
        else:
            Check_outrigger.append(False)

        delta_belt = belt.t[posOut] - t_check_belt[posOut]
        if abs(delta_belt) < fehler:
            Check_belt.append(True)
        else:
            Check_belt.append(False)



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
        if Check_belt[i] == False:
            Test = False

    
    return Test

def checkChangeGebäudenachweis(materialProp, t_check_kern, kern):
    """Überprüft ob sich Dimensionen des Kerns seit der letzten Iteration verändert haben
        Wenn nicht, Abbruch der Schleife
        """
    length = len(kern.t)
    Check_kern = [False]*length
    fehler = 0.9*materialProp.delta_t                                                          # in cm, ist das sinnvoll?
    
    # Vergleich t_alt zu t_neu
    for i in range(0, len(kern.t)):
        delta_kern = kern.t[i] - t_check_kern[i]
        if delta_kern < fehler:
            Check_kern[i] = True
    
    Test = True
    for i in range(0, length):
        if Check_kern[i] == False:
            Test = False
    
    return Test


def buildingStiffness(buildingProp,materialProp,kern,element2,element3,element4):                                           # wird in calculations.buildingDeflection und interstoryDrift aufgerufen
    """Berechnet die Steifigkeit des Kerns (bei Outrigger) je Abschnitt
        """
    I=[]
    GA=[]

    for i in range (0,len(kern.t)):                                                                                         # Länge = Anzahl Abschnitte
        calculations.calcProfileProp(kern, buildingProp, materialProp, kern.t[i])
            
        I.append(kern.I) 
        GA.append(kern.A_min/2*5/6*materialProp.G)                                                                          # Halbe Kernfläche, da nur die Wände in Kraftrichtung ("Stegwände") die Schubkraft aufnehmen. Abmidnerung mit Kappa (s. Engelke S.14)

    buildingProp.I=I
    buildingProp.GA=GA
   
def bendingStiffnessModification(buildingProp, materialProp, kern, randStütze, eckStütze, outrigger, delta_t):
    """Erhöht die Querschnittsdimension der Tragwerkselemente die an der Aussteifung beteiligt sind, da diese maßgebend für die Biegeverformung ist
        """
    # Erhöhung mit Berücksichtigung der vorgegebenen Steifigkeitsverhältnisse
    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger#*10**(-3)
    b_total = buildingProp.b_total
    h = buildingProp.h_total
    h_geschoss = buildingProp.h_geschoss
    minderung_I=materialProp.minderung_I
    buildingProp.t_stütze = []
    b_kern = 2*buildingProp.b_raster

    # Zunächst erhöhen des Outriggers (dann aber große Erhöhung direkt, zu großer Schritt)
    #for stelle in buildingProp.posOut:
        #outrigger.t[stelle] += delta_t
        #calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[stelle])
        #stelle_abschnitt = stelle//buildingProp.n_abschnitt
        # Anpassen Kern
        #I_kern_i = materialProp.E/materialProp.E*beta_outrigger*4*outrigger.I*h/b_total
        #I_stern = I_kern_i/minderung_I*12*100**4/(8*b_kern)
        #t_soll = 0.30285*(1.2599*((12*b_kern**6+81*I_stern**2)**(1/2)+9*I_stern)**(2/3)-2.8845*b_kern**2)\
                                #/((12*b_kern**6+81*I_stern**2)**(1/2)+9*I_stern)**(1/3)
        #while kern.t[stelle_abschnitt] < t_soll:
            #kern.t[stelle_abschnitt] += delta_t
        #buildingStiffness(buildingProp, materialProp, kern,element2=None, element3=None, element4=None)
        # Anpassen Stützen
        #A_stütze = 2*materialProp.E/materialProp.E*buildingProp.I[stelle_abschnitt]/(alpha_outrigger*b_total**2*10)
        #t_stütze = numpy.sqrt(A_stütze)*100                                         # cm
        #buildingProp.t_stütze.append(t_stütze)
        
        #while t_stütze > randStütze.t[stelle_abschnitt]:                                           # Vergrößern von t wenn nötig
            #randStütze.t[stelle_abschnitt] += delta_t                                                    # Ersetzen von t für ganzen Abschnitt
        #while t_stütze > eckStütze.t[stelle_abschnitt]:                                            # Vergrößern von t wenn nötig
            #eckStütze.t[stelle_abschnitt] += delta_t

    # Aktualisierung K_Feder für fea
    #buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)
    # Müssen Stützen etc in Geschossen darunter wieder neu bemessen werden
    # deswegen GzG mit in Iterationsschleife von GzT?

    #-------------------------------------------------------------------------------------------------------
    ## oder
    # Anpassen des Kerns
    kern.t = [element + delta_t for element in kern.t]

    # Aktualisieren Steifigkeit des Kerns
    buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)

    # Aktualisieren Outrigger und Stützen
    ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

    # Aktualisierung K_Feder für fea
    buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)




def shearStiffnessModification(buildingProp, materialProp, kern, randStütze, eckStütze, outrigger, delta_t):
    """Erhöht die Querschnittsdimension des Kerns, da dieser maßgebend für die Schubverformung ist
        """

    kern.t = [element + delta_t for element in kern.t]

    # Aktualisieren Steifigkeit des Kerns
    buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)

    # Aktualisieren Outrigger und Stützen
    ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

    # Aktualisierung K_Feder für fea
    buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)




def design(buildingProp,loads,materialProp,DataProp):                                                                                # Aufruf in app.py.MainWondow.MainCalculation
    """Führt die Bemessung  und Analyse des Outrigger-Tragwerks durch
        """
    # Anordnung der Outrigger über die Höhe nach Eingabe gui bestimmen
    arrangementOutrigger(buildingProp)

    # Gebäudegeometrie übergeben:
    b_raster = buildingProp.b_raster
    h_geschoss = buildingProp.h_geschoss
    A = b_raster**2

    # Tragwerksteile initiieren
    # mit Aussteifung
    innenStütze=building.elements(A,0,'Stütze','Vollprofil', 'innenStütze')                # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStütze=building.elements(1/2*A,b_raster,'Stütze','Vollprofil', 'randStütze')
    eckStütze=building.elements(1/4*A,b_raster,'Stütze','Vollprofil', 'eckStütze')
    kern=building.elements(8*A,0,'Kern','Kern', 'kern')
    outrigger=building.elements(0,0,'Outrigger','Rechteckiges Vollprofil', 'outrigger')
    belt = building.elements(0,0,'Belt Truss','Rechteckiges Vollprofil', 'belt')

    # ohne Aussteifung
    randStützeOhnePFH=building.elements(1/2*A,b_raster,'Stütze ohne PFH','Vollprofil', 'randStützeOhnePFH')
    eckStützeOhnePFH=building.elements(1/4*A,b_raster,'Stütze ohne PFH','Vollprofil', 'eckStützeOhnePFH')
    kernOhnePFH=building.elements(8*A,0,'Kern ohne PFH','Kern', 'kernOhnePFH')

    # Aktuelles Tragsystem für calculation zuordnen
    str_=str_outrigger

    # Übergabe an Bemessung, dass noch keine Iteration (ohne Outrigger)
    buildingProp.Iteration = False

    # Übergabe Materialeigenschaften
    a = materialProp.delta_t
    t = materialProp.t_min
    f = materialProp.f

    # Tragfähigkeitsnachweise Vertikallasten (Für Elemente ohne Aussteifung):
    if buildingProp.Nachweis_GZT == True:

        # Für Elemente ohne Aussteifung
        calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(randStützeOhnePFH,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(eckStützeOhnePFH,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(kernOhnePFH,buildingProp,loads,materialProp,str_)
        
        # Für Elemente mit Aussteifung als Startwert
        calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)
        #randStütze.t = [materialProp.t_min]*len(innenStütze.t)
        #eckStütze.t = [materialProp.t_min]*len(innenStütze.t)
        calculations.calcElementWidth(kern,buildingProp,loads,materialProp,str_)
        #kern.t = [materialProp.t_min]*len(innenStütze.t)

    else:
        if buildingProp.x*buildingProp.n_abschnitt < buildingProp.n:   #ein weiterer Abschnitt nötig, da es unten noch einen zusätzlichen Abschnitt gibt
            z = buildingProp.x+1
        
    
        else:               #wenn x*n_abschnitt = n
            z = buildingProp.x

        randStütze.t = [t]*z
        eckStütze.t = [t]*z
        kern.t = [t]*z
        innenStütze.t = [t]*z
        belt.t = [0] * buildingProp.n



    

    # Erstellung Liste Outriggerdimension über n = Anzahl Geschosse
    outrigger.t = [0] * buildingProp.n

    # Iterationszähler initialisieren
    buildingProp.Iteration_counter = 0

    # Schleifenbedingung initialisieren
    Tragfähigkeitsnachweis_erbracht = False
    #t_check_eckstütze = []
     

    while Tragfähigkeitsnachweis_erbracht == False:

        # Angabe nun Iteration aktiv
        buildingProp.Iteration = True

        # Erhöhen des Iterationszählers
        buildingProp.Iteration_counter += 1

        
        # Berechnung Steifigkeit des Kerns
        buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)
        
        # Berechnung/ Aktualisierung des Outriggers und der Stützen bezogen auf Steifigkeitsverhältnis zu Kern
        ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

        # Federsteifigkeiten der Outriggergeschosse:
        buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)

        # Initialisierung fea Model für Kalkulation Momentenverlauf und Auflagerkräfte
        buildingProp.mue=[]
        for i in range (0,len(kern.t)): 
            buildingProp.mue.append(0)

        feModel = fea.feModel(buildingProp, loads, materialProp)

        if buildingProp.Nachweis_GZT == True:
        
            # Auslesen Biegemomente über die Höhe
            moment = fea.feModel.calcBendingMoment(feModel)
            
            buildingProp.moment_ges = moment                                                # kNm                 Array mit 2 Momenten je Element bzw Geschoss
            
            buildingProp.moment = [0]                                                                           # Array mit maximalem Moment je Geschoss
            for i in range(0,buildingProp.n):
                bm_max = max(abs(moment[2*i]), abs(moment[2*i+1]))
                buildingProp.moment.append(bm_max)

            # Auslesen Auflagerkräfte (3x Auflager Fußpunkt + Feder(n))
            # reactions = [H_Fußpunkt, V_Fußpunkt, M_Fußpunkt, M_Feder_1 (oberste), ...]
            reactions = fea.feModel.calcreactions(feModel)
            print(reactions)
        
            # Auslesen Auflagerkraft an Feder
            buildingProp.moment_spring = []
            buildingProp.kraftStütze_outrigger = []                                         # Kraft durch Outrigger von oben nach unten
            I_outEff  = [None]*len(buildingProp.posOut)

            for n, posOut in enumerate(buildingProp.posOut):
                calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[posOut])
                I_outEff[n] = outriggerStiffness(outrigger.I, buildingProp)              # m**4

            for i in range(0, len(buildingProp.posOut)):
                moment_spring = reactions[i+3]                                              # kNm
                buildingProp.moment_spring.append(moment_spring)

                n = buildingProp.posOut_abschnitt[i]                                        # Positionen von oben nach unten
                
                A_stützen_ges = 6*(randStütze.t[n]/100)**2+4*(eckStütze.t[n]/100)**2
                H = buildingProp.H_outrigger[i]
                Kraft_stütze = 3*moment_spring/buildingProp.K[i]*materialProp.E*A_stützen_ges*4*I_outEff[i]*2*b_raster\
                    /(3*H*4*I_outEff[i]+A_stützen_ges*(2*b_raster)**3)                       # kN        aus KV, ergibt höhere Last als I_outrigger und l = b_raster
                
                buildingProp.kraftStütze_outrigger.append(Kraft_stütze/10)                   # Kraft je Stütze

            
            #---------------------------------------------------------------------------------------------------
            # Tragfähigkeitsnachweise Aussteifung:
            #---------------------------------------------------------------------------------------------------

            # Nachweis Outrigger

            # Weglassen, da Beanspruchung nicht einfach betrachtbar
            # Wird als nachgewiesen angenommen
            
            
            # Nachweis Belt Truss

            # System: Durchlaufträger von Eck- zu Eckstütze, beansprucht durch Moment aus Auflagerkräften Stützen und Outrigger

            belt.t = [0] * buildingProp.n
            A_Out = []
            buildingProp.sigma_belt =[]
        
            for n,posOut in enumerate(buildingProp.posOut):
                A_Out.append(3/8*outrigger.t[posOut]/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster+\
                    2.5*buildingProp.kraftStütze_outrigger[n]*loads.gamma_w)                                                       # Normalkraft auf 2,5 Stützen, da es zwei Outrigger je Seite gibt

            for n,posOut in enumerate(buildingProp.posOut):
                t = materialProp.t_min
                sigma = 2*f

                M_belt_Out1 = 1/5*b_raster*A_Out[n]                                                             # M_max Feld 1
                M_belt_Out2 = 3/10*b_raster*A_Out[n]                                                            # M_max Feld 2
                
                while sigma > f:
                    M_belt_EG1 = 0.077*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2            # kNm
                    M_belt_EG2 = 0.036*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2            # kNm
                    M_max_belt = max(M_belt_EG1+M_belt_Out1, M_belt_EG2+M_belt_Out2)
                    calcProfileProp(belt, buildingProp, materialProp, t)
                    sigma = M_max_belt/belt.W                                                                   # kN/m**2
                    t += a
                buildingProp.sigma_belt.append(sigma)
                belt.t[posOut] = t

            # Gleichsetzen der Beltdimension mit maximalem t des Abschnitts
            for abschnitt in range(0, len(kern.t)):
                t_max = 0
                for i in range(0,len(buildingProp.posOut)):
                    if buildingProp.posOut_abschnitt[i] == abschnitt and belt.t[buildingProp.posOut[i]]>t_max:
                        t_max = belt.t[buildingProp.posOut[i]]

                for i in range(0,len(buildingProp.posOut)):
                    if buildingProp.posOut_abschnitt[i] == abschnitt:
                        belt.t[buildingProp.posOut[i]] = t_max
        
            buildingProp.t_belt_now = belt.t
            buildingProp.t_outrigger_now = outrigger.t


            # Nachweis Stützen
            # Beanspruchung: Auflagerkraft Feder + Vertikallasten inkl aus Federn darüber
            calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
            calculations.calcElementWidth(eckStütze, buildingProp,loads,materialProp,str_)
                
            # Nachweis Kern 
            # Beanspruchung: Momente aus fea + Vertikallasten
            calculations.calcElementWidth(kern, buildingProp,loads,materialProp,str_)
        
        # Kontrolle Querschnittsänderung für Iterationsabbruch
        if buildingProp.Iteration_counter == 1:
            Test = False
            #Test = True #Für Vergleich ohne Iteration
        elif buildingProp.Nachweis_GZT == False:
            Test = True

        else:
            Test = checkChange(buildingProp, materialProp, t_check_kern, t_check_eckstütze, t_check_randstütze, t_check_outrigger, t_check_belt, kern, eckStütze, randStütze, outrigger, belt)    

        # Speicherung t-Werte für Check in nächster Iteration
        t_check_kern = []
        t_check_eckstütze = []
        t_check_randstütze = []
        for i in range(len(eckStütze.t)):
            t_check_eckstütze.append(eckStütze.t[i])
            t_check_randstütze.append(randStütze.t[i])
            t_check_kern.append(kern.t[i])

        t_check_outrigger = [0] * buildingProp.n
        t_check_belt = [0] * buildingProp.n
        for posOut in buildingProp.posOut:
            t_check_outrigger[posOut] = outrigger.t[posOut]
            t_check_belt[posOut] = belt.t[posOut]


        if Test == True:
            Tragfähigkeitsnachweis_erbracht = True
        
        
    print('Iterationsanzahl = ', buildingProp.Iteration_counter)
        
    
    # Gebäudenachweise:
    t0=kern.t[-1]                                                                   #Kerndicke ganz unten

    #Tragfähigkeitsnachweis = False
    #Iteration_counter = 0

    #while Tragfähigkeitsnachweis == False:
        #Iteration_counter += 1
    calculations.buildingDeflection(buildingProp,loads,materialProp,str_,kern, randStütze, eckStütze, outrigger)
    calculations.interstoryDrift(buildingProp,loads,materialProp,str_,kern, randStütze, eckStütze, outrigger)
    buildingProp.t_kern=kern.t

        #if Iteration_counter == 1:
            #Test = False
        #else:
            #Test = checkChangeGebäudenachweis(materialProp, t_check_kern, kern)    

        #t_check_kern = kern.t                                                       #Speicherung t-Werte für Check in nächster Iteration

        #if Test == True:
            #Tragfähigkeitsnachweis = True

    # Vergleich Tragfähigkeit und Gebäudenachweise:
    if kern.t[-1] > t0:     
        buildingProp.NW_maßgebend='Verformung'
        print(kern.t[-1],t0)
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'

    print('maßgebender NW', buildingProp.NW_maßgebend)

    # Massenkalkulation:                                                            # in t
    if buildingProp.Nachweis_GZT == True:
        calculations.calcWeight(buildingProp,materialProp,randStützeOhnePFH)
        calculations.calcWeight(buildingProp,materialProp,eckStützeOhnePFH)
        calculations.calcWeight(buildingProp,materialProp,kernOhnePFH)

    
    else:
        randStützeOhnePFH.G = [0]*z
        eckStützeOhnePFH.G = [0]*z
        kernOhnePFH.G = [0]*z
        
        for posOut in buildingProp.posOut:
            belt.t[posOut] = outrigger.t[posOut]

    calculations.calcWeight(buildingProp,materialProp,innenStütze)
    calculations.calcWeight(buildingProp,materialProp,randStütze)
    calculations.calcWeight(buildingProp,materialProp,eckStütze)
    calculations.calcWeight(buildingProp,materialProp,kern)
    calculations.calcWeight(buildingProp,materialProp,outrigger)
    calculations.calcWeight(buildingProp,materialProp,belt)


    buildingProp.G_innenStützen=innenStütze.G
    buildingProp.G_randStützenOhnePFH=[element*12 for element in randStützeOhnePFH.G]
    buildingProp.G_randStützen=[element*12 for element in randStütze.G]
    buildingProp.G_eckStützenOhnePFH=[element*4 for element in eckStützeOhnePFH.G]
    buildingProp.G_eckStützen=[element*4 for element in eckStütze.G]
    buildingProp.G_kernOhnePFH=kernOhnePFH.G
    
    # Outrigger in Abschnitte umwandeln
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
    
    buildingProp.t_innenStützen=innenStütze.t                   # cm
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

    # Querkraftverlauf bestimmen
    bm = buildingProp.moment_ges
    bm_v = []
    for i in range(0, buildingProp.n):
        bm_v.append(bm[2*i])
    bm_v.append(bm[2*buildingProp.n-1])
    v = numpy.diff(bm_v, n = 1)
    buildingProp.querkraft = v                                  # kN


    # Berechnen des Steifigkeitsverhältnis zwischen Kern und Stützen am Ende
    buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)
    buildingProp.alpha_ist_rand = []
    buildingProp.alpha_ist_eck = []
    for abschnittOut in buildingProp.posOut_abschnitt:
        A_randstütze = (randStütze.t[abschnittOut]/100)**2
        A_eckstütze = (eckStütze.t[abschnittOut]/100)**2
        alpha_ist_rand = 2*materialProp.E*buildingProp.I[abschnittOut]/(materialProp.E*A_randstütze*buildingProp.b_total**2)
        buildingProp.alpha_ist_rand.append(round(alpha_ist_rand,2))
        alpha_ist_eck = 2*materialProp.E*buildingProp.I[abschnittOut]/(materialProp.E*A_eckstütze*buildingProp.b_total**2)
        buildingProp.alpha_ist_eck.append(round(alpha_ist_eck,2))


    # Überprüfen ob Stütze größer als erf für alpha ist (keine Erhöhung)
    if buildingProp.Nachweis_GZT == True:

        buildingProp.alpha_outriggerFulfilled = True

        for n,abschnittOut in enumerate(buildingProp.posOut_abschnitt):
            if randStütze.t[abschnittOut] - a > buildingProp.t_stütze[n]:
                buildingProp.alpha_outriggerFulfilled = False
            elif eckStütze.t[abschnittOut] - a > buildingProp.t_stütze[n]:
                buildingProp.alpha_outriggerFulfilled = False

        print('Steifigkeitsverhältnis alpha zwischen Kern und Stützen ist ', buildingProp.alpha_outriggerFulfilled)

        print('M_Feder', buildingProp.moment_spring)
        print('N_Feder', buildingProp.kraftStütze_outrigger)
        print('M_Fuss', buildingProp.moment[-1])
 