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
            
            #buildingProp.sigmaStütze_outrigger[j]*element.A
            deltaNOut = 0
            for outrigger_i in range(0,buildingProp.anzahlOutInAbschnitt):
                deltaNOut += gamma_w*buildingProp.kraftStütze_outrigger[ j + outrigger_i]
            deltaN = deltaNOut\
                +3*buildingProp.t_outrigger_now[buildingProp.posOut[j]]/100*h_geschoss*gamma*gamma_g*buildingProp.b_raster/8*2/5\
                +h_geschoss*b_raster*buildingProp.t_belt_now[buildingProp.posOut[j]]/100*gamma*gamma_g                          # deltaN = N aus Outriggerwirkung + N_EG_Outrigger + N_EG_belt
            buildingProp.N_max_delta += Psi_w*deltaN
            buildingProp.N_kombi_delta += deltaN
            buildingProp.j += buildingProp.anzahlOutInAbschnitt

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

    # Bestimmen der Outrigger Geometrie (nach Park, Lee et all 2016)
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

    # Umwandlung in Abschnitte
    buildingProp.posOut_abschnitt = []
    for i in range(0, len(buildingProp.posOut)):
        buildingProp.posOut_abschnitt.append(buildingProp.posOut[i]//buildingProp.n_abschnitt)              # Liste mit Liste der Abschnitte in denen Outrigger sind
 
def ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze):
    """Berechnet die Dimension des Outriggers über das Steifigkeitsverhältnis zum Kern 
    und überprüft das Steifigkeitsverhältnis zwischen Stützen und Kern, bei Bedarf wird der Stützen QS erhöht
        """

    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger#*10**(-3)
    b_total = buildingProp.b_total
    b_raster = buildingProp.b_raster
    h_geschoss = buildingProp.h_geschoss
    h = buildingProp.h_total
    a = materialProp.delta_t
    # Berechnung Steifigkeit des Kerns
    buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)   # nur core, passt das?

    # Dimensionierung des Outriggers über Steifigkeitsverhältnis zu Kern
    #for i in buildingProp.posOut:
        #outrigger.t[i] = materialProp.t_min                                         # Dicke t auf t_min zurücksetzen

    for i in buildingProp.posOut:
        I_outrigger_eff_i = materialProp.E/materialProp.E*buildingProp.I[i//buildingProp.n_abschnitt]*b_total/(beta_outrigger*h)          # m**4
        I_outrigger_i = I_outrigger_eff_i/(1+b_raster/b_raster)**3
        b_erf=I_outrigger_i*12/(h_geschoss)**3*100                                  # cm
        b_vorh = outrigger.t[i]
        while b_vorh < b_erf:
            b_vorh += a
        outrigger.t[i] = b_vorh

    #print('Rand vor alpha', randStütze.t)
    #print('Eck vor alpha', eckStütze.t)

    buildingProp.t_stütze = []

    # Erhöhen der Stützenquerschnitte, wenn kleiner als erforderliches Steifigkeitsverhältnisses
    for i in buildingProp.posOut_abschnitt:
        A_stütze = 2*materialProp.E/materialProp.E*buildingProp.I[i]/(alpha_outrigger*b_total**2)                # m**2
        t_stütze = numpy.sqrt(A_stütze)*100                                         # cm
        buildingProp.t_stütze.append(t_stütze)

        while t_stütze > randStütze.t[i]:                                           # Vergrößern von t wenn nötig
            randStütze.t[i] += a                                                    # Ersetzen von t für ganzen Abschnitt
        while t_stütze > eckStütze.t[i]:                                            # Vergrößern von t wenn nötig
            eckStütze.t[i] += a                                                     # Ersetzen von t für ganzen Abschnitt
    
    #print('t_stütze', buildingProp.t_stütze)
    #print('Rand nach alpha', randStütze.t)
    #print('Eck nach alpha', eckStütze.t)


def outriggerStiffness(I_outriggerExakt, buildingProp):
    """Berechnet die effektive Biegesteifigkeit des aktuellen Outriggers
        :param I_outriggerExakt: Flächenträgheitsmoment des Outriggerquerschnitts [m**4]
        :type I_outriggerExakt: float
        """
    I_eff = I_outriggerExakt*(1+buildingProp.b_raster/buildingProp.b_raster)**3

    return I_eff

def outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp):
    """Berechnet die Federsteifigkeiten der Outrigger
        """
    # Übergabe der Gebäudebreite
    b_total = buildingProp.b_total

    I_effOutrigger = []
    K = []
    buildingProp.H_outrigger = []

    # Summen von A und I jeweils von einer Hälfte (2 Outrigger, 2 Eckstützen, 3 Randstützen)
    A_randStütze = []                                                                                   # m**2
    A_eckStütze = []
    t_randStützeOutrigger = []                                                                          # cm
    t_eckStützeOutrigger = []
    I_outriggerExakt = []                                                                               # m**4

    # Aktualisierung Querschnittswerte
    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        t_randStützeOutrigger.append(randStütze.t[n])
        t_eckStützeOutrigger.append(eckStütze.t[n])
        t_outrigger = outrigger.t[i]
        calcProfileProp(outrigger, buildingProp, materialProp, t_outrigger)
        I_aktuell = outrigger.I
        I_outriggerExakt.append(I_aktuell)
    
    # Berechnung einzelne Querschnitssflächen und Steifigkeiten
    for i in range(0,len(buildingProp.posOut)):
        #I_effOutrigger_i = outriggerStiffness(I_outriggerExakt[i], buildingProp)
        I_effOutrigger_i = I_outriggerExakt[i]      #Änderung
        I_effOutrigger.append(I_effOutrigger_i)
        buildingProp.I_effOutrigger = I_effOutrigger
        calcProfileProp(randStütze, buildingProp, materialProp, t_randStützeOutrigger[i])
        A_randStütze.append(randStütze.A_min)
        calcProfileProp(eckStütze, buildingProp, materialProp, t_eckStützeOutrigger[i])
        A_eckStütze.append(eckStütze.A_min)
    
    for i in buildingProp.posOut:
        H = (buildingProp.n-i)*buildingProp.h_geschoss                                                 # Höhe des Outriggers
        buildingProp.H_outrigger.append(H)

    for i in range(0, len(buildingProp.posOut)):
        I_outriggerGesamt = 2*I_effOutrigger[i]                                                         # 4 Outrigger in Windrichtung
        A_stützen = 3*A_randStütze[i] + 2*A_eckStütze[i]
        K_i = 12*materialProp.E*I_outriggerGesamt*A_stützen*b_total**2/(A_stützen*b_total**3+24*buildingProp.H_outrigger[i]*I_outriggerGesamt)
        K.append(K_i)

    return K                                                                                            # kNm

def calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n):
    """Berechnet die aktuellen Steifigkeitsverhältnisse zwischen Kern und Outrigger und Kern und Stütze
        :param i: aktuell zu berechnende Outriggerposition
        :type i: int
        :param n: Abschnitt der aktuell zu berechnenden Outriggerposition
        :type n: int
        """
    b_total = buildingProp.b_total
    h = buildingProp.h_total

    t=kern.t[n]
    calculations.calcProfileProp(kern, buildingProp, materialProp, t)
    I_core = kern.I

    t = outrigger.t[i]
    calculations.calcProfileProp(outrigger, buildingProp, materialProp, t)
    I_out = outrigger.I

    #beta_ist = materialProp.E/materialProp.E*I_core*b_total/(4*I_out*h)

    t = randStütze.t[n]
    calculations.calcProfileProp(randStütze, buildingProp, materialProp, t)
    A_rand = randStütze.A_min

    t = eckStütze.t[n]
    calculations.calcProfileProp(eckStütze, buildingProp, materialProp, t)
    A_eck = eckStütze.A_min

    #alpha_ist = 2*materialProp.E*I_core/(b_total**2*materialProp.E*(6*A_rand+4*A_eck))

    return beta_ist, alpha_ist

def controllStiffnesRatios(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze):
    """Überprüft die aktuellen Steifigkeitsverhältnisse zwischen Kern und Outrigger und Kern und Stütze
        Erhöht Dimension Outrigger/ Stütze / Kern bis gewolltes Verhältniss erreicht ist
        """
    a = materialProp.delta_t                                                        # cm
    Verhältnis1 = [None]*buildingProp.x
    Verhältnis2 = [None]*buildingProp.x
    beta_outrigger = buildingProp.beta_outrigger
    alpha_outrigger = buildingProp.alpha_outrigger#*10**(-3)
    fehler = 0.4

    for i in buildingProp.posOut:
        n = i//buildingProp.n_abschnitt
        beta_ist, alpha_ist = calculateRatio(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze, i, n)
        delta_beta = beta_outrigger-beta_ist
        delta_alpha = alpha_outrigger-alpha_ist

        if abs(delta_beta) <= fehler:
            Verhältnis1[n] = True
        elif delta_beta > fehler:                                                   # Erhöhung Kern
            kern.t[n] += a 
            Verhältnis1[n] = False
        elif delta_beta < -fehler:                                                  # Erhöhung Outrigger
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

    for i in buildingProp.posOut:
        delta_outrigger = outrigger.t[i] - t_check_outrigger[i]
        if abs(delta_outrigger) < fehler:
            Check_outrigger.append(True)
        else:
            Check_outrigger.append(False)

        delta_belt = belt.t[i] - t_check_belt[i]
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

    # Aktualisieren Outrigger und Stützen
    ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

    # Erneut GzT für aktualisierte QS?



def shearStiffnessModification(buildingProp, kern, randStütze, eckStütze, outrigger, delta_t):
    """Erhöht die Querschnittsdimension des Kerns, da dieser maßgebend für die Schubverformung ist
        """

    kern.t = [element + delta_t for element in kern.t]


def design(buildingProp,loads,materialProp,DataProp):                                                                                # Aufruf in app.py.MainWondow.MainCalculation
    """Führt die Bemessung  und Analyse des Outrigger-Tragwerks durch
        """
    # Anordnung der Outrigger über die Höhe nach Eingabe gui bestimmen
    arrangementOutrigger(buildingProp)

    # Gebäudegeometrie übergeben:
    b_raster=buildingProp.b_raster
    h_geschoss=buildingProp.h_geschoss
    A=b_raster**2

    # Tragwerksteile initiieren
    # mit Aussteifung
    innenStütze=building.elements(A,0,'Stütze','Vollprofil')                # Einzugsfläche, Fassadenlänge, Typ, Profil
    randStütze=building.elements(1/2*A,b_raster,'Stütze','Vollprofil')
    eckStütze=building.elements(1/4*A,b_raster,'Stütze','Vollprofil')
    kern=building.elements(8*A,0,'Kern','Kern')
    outrigger=building.elements(0,0,'Outrigger','Rechteckiges Vollprofil')
    belt = building.elements(0,0,'Belt Truss','Rechteckiges Vollprofil')

    # ohne Aussteifung
    randStützeOhnePFH=building.elements(1/2*A,b_raster,'Stütze ohne PFH','Vollprofil')
    eckStützeOhnePFH=building.elements(1/4*A,b_raster,'Stütze ohne PFH','Vollprofil')
    kernOhnePFH=building.elements(8*A,0,'Kern ohne PFH','Kern')

    # Aktuelles Tragsystem für calculation zuordnen
    str_=str_outrigger

    # Übergabe an Bemessung, dass noch keine Iteration (ohne Outrigger)
    buildingProp.Iteration = False

    # Tragfähigkeitsnachweise Vertikallasten (Für Elemente ohne Aussteifung):
    # Für Elemente ohne Aussteifung
    calculations.calcElementWidth(innenStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(randStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStützeOhnePFH,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kernOhnePFH,buildingProp,loads,materialProp,str_)
    
    # Für Elemente mit Aussteifung als Startwert
    calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(eckStütze,buildingProp,loads,materialProp,str_)
    calculations.calcElementWidth(kern,buildingProp,loads,materialProp,str_)

    # Übergabe Materialeigenschaften
    a = materialProp.delta_t
    f = materialProp.f

    # Erstellung Liste Outriggerdimension über n = Anzahl Geschosse
    outrigger.t = [0] * buildingProp.n

    # Iterationszähler initialisieren
    buildingProp.Iteration_counter = 0

    # Schleifenbedingung initialisieren
    Tragfähigkeitsnachweis = False
    t_check_eckstütze = []
     

    while Tragfähigkeitsnachweis == False:

        # Angabe nun Iteration aktiv
        buildingProp.Iteration = True

        # Erhöhen des Iterationszählers
        buildingProp.Iteration_counter += 1

        
        # Berechnung Steifigkeit des Kerns
        buildingStiffness(buildingProp,materialProp,kern,element2=None,element3=None,element4=None)   # nur core, passt das?
        
        # Berechnung/ Aktualisierung des Outriggers und der Stützen bezogen auf Steifigkeitsverhältnis zu Kern
        ratiosOutrigger(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

        # Federsteifigkeiten der Outriggergeschosse:
        buildingProp.K = outriggerSystemStiffness(randStütze, eckStütze, outrigger, buildingProp, materialProp)

        # Initialisierung fea Model für Kalkulation Momentenverlauf und Auflagerkräfte
        buildingProp.mue=[]
        for i in range (0,len(kern.t)): 
            buildingProp.mue.append(0)

        feModel = fea.feModel(buildingProp, loads, materialProp)
        
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
        #buildingProp.sigmaStütze_outrigger = []
        buildingProp.kraftStütze_outrigger = []                                         # Kraft durch Outrigger von oben nach unten
        I_outrigger  = [None]*len(buildingProp.posOut)
        I_outrigger2 = [None]*len(buildingProp.posOut)

        for n, i in enumerate(buildingProp.posOut):
            calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
            I_outrigger2[n] = outrigger.I
            I_outrigger[n] = outrigger.I#outriggerStiffness(outrigger.I, buildingProp)              # m**4

        for i in range(0, len(buildingProp.posOut)):
            moment_spring = reactions[i+3]                                              # kNm
            buildingProp.moment_spring.append(moment_spring)
            n = buildingProp.posOut_abschnitt[i]                                        # Positionen von oben nach unten
            A_stützen_ges = 6*(randStütze.t[n]/100)**2+4*(eckStütze.t[n]/100**2)
            H = buildingProp.H_outrigger[i]
            Kraft_stütze = 3*moment_spring/buildingProp.K[i]*materialProp.E*A_stützen_ges*4*I_outrigger[i]*2*b_raster\
                /(3*H*4*I_outrigger[i]+A_stützen_ges*(2*b_raster)**3)                       # kN        aus KV, ergibt höhere Last als I_outrigger,exakt und l = b_raster
            #Kraft_stütze2 = 3*moment_spring/buildingProp.K[i]*materialProp.E*A_stützen_ges*4*I_outrigger2[i]\
                #/(3*4*I_outrigger2[i]+A_stützen_ges*(b_raster)**2)
            #Kraft_stütze_alternativ = moment_spring*2/buildingProp.b_total
            #spannung_stütze = Kraft_stütze/A_stützen_ges                                # kN/m**2   Kraft/Fläche
            #buildingProp.sigmaStütze_outrigger.append(spannung_stütze)
            buildingProp.kraftStütze_outrigger.append(Kraft_stütze/10)                   # Kraft je Stütze

        
        #---------------------------------------------------------------------------------------------------
        # Tragfähigkeitsnachweise Aussteifung:
        #---------------------------------------------------------------------------------------------------

        # Nachweis Outrigger
        # Weglassen, da Beanspruchung nicht einfach betrachtbar
        # Wird als nachgewiesen angenommen

        # Moment aus Feder in fea
        #f_outrigger = 100*f
        #for n,i in enumerate(buildingProp.posOut):
            #M_max_outrigger = loads.gamma_w*reactions[n+3]
            #sigma_outrigger = 2*f_outrigger
            #while sigma_outrigger > f_outrigger:
                #M_g = outrigger.t[i]/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2/8
                #calcProfileProp(outrigger, buildingProp, materialProp, outrigger.t[i])
                #sigma_outrigger = (M_max_outrigger/4+M_g)/outrigger.W   #da 4 Outrigger je Richtung
                #if sigma_outrigger > f_outrigger:
                    #outrigger.t[i] += a
        
        
        # Nachweis Belt Truss
        # Durchlaufträger von Eck- zu Eckstütze
        # beansprucht durch Moment aus Auflagerkräften Stützen und Outrigger
        # Mmax in 1. oder 2. Feld
        belt.t = [0] * buildingProp.n
        A_Out = []
    
        for n,i in enumerate(buildingProp.posOut):
            verdrehung = reactions[2+buildingProp.n_outrigger-n]//buildingProp.K[n]                         # Phi = M/K
            A_Out.append(3/8*outrigger.t[i]/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster+\
                3*loads.gamma_w*materialProp.E*buildingProp.I_effOutrigger[n]*verdrehung/b_raster**2)       # für EG Outrigger

        for n,i in enumerate(buildingProp.posOut):
            t = materialProp.t_min
            M_belt_Out1 = 1/5*b_raster*A_Out[n]                                                             # M_max Feld 1
            M_belt_Out2 = 3/10*b_raster*A_Out[n]                                                            # M_max Feld 2
            sigma = 2*f
            while sigma > f:
                M_belt_EG1 = 0.077*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2            # kNm
                M_belt_EG2 = 0.036*t/100*h_geschoss*materialProp.gamma*loads.gamma_g*b_raster**2            # kNm
                M_max_belt = max(M_belt_EG1+M_belt_Out1, M_belt_EG2+M_belt_Out2)
                calcProfileProp(belt, buildingProp, materialProp, t)
                sigma = M_max_belt/belt.W                                                                   # kN/m**2
                t += a
            belt.t[i] = t
        
        buildingProp.t_belt_now = belt.t
        buildingProp.t_outrigger_now = outrigger.t

        t_rand_preCalc = []
        t_eck_preCalc = []

        for i in buildingProp.posOut_abschnitt:
            t_rand_preCalc.append(randStütze.t[i])
            t_eck_preCalc.append(eckStütze.t[i])

        #print('Rand vor Bemessung', t_rand_preCalc)
        #print('Eck vor Bemessung', t_eck_preCalc)

        # Nachweis Stützen
        # Beanspruchung: Auflagerkraft Feder + Vertikallasten inkl aus Federn darüber
        calculations.calcElementWidth(randStütze,buildingProp,loads,materialProp,str_)
        calculations.calcElementWidth(eckStütze, buildingProp,loads,materialProp,str_)

        t_rand_postCalc = []
        t_eck_postCalc = []
        for i in buildingProp.posOut_abschnitt:
            t_rand_postCalc.append(randStütze.t[i])
            t_eck_postCalc.append(eckStütze.t[i])

        #print('Rand nach Bemessung', t_rand_postCalc)
        #print('Eck nach Bemessung', t_eck_postCalc)
        
        

        
        # Nachweis Kern 
        # Beanspruchung: Momente aus fea + Vertikallasten
        calculations.calcElementWidth(kern, buildingProp,loads,materialProp,str_)
        
        # Kontrolle beta_ist = beta_outrigger/alpha_ist = alpha_outrigger
        # inkl. Erhöhen t mit Berücksichtigung Steifigkeitsverhältnis
        #Test = controllStiffnesRatios(buildingProp, materialProp, kern, outrigger, randStütze, eckStütze)

        if buildingProp.Iteration_counter == 1:
            Test = False
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
        for i in buildingProp.posOut:
            t_check_outrigger[i] = outrigger.t[i]
            t_check_belt[i] = belt.t[i]


        if Test == True:
            Tragfähigkeitsnachweis = True
        
        
    print('Iterationsanzahl = ', buildingProp.Iteration_counter)
        
    
    # Gebäudenachweise:
    t0=kern.t[-1]                                                                   #Kerndicke ganz unten

    Tragfähigkeitsnachweis = False
    Iteration_counter = 0

    while Tragfähigkeitsnachweis == False:
        Iteration_counter += 1
        calculations.buildingDeflection(buildingProp,loads,materialProp,str_,kern, randStütze, eckStütze, outrigger)
        calculations.interstoryDrift(buildingProp,loads,materialProp,str_,kern, randStütze, eckStütze, outrigger)
        buildingProp.t_kern=kern.t

        if Iteration_counter == 1:
            Test = False
        else:
            Test = checkChangeGebäudenachweis(materialProp, t_check_kern, kern)    

        t_check_kern = kern.t                                                       #Speicherung t-Werte für Check in nächster Iteration

        if Test == True:
            Tragfähigkeitsnachweis = True


    # Vergleich Tragfähigkeit und Gebäudenachweise:
    if kern.t[-1] > t0:     
        buildingProp.NW_maßgebend='Verformung'
        print(kern.t[-1],t0)
    else:
        buildingProp.NW_maßgebend='Tragfähigkeit'

    print('maßgebender NW', buildingProp.NW_maßgebend)

    # Massenkalkulation:                                                            # in t
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

 
    # Überprüfen ob Stütze größer als erf für alpha ist (keine Erhöhung)
    buildingProp.alpha_outriggerFulfilled = True
    for n,i in enumerate(buildingProp.posOut_abschnitt):
        if randStütze.t[i] - a > buildingProp.t_stütze[n]:
            buildingProp.alpha_outriggerFulfilled = False
        elif eckStütze.t[i] - a > buildingProp.t_stütze[n]:
            buildingProp.alpha_outriggerFulfilled = False

    print('Steifigkeitsverhältnis alpha zwischen Kern und Stützen ist ', buildingProp.alpha_outriggerFulfilled)

    print('M', buildingProp.moment_spring)
    print('N_Feder', buildingProp.kraftStütze_outrigger)
 