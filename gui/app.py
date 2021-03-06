# ------------------------------------------------------------------------------
# Description:  Startet die GUI, andere Funktionen werden in der GUI aufgerufen
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-23      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke

# Co-Author:    Carola Grupp
# Created:      2021-11-02  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:      https://www.ilek.uni-stuttgart.de/
# ------------------------------------------------------------------------------
# Imports:
# PyQt5 gui
from matplotlib.pyplot import vlines, xlim, ylim
from numpy.lib.function_base import append
import numpy as np
from Masterarbeit.pfh.building import buildingProp, loads, materialProp
from PyQt5 import QtWidgets, uic
from gui.mainWindow import Ui_MainWindow

# python modules
import sys
import os
import numpy

# pyLEK/helpers
import pyLEK.plotters.plot2D as plt
import pyLEK.plotters.plot2D as plot2D
from pyLEK.sampleCode.gui.pdHelpers import readDataframe
from pyLEK.sampleCode.gui.pdHelpers import getAvailStrengthClasses
from pyLEK.sampleCode.gui.pdHelpers import getAvailMaterials
from pyLEK.sampleCode.gui.pdHelpers import getMaterialProperties

from pfh import building
from WindData import data
from pfh import str_core
from pfh import str_frame
from pfh import str_bracedtube
from pfh import str_framedTube
from pfh import str_outrigger

# ------------------------------------------------------------------------------


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, *args, **kwargs):
        # Here we declare that the MainWindow class inherits from
        # QtWidgets.QMainWindow, Ui_MainWindow
        super().__init__(*args, **kwargs)

        # Initialize gui
        self.gui = Ui_MainWindow()

        # Erstellen der Klassen
        buildingProp = building.buildingProp()
        materialProp = building.materialProp()
        loads = building.loads()
        DataProp = data.DataProp()

        # Setup gui
        self.gui.setupUi(self)

        # Abschnitt Multiberechnung
        self.gui.comboBox_parameter.currentTextChanged.connect(self.fitUnits)
        self.gui.comboBox_parameter_2.currentTextChanged.connect(self.fitUnits)

        self.gui.comboBox_Nachweis.currentTextChanged.connect(lambda: self.checkNachweis(buildingProp)) #pr??ft ob der Tragf??higkeitsnachweis bei dem Outriggertragwerk gef??hrt werden soll

        # Abschnitt Geb??udegeometrie (Zusammenh??nge des Tabs)
        self.gui.spinBox_h_geschoss.valueChanged.connect(self.calcTotalBuildingHeight) #aktualisiert die Gesamth??he bei ??nderung der Geschossh??he
        self.gui.spinBox_n.valueChanged.connect(self.calcTotalBuildingHeight) #aktualisiert die Gesamth??he bei ??nderung der Geschossanzahl

        self.gui.spinBox_b_raster.valueChanged.connect(self.calcTotalBuildingWidth) #aktualisiert die Gesamtbreite bei ??nderung des St??tzenabstands

        self.gui.spinBox_h_total.valueChanged.connect(self.calcBuildingSlenderness) #aktualisiert Schlankheit bei ??nderung der Gesamth??he
        self.gui.spinBox_b_total.valueChanged.connect(self.calcBuildingSlenderness) #aktualisiert Schlankheit bei ??nderung der Gesamtbreite

        self.gui.comboBox_tragwerk.currentTextChanged.connect(self.activateStielanzahl) #aktualisiert Anzahl Stiele pro Raster bei ??nderung des Tragwerks
        self.gui.comboBox_tragwerk.currentTextChanged.connect(self.activateOutriggeranzahl)    #aktualisiert Anzahl Outrigger bei ??nderung des Tragwerks
        
        self.gui.comboBox_tragwerk.currentTextChanged.connect(self.activateOutriggerVerhaeltnisse)

        #if self.gui.checkBox_170_02.isChecked:
            #self.gui.spinBox_n.valueChanged.connect(self.updateGeometrie)
            #self.gui.spinBox_b_raster.valueChanged.connect(self.updateGeometrie)



        # Abschnitt Vertikallasten
        # Aktualisierung der Designlast bei ??nderungen der charakteristischen Teillasten
        #Eigengewicht
        self.gui.spinBox_gk1.valueChanged.connect(self.calcDesignLoadG)
        self.gui.spinBox_gk2.valueChanged.connect(self.calcDesignLoadG)
        self.gui.spinBox_gamma_g.valueChanged.connect(self.calcDesignLoadG)

        #Fassade
        self.gui.spinBox_gk_fassade.valueChanged.connect(self.calcDesignLoadFassade)
        self.gui.spinBox_gamma_g.valueChanged.connect(self.calcDesignLoadFassade)

        #Nutzlast
        self.gui.spinBox_qk.valueChanged.connect(self.calcDesignLoadQ)
        self.gui.spinBox_gamma_q.valueChanged.connect(self.calcDesignLoadQ)
        
        # Abschnitt Horizontallasten
        self.gui.spinBox_GK.valueChanged.connect(self.setExponent) #akualisiert Profilexponent bei ??nderung Gel??ndekategorie

        self.gui.spinBox_cf0.valueChanged.connect(self.calcCF) #akt. c_f bei ??nderung c_f,0
        self.gui.spinBox_psi_r.valueChanged.connect(self.calcCF)
        self.gui.spinBox_psi_lambda.valueChanged.connect(self.calcCF)

        self.gui.spinBox_qb_k.valueChanged.connect(self.calcWindLoad10) #Aktualisierung der charakteristischen Windlast auf 10m  H??he
        self.gui.spinBox_cscd.valueChanged.connect(self.calcWindLoad10)
        self.gui.spinBox_cf.valueChanged.connect(self.calcWindLoad10)
        self.gui.spinBox_GK.valueChanged.connect(self.calcWindLoad10)

        self.gui.spinBox_wk.valueChanged.connect(self.calcWindLoad10Design) #Aktualisierung Designwert Windlast
        self.gui.spinBox_gamma_w.valueChanged.connect(self.calcWindLoad10Design)

        # Abschnitt Windkanaldaten

        #self.gui.checkBox_170_02.stateChanged.connect(self.activateWindkanal)
        self.gui.comboBox_170_23.currentTextChanged.connect(self.calcAlpha_v)
        self.gui.comboBox_170_13.currentTextChanged.connect(self.updateGeometrie)

        # Zum berechnen: Push button zum Starten
        self.gui.pushButton_einzelberechnung.clicked.connect(lambda: self.pushButton_einzelberechnung(buildingProp,materialProp,loads,DataProp))
        self.gui.pushButton_multiberechnung.clicked.connect(lambda: self.pushButton_multiberechnung(buildingProp,materialProp,loads,DataProp))
        self.gui.pushButton_multiberechnungParameter.clicked.connect(lambda: self.pushButton_multiberechnungParameter(buildingProp,materialProp,loads,DataProp))
        self.gui.pushButton_multiberechnungParameter_2D.clicked.connect(lambda: self.pushButton_multiberechnungParameter_2D(buildingProp,materialProp,loads,DataProp))

        # Zum plotten: Plot bei ??nderung
        # Einzelberechnung
        self.gui.comboBox_yAxis_222.currentTextChanged.connect(lambda: self.plotMassDistribution(buildingProp))

        self.gui.comboBox_xAxis_231.currentTextChanged.connect(lambda: self.plotStiffnessDistribution(buildingProp,materialProp))
        self.gui.comboBox_yAxis_232.currentTextChanged.connect(lambda: self.plotStiffnessDistribution(buildingProp,materialProp))

        self.gui.comboBox_xAxis_241.currentTextChanged.connect(lambda: self.plotInternalForcesCurve(buildingProp,loads))
        self.gui.comboBox_yAxis_242.currentTextChanged.connect(lambda: self.plotInternalForcesCurve(buildingProp,loads))
        
        self.gui.comboBox_xAxis_251.currentTextChanged.connect(lambda: self.plotDeformationCurve(buildingProp))
        self.gui.comboBox_yAxis_252.currentTextChanged.connect(lambda: self.plotDeformationCurve(buildingProp))

        self.gui.comboBox_xAxis_271.currentTextChanged.connect(lambda: self.plotAusnutzung(buildingProp, materialProp))
        self.gui.comboBox_yAxis_272.currentTextChanged.connect(lambda: self.plotAusnutzung(buildingProp, materialProp))

        # H??henvariation
        self.gui.comboBox_xAxis_311.currentTextChanged.connect(lambda: self.plotResourceAnalysis(buildingProp,materialProp))
        self.gui.comboBox_yAxis_312.currentTextChanged.connect(lambda: self.plotResourceAnalysis(buildingProp,materialProp))

        self.gui.comboBox_xAxis_321.currentTextChanged.connect(lambda: self.plotResourceAnalysis_m2(buildingProp,materialProp))
        self.gui.comboBox_yAxis_322.currentTextChanged.connect(lambda: self.plotResourceAnalysis_m2(buildingProp,materialProp))

        self.gui.comboBox_xAxis_331.currentTextChanged.connect(lambda: self.plotPFHAnalysis(buildingProp,materialProp))
        self.gui.comboBox_yAxis_332.currentTextChanged.connect(lambda: self.plotPFHAnalysis(buildingProp,materialProp))

        self.gui.comboBox_xAxis_341.currentTextChanged.connect(lambda: self.plotPFHAnalysis_m2(buildingProp,materialProp))
        self.gui.comboBox_yAxis_342.currentTextChanged.connect(lambda: self.plotPFHAnalysis_m2(buildingProp,materialProp))

        self.gui.comboBox_xAxis_351.currentTextChanged.connect(lambda: self.plotEigenFrequency(buildingProp))
        self.gui.comboBox_yAxis_352.currentTextChanged.connect(lambda: self.plotEigenFrequency(buildingProp))

        self.gui.comboBox_xAxis_361.currentTextChanged.connect(lambda: self.plotMoment(buildingProp))
        self.gui.comboBox_yAxis_362.currentTextChanged.connect(lambda: self.plotMoment(buildingProp))

        # Parametervariation
        self.gui.comboBox_xAxis_411.currentTextChanged.connect(lambda: self.plotParameterAnalysis(buildingProp, materialProp))
        self.gui.comboBox_yAxis_412.currentTextChanged.connect(lambda: self.plotParameterAnalysis(buildingProp, materialProp))

        self.gui.comboBox_xAxis_414.currentTextChanged.connect(lambda: self.plotDeformation(buildingProp, materialProp))
        self.gui.comboBox_yAxis_414.currentTextChanged.connect(lambda: self.plotDeformation(buildingProp, materialProp))

        self.gui.comboBox_xAxis_516.currentTextChanged.connect(lambda: self.plotAlphaIst(self,buildingProp,materialProp))
        self.gui.comboBox_yAxis_517.currentTextChanged.connect(lambda: self.plotAlphaIst(self,buildingProp,materialProp))
        
        self.gui.comboBox_xAxis_421.currentTextChanged.connect(lambda: self.plotOptimalParameter(buildingProp,materialProp))
        self.gui.comboBox_yAxis_422.currentTextChanged.connect(lambda: self.plotOptimalParameter(buildingProp,materialProp))

        self.gui.comboBox_xAxis_431.currentTextChanged.connect(lambda: self.plotPFHAnalysisOptimalParameter_m2(buildingProp,materialProp))
        self.gui.comboBox_yAxis_432.currentTextChanged.connect(lambda: self.plotPFHAnalysisOptimalParameter_m2(buildingProp,materialProp))

        self.gui.comboBox_xAxis_441.currentTextChanged.connect(lambda: self.plotPFHAnalysisParameterComparisson_m2(buildingProp,materialProp))
        self.gui.comboBox_yAxis_442.currentTextChanged.connect(lambda: self.plotPFHAnalysisParameterComparisson_m2(buildingProp,materialProp))

        self.gui.comboBox_xAxis_451.currentTextChanged.connect(lambda: self.plotInfluenceOptimalParameter_m2(buildingProp,materialProp))
        self.gui.comboBox_yAxis_452.currentTextChanged.connect(lambda: self.plotInfluenceOptimalParameter_m2(buildingProp,materialProp))
        
        self.gui.comboBox_parameter.currentTextChanged.connect(lambda: self.checkIfPosOutVariing(buildingProp))
        self.gui.comboBox_parameter_2.currentTextChanged.connect(lambda: self.checkIfPosOutVariing_2D(buildingProp))

        self.gui.comboBox_darstellung_473.currentTextChanged.connect(lambda: self.plot3DParameter(buildingProp,materialProp))
        self.gui.comboBox_zAxis_472.currentTextChanged.connect(lambda: self.plot3DParameter(buildingProp,materialProp))
        
        
        # Save values from selected PyQt5.QtWidgets
        # See qtHelpers for further details
        self.gui.actionSpeichern.triggered.connect(self.guiSaveState)

        # Restore values from selected PyQt5.QtWidgets
        # See qtHelpers for further details
        self.gui.actionLaden.triggered.connect(self.guiRestoreState)

        # Example: Linking with pandas-dataframe (Using .csv as a database)
        # 1) Load dataframe
        # Import from pdHelpers

        # Change to current file location, avoiding abs-path for link to database
        os.chdir(os.path.dirname(sys.argv[0]))

        # Reading from database
        self.df = readDataframe('gui/materials.csv')

        # 2) Initialize at startup
        self.addMaterials()
        self.modifyStrengthClasses()
        self.modifyMaterialProperties()

        # 3) Link signals
        self.gui.comboBoxMaterial.textActivated.connect(
            self.modifyStrengthClasses)
        self.gui.comboBoxStrengthClass.textActivated.connect(
            self.modifyMaterialProperties)

    # Speichern und Laden:

    def guiSaveState(self):
        # Import from qtHelpers
        from pyLEK.sampleCode.gui import qtHelpers

        # Execute imported method to save state
        qtHelpers.guiSaveState(self.gui)

    def guiRestoreState(self):
        # Import from qtHelpers
        from pyLEK.sampleCode.gui import qtHelpers

        # Execute imported method to restore state
        qtHelpers.guiRestoreState(self.gui)


    # Materialeigenschaften mittels Dropdownauswahl:

    def addMaterials(self):
        # Import from pdHelpers
        #from gui.pdHelpers import getAvailMaterials

        # Filtering available materials
        materials = getAvailMaterials(self.df)

        # Set up ComboBox
        self.gui.comboBoxMaterial.clear()
        self.gui.comboBoxMaterial.addItems(materials)

    def modifyStrengthClasses(self):
        # Import from pdHelpers
        #from gui.pdHelpers import getAvailStrengthClasses

        # Get selected material
        material = self.gui.comboBoxMaterial.currentText()

        # Get strengthClasses from selected materials
        strengthClasses = getAvailStrengthClasses(self.df, material)

        # Set up ComboBox
        self.gui.comboBoxStrengthClass.clear()
        self.gui.comboBoxStrengthClass.addItems(strengthClasses)

        # Change material properties as well
        self.modifyMaterialProperties()

    def modifyMaterialProperties(self):
        # Import from pdHelpers
        #from gui.pdHelpers import getMaterialProperties

        # Get selected strengthClass
        strengthClass = self.gui.comboBoxStrengthClass.currentText()

        # Get material properties from selected strengthClass
        matProperties = getMaterialProperties(self.df, strengthClass)

        # Set up spinBoxes
        self.gui.spinBox_f.setValue(matProperties["fd"]/10)
        self.gui.spinBox_E.setValue(matProperties["E-Modul"])
        self.gui.spinBox_G.setValue(matProperties["G-Modul"])
        self.gui.spinBox_gamma.setValue(matProperties["Wichte"]/10)

    #-------------------------------------------------------------------------------------------------------------------------------------------------
    # Automatische Berechnungen der Angaben im Fenster (ben??tigte Funktionen):
    
    # Multiberechnung (verschiedene Parametervariationen je nach Thematik):
    def fitUnits(self):
        parameter=self.gui.comboBox_parameter.currentText()
        parameter2D = self.gui.comboBox_parameter_2.currentText()

        if parameter == 'Keiner':
            label=''
            self.gui.pushButton_multiberechnungParameter.setEnabled(False)
        
        else:
           self.gui.pushButton_multiberechnungParameter.setEnabled(True) 
        
        if parameter2D == 'Keiner':
            label=''
            self.gui.pushButton_multiberechnungParameter_2D.setEnabled(False)
        
        else:
           self.gui.pushButton_multiberechnungParameter_2D.setEnabled(True)

        if parameter == 'Staffelung n_abschnitt' or parameter == 'Anzahl Stiele pro Raster (Framed Tube) n_stiele/b_raster':
            label='[-]'
            self.gui.spinBox_p_min.setMinimum(1)
            self.gui.spinBox_p_max.setMinimum(1)

        elif parameter == 'Outriggeranordnung ??ber die H??he':
            self.gui.spinBox_p_min.setMinimum(0)
            self.gui.spinBox_p_max.setMaximum(1)

        else:
            self.gui.spinBox_p_min.setMinimum(0) #keine Variation der Staffelung
            self.gui.spinBox_p_max.setMinimum(0)

        if parameter == 'Deckengewicht g_k1':
            label='[kN/m??]' #f??gt die Einheit hinter Eingabebereich hinzu

        if parameter == 'Windlast q_b':
            label='[kN/m??]'
           
        if parameter == 'Geb??udeform psi_r':
            label='[-]'

        if parameter == 'Verformungsverh??ltnis w_EI/w_GA':
            label='[-]'

        if parameter == 'Mindestwirksamkeit Framed Tube tube_wirksam':
            label='[%]'
        
        if parameter == 'Anzahl Outrigger':
            label = '[-]'

        if parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
            label = '[-]'

        if parameter == 'Steifigkeitsverh??ltnis Beta (Outrigger)':
            label = '[-]'
        
        if parameter == 'Outriggeranordnung ??ber die H??he':
            label = '[Vorfaktor f??r H??he ab Fu??punkt]'

        self.gui.label_p_min.setText(label)
        self.gui.label_p_max.setText(label)

        if parameter2D == 'Steifigkeitsverh??ltnis Alpha und Beta (Outrigger)':
            label = '[-]'
        if parameter2D == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
            label = '[Vorfaktor f??r H??he ab Fu??punkt]'

        self.gui.label_p_min_1.setText(label)
        self.gui.label_p_min_2.setText(label)
        self.gui.label_p_max_1.setText(label)
        self.gui.label_p_max_2.setText(label)

    # Geometrie         Funktionen f??r Geometrieberechnung

    def checkNachweis(self, buildingProp):

        if self.gui.comboBox_Nachweis.currentText() == 'Ja':
            buildingProp.Nachweis_GZT = False

        else:
            buildingProp.Nachweis_GZT = True


    def calcTotalBuildingHeight(self):
        """Calculates the total building height
        """                                     #Mehrzeiliger String
        #if self.gui.checkBox_170_02.isChecked:
            #self.updateGeometrie

        n = self.gui.spinBox_n.value()
        h_geschoss = self.gui.spinBox_h_geschoss.value()
        h_total = n * h_geschoss
        self.gui.spinBox_h_total.setValue(h_total)

    def calcTotalBuildingWidth(self):
        """Calculates the total building height
        """
        b_raster = self.gui.spinBox_b_raster.value()
        b_total = 4 * b_raster                  #St??tzenraster 5x5, also 4 Felder mit b_raster
        self.gui.spinBox_b_total.setValue(b_total)

    def calcBuildingSlenderness(self):
        """Calculates the building's slenderness
        """
        h_total = self.gui.spinBox_h_total.value()
        b_total = self.gui.spinBox_b_total.value()
        if b_total==0:              #damit nicht /0
            lambda_slenderness=0
        else:    
            lambda_slenderness = h_total / b_total
        self.gui.spinBox_lambda.setValue(lambda_slenderness)

    def activateStielanzahl (self):
        ''''Aktiviert Stielanzahl f??r Framed Tube'''
        if self.gui.comboBox_tragwerk.currentText() == 'Framed Tube':
            self.gui.spinBox_n_stiele.setEnabled(True)

        else:
            self.gui.spinBox_n_stiele.setEnabled(False)
            self.gui.spinBox_n_stiele.setValue(1)

    def activateOutriggeranzahl (self):
        ''''Aktiviert Outriggeranzahl f??r Outriggersystem'''
        if self.gui.comboBox_tragwerk.currentText() == 'Outrigger':
            self.gui.spinBox_n_outrigger.setEnabled(True)
        else:
            self.gui.spinBox_n_outrigger.setEnabled(False)
            self.gui.spinBox_n_outrigger.setValue(1)

    def activateOutriggerVerhaeltnisse(self):
        ''''Aktiviert Steifgkeitsverh??ltnisse f??r Outriggersystem'''
        if self.gui.comboBox_tragwerk.currentText() == 'Outrigger':
            self.gui.spinBox_alpha_outrigger.setEnabled(True)
            self.gui.spinBox_beta_outrigger.setEnabled(True)
        else:
            self.gui.spinBox_alpha_outrigger.setEnabled(False)
            self.gui.spinBox_beta_outrigger.setEnabled(False)
            self.gui.spinBox_alpha_outrigger.setValue(20)
            self.gui.spinBox_beta_outrigger.setValue(20)

    # Vertikallasten

    def calcDesignLoadG(self):
        """Calculates the design load of g
        """
        g_k1 = self.gui.spinBox_gk1.value()
        g_k2 = self.gui.spinBox_gk2.value()
        gamma_g = self.gui.spinBox_gamma_g.value()
        g_d = (g_k1+g_k2)*gamma_g
        self.gui.spinBox_gd.setValue(g_d)

    def calcDesignLoadFassade(self):
        """Calculates the design load of the fassade
        """
        g_k_fassade = self.gui.spinBox_gk_fassade.value()
        gamma_g = self.gui.spinBox_gamma_g.value()
        g_d_fassade = g_k_fassade*gamma_g
        self.gui.spinBox_gd_fassade.setValue(g_d_fassade)

    def calcDesignLoadQ(self):
        """Calculates the design load of q
        """
        q_k = self.gui.spinBox_qk.value()
        gamma_q = self.gui.spinBox_gamma_q.value()
        q_d = q_k*gamma_q
        self.gui.spinBox_qd.setValue(q_d)

    # Horizontallasten

    def setExponent(self): #Profilexponent 
        """Sets the Exponent
        """
        GK = self.gui.spinBox_GK.value()
        if GK == 0:
            exp=1   
        
        if GK == 1:
            exp=0.19            #aus B??engeschwindigkeitsdruck (EC1-1-4 Tab Na B.2)
        
        if GK == 2:
            exp=0.24

        if GK == 3:
            exp=0.31

        if GK == 4:
            exp=0.4

        self.gui.spinBox_exp.setValue(exp)

    def calcCF(self):               #Verbessern durch implementieren der Beiwerte nach Norm, im Moment h??ndische Eingabe, Schneider 8.2.1
        """Calculates the "Kraftbeiwert"
        """
        c_f0= self.gui.spinBox_cf0.value()
        psi_r= self.gui.spinBox_psi_r.value()
        psi_lambda= self.gui.spinBox_psi_lambda.value()
        c_f=c_f0*psi_r*psi_lambda
        self.gui.spinBox_cf.setValue(c_f)

    def calcWindLoad10(self):
        """Calculates the load wk at 10m height
        """
        c_f = self.gui.spinBox_cf.value()
        cscd = self.gui.spinBox_cscd.value()
        q_b = self.gui.spinBox_qb_k.value()

        GK = self.gui.spinBox_GK.value()        #x-Werte nach EC1 NA Tabelle NA.B.2
        if GK == 0:
            x=1
        
        if GK ==1:
            x=2.6

        if GK ==2:
            x=2.1

        if GK ==3:
            x=1.6

        if GK ==4:
            x=1.1

        w_k = cscd*c_f*x*q_b    #EC1 5.3 z_e = 10m vereinfacht erstmal
        self.gui.spinBox_wk.setValue(w_k)

    def calcWindLoad10Design(self):     #F??r Anzeige w_k und w_d in GUI, dann ??bergabe w_k an wk (s. building.py)
        """Calculates the design load of wk
        """
        w_k = self.gui.spinBox_wk.value()
        gamma_w = self.gui.spinBox_gamma_w.value()
        w_d=w_k*gamma_w
        self.gui.spinBox_wd.setValue(w_d)

    # Windkanaldaten

    def calcAlpha_v(self):
        """Updates Profilexponent bei ??nderung der Gel??ndekategorie
        """
        GK = self.gui.comboBox_170_23
        if GK == 2:
            alpha_v = 0.16
        elif GK == 4:
            alpha_v = 0.30
        self.gui.doubleSpinBox_170_33.setValue(alpha_v)

    def activateWindkanal(self):
        ''''Aktiviert Eingabe Windkanaldaten bei ausgew??hlter Checkbox'''
        if self.gui.checkBox_170_02.isChecked() == True:
            self.gui.comboBox_170_13.setEnabled(True)
            self.gui.comboBox_170_23.setEnabled(True)
            self.gui.doubleSpinBox_170_43.setEnabled(True)
            self.gui.doubleSpinBox_170_53.setEnabled(True)
            #self.gui.spinBox_h_geschoss.setEnabled(False)
            #self.gui.spinBox_lambda.setValue(self.gui.comboBox_170_13)
            #self.updateGeometrie()
        else:
            self.gui.comboBox_170_13.setEnabled(False)
            self.gui.comboBox_170_23.setEnabled(False)
            self.gui.doubleSpinBox_170_43.setEnabled(False)
            self.gui.doubleSpinBox_170_43.setValue(20.00)
            self.gui.doubleSpinBox_170_53.setEnabled(False)
            self.gui.doubleSpinBox_170_53.setValue(0.00)
            self.gui.spinBox_h_geschoss.setEnabled(True)

    def updateGeometrie(self):
        ''''Aktualisiert Geb??udegeometrie f??r Windkanaldaten bei ??nderung Schlankheit etc'''
        b_total = self.gui.spinBox_b_total.value()
        schlankheit = int(self.gui.comboBox_170_13.currentText())
        n = self.gui.spinBox_n.value()
        h_total = b_total*schlankheit
        h_geschoss = h_total/n
        self.gui.spinBox_h_geschoss.setValue(h_geschoss)
        # auch aufrufen, wenn sich b_raster, schlankheit, n ??ndert



    # Profilausgabe (Statische Kennwerte -> Querschnitte):

    def showProfiles(self,buildingProp):
        self.gui.label_212_01.setText(str(buildingProp.t_innenSt??tzen))
        self.gui.label_212_02.setText(str(buildingProp.t_randSt??tzen))
        self.gui.label_212_03.setText(str(buildingProp.t_eckSt??tzen))
        self.gui.label_212_04.setText(str(buildingProp.t_kern))
        self.gui.label_212_05.setText(str(buildingProp.t_diagonale))
        self.gui.label_212_06.setText(str(buildingProp.t_querstrebe))
        self.gui.label_212_07.setText(str(buildingProp.t_riegel))
        self.gui.label_212_08.setText(str(buildingProp.t_outrigger))
        self.gui.label_212_10.setText(str(buildingProp.t_belt))
        self.gui.label_212_09.setText(str(buildingProp.eigenFrequenz))
        self.gui.label_212_11.setText(str(buildingProp.w*1000))
        self.gui.label_212_12.setText(str(buildingProp.NW_ma??gebend))
        
    
    # Plots:
        # In the class MplWidget pyLEK-plotters are imported. Therefore all objects
        # of the class MplWidget can use the pyLEK-plotters methods
        # Class MplWidget
        # -> Object MplWidget_XXX

        # See also mplwidget.py for more information

        # Using the pyLEK-plotter functions
        # Die Plot2D und barChart Funktion in mplwidget.py ist schon angepasst, pieChart nicht


    def plotMassDistribution(self,buildingProp): #Statische Kennwerte -> Massenverteilung: Balkendiagramm

        # X-Achse:
        x1='Aussteifungselement'    

        if buildingProp.tragwerk == 'Kerntragwerk':
            x1='Kern'
        if buildingProp.tragwerk == 'Braced Tube':
            x1='Diagonalen'
        if buildingProp.tragwerk == 'Rahmentragwerk':
            x1='Riegel' 
        if buildingProp.tragwerk == 'Outrigger':
            x1 ='Kern, Outrigger u. Belt'
        
        # Y-Achse:
        y1=buildingProp.G_aussteifung[-1]
        y2=buildingProp.G_au??enSt??tzen[-1]
        y3=buildingProp.G_innenSt??tzen[-1]
        y4=buildingProp.G_total[-1]
        y5=buildingProp.G_decken       
        
        if self.gui.comboBox_yAxis_222.currentText() == 'Masse':
            y=[[y1,y2,y3,y4,y5]]
            ylabel='Masse [t]'
            ymax=1.1*max(y4,y5)
        
        if self.gui.comboBox_yAxis_222.currentText() == 'Anteil an Gesamtmasse':
            gesamtMasse=y4+y5
            if gesamtMasse == 0:
                y=[[0,0,0,0,0]]
            else:
                y=[[100*y1/(gesamtMasse),100*y2/(gesamtMasse),100*y3/(gesamtMasse),100*y4/(gesamtMasse),100*y5/(gesamtMasse)]]

            ymax=1.1*max(y4,y5)*100/gesamtMasse
            ylabel='Anteil an Gesamtmasse [%]'

        # Legende:
        xticklabels=['',x1,'Au??enst??tzen','Innenst??tzen', 'Tragwerk o. Decken','Decken']
        xlabel=None
        annotations='individual'
        annotations_position='center'
        title=None#'Massenverteilung'
        legend=['Masse']

        xlim=[-0.5,4.5]
        ylim=[0,ymax]

        # Speichern:        
        dir_fileName = "Massenverteilung"
        saveTex=False
        savePlt=True
        

        self.gui.MplWidget_220.plotBarChart(y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend, xticklabels=xticklabels, 
                annotations=annotations, annotations_position=annotations_position, xlim=xlim, ylim=ylim, dir_fileName=dir_fileName, vLines=None, savePlt=savePlt, saveTex=saveTex)
    

    def plotStiffnessDistribution(self,buildingProp,materialProp):      #Statische Kennwerte -> Steifigkeitsverteilung: Kurvendiagramm

        # Y-Achse:
        y=[]
        if self.gui.comboBox_yAxis_232.currentText() == 'Geb??udeh??he':
            #y.append(buildingProp.h_total)
            for i in range (0,len(buildingProp.I)):
                a=max(buildingProp.h_total-i*buildingProp.n_abschnitt*buildingProp.h_geschoss,0)
                y.append(a)
                a=max(buildingProp.h_total-(i+1)*buildingProp.n_abschnitt*buildingProp.h_geschoss,0)
                y.append(a)
            ylabel = 'H??he [m]'


        if self.gui.comboBox_yAxis_232.currentText() == 'Anzahl Stockwerke':
            #y.append(buildingProp.n)
            for i in range (0,len(buildingProp.I)):
                a=max(buildingProp.n-i*buildingProp.n_abschnitt,0)
                y.append(a)
                a=max(buildingProp.n-(i+1)*buildingProp.n_abschnitt,0)
                y.append(a)
            ylabel = 'Anzahl Stockwerke [-]'
        
        # X-Achse:
        x=[]
        if self.gui.comboBox_xAxis_231.currentText() == 'Biegesteifigkeit':
            x = 2*[i * materialProp.E/1000 for i in buildingProp.I] #2x: um die Werte zu verdoppeln, da Plot alle Werte doppelt braucht, buildingProp.I = I f??r jeden Abschnitt ??ber die H??he aus berechnetem Querschnitt
            #x.append(0)
            x.sort()
            if buildingProp.tragwerk == 'Outrigger':
                xlabel = 'Biegesteifigkeit Kern'
            else:
                xlabel = 'Biegesteifigkeit'
            legend = ['EI [MNm??]']
            title = 'Biegesteifigkeitsverteilung'
        
        if self.gui.comboBox_xAxis_231.currentText() == 'Schubsteifigkeit':
            x = 2*buildingProp.GA
            #x.append(0)
            x.sort()
            if buildingProp.tragwerk == 'Outrigger':
                xlabel = 'Schubsteifigkeit Kern'
            else:
                xlabel = 'Schubsteifigkeit'
            legend = ['GA [kN]']
            title = 'Schubsteifigkeitsverteilung'

        title=None
        #mpl='plotStyle_legendeRechts'
        maxWert=max(x)
        minWert=min(x)
        bereich=maxWert-minWert
        xlim=[minWert-0.1*bereich,maxWert+0.1*bereich]

        # Speichern als pdf_tex und als png:
        # Change to current file location
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Steifigkeitsverteilung"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_230.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, xlim=xlim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)


    def plotInternalForcesCurve(self, buildingProp, loads):    
        
        # Y-Achse:
        if self.gui.comboBox_yAxis_242.currentText() == 'Geb??udeh??he':
            if buildingProp.tragwerk == 'Outrigger':
                y = [buildingProp.h_total]
                for i in range(1, buildingProp.n):
                    y.append(buildingProp.h_total - buildingProp.h_geschoss*i)
                    y.append(buildingProp.h_total - buildingProp.h_geschoss*i)
                y.append(buildingProp.h_total - buildingProp.h_geschoss*buildingProp.n)
            else:
                y=numpy.arange(buildingProp.h_total,-buildingProp.h_geschoss,-buildingProp.h_geschoss)  #generiert  ndarray mit Staffelung ??ber Geb??udeh??he, beginnt oben
            ylabel = 'H??he [m]'
            
        if self.gui.comboBox_yAxis_242.currentText() == 'Anzahl Stockwerke':
            if buildingProp.tragwerk == 'Outrigger':
                y = [buildingProp.n]
                for i in range(1, buildingProp.n):
                    y.append(buildingProp.n - i)
                    y.append(buildingProp.n - i)
                y.append(buildingProp.n - buildingProp.n)
            else:
                y=numpy.arange(buildingProp.n,-1,-1)
            ylabel = 'Anzahl Stockwerke [-]'

        # X-Achse:
        if self.gui.comboBox_xAxis_241.currentText() == 'Momentenverlauf':
            if buildingProp.tragwerk == 'Outrigger':
                x = [-i/1000 for i in buildingProp.moment_ges]
            else:
                x = [i/1000 for i in loads.M]
            xlabel = 'Moment [MNm]'
            title = 'Momentenverlauf'
            legend = ['M']

        if self.gui.comboBox_xAxis_241.currentText() == 'Querkraftverlauf':
            if buildingProp.tragwerk == 'Outrigger':
                x = []
                for i in buildingProp.querkraft:
                    x.append(i/1000)
                    x.append(i/1000)
                #x = 2*[i/1000 for i in buildingProp.querkraft]
            else:
                x = [i/1000 for i in loads.V]
            xlabel = 'Querkraft [MN]'
            title = 'Querkraftverlauf'
            legend = ['V']
       
        #mpl='plotStyle_legendeRechts'
        title=None

        # Speichern als pdf_tex und als png:
        # Change to current file location
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Schnittgr????enverlauf"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_240.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)

    def plotDeformationCurve(self,buildingProp):
        
        # Y-Achse:
        if self.gui.comboBox_yAxis_252.currentText() == 'Geb??udeh??he':
            y1=numpy.arange(buildingProp.h_total,-buildingProp.h_geschoss,-buildingProp.h_geschoss)
            ylabel = 'H??he [m]'
            
        if self.gui.comboBox_yAxis_252.currentText() == 'Anzahl Stockwerke':
            y1=numpy.arange(buildingProp.n,-1,-1)
            ylabel = 'Anzahl Stockwerke [-]'

        # X-Achse:
        if self.gui.comboBox_xAxis_251.currentText() == 'Verformungen':
            x1=[]
            for i in range (0, len(buildingProp.w_EI)): 
                x1.append((buildingProp.w_EI[i]+buildingProp.w_GA[i])*1000)
            x2=[w_EI *1000 for w_EI in buildingProp.w_EI]   #mm, da m*1000
            x3=[w_GA *1000 for w_GA in buildingProp.w_GA]

            xlabel = 'Verformung [mm]'
            title = 'Verformungen'

        if self.gui.comboBox_xAxis_251.currentText() == 'Interstory-Drift':
            x1=[]
            for i in range (0, len(buildingProp.Teta_i_EI)):
                x1.append((buildingProp.Teta_i_EI[i]+buildingProp.Teta_i_GA[i])*1000)
            x2=[Teta_i_EI *1000 for Teta_i_EI in buildingProp.Teta_i_EI]
            x3=[Teta_i_GA *1000 for Teta_i_GA in buildingProp.Teta_i_GA]

            xlabel = 'Verformung pro Geschoss [mm]'
            title = 'Interstory-Drift'       
        
        print('w_EI=',(buildingProp.w_EI[0])*1000)
        print('w_GA=',(buildingProp.w_GA[0])*1000)
        print('w_max=',(buildingProp.w_EI[0]+buildingProp.w_GA[0])*1000)

        # Legende:
        
        legend = ['w','w_{EI}','w_{GA}']       

        x=[x1,x2,x3]
        y=[y1,y1,y1]

        mpl='plotStyle_plot2D_small_legendLowerRight'
        title=None

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Verformungen"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_250.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, mpl=mpl, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)


    def plotalphaOutrigger(self, buildingProp):

        self.gui.label_264.setText(str(buildingProp.alpha_ist_rand))
        self.gui.label_266.setText(str(buildingProp.alpha_ist_eck))

        # Y-Achse:
        if self.gui.comboBox_yAxis_262.currentText() == 'Querschnitssabmessungen St??tzen in Abschnitten mit Outriggern':
            if buildingProp.tragwerk == 'Outrigger':
                y = []
                for n,i in enumerate(buildingProp.posOut_abschnitt):
                    y.append(buildingProp.t_st??tze[n])
                    y.append(buildingProp.t_randSt??tzen[i])
                    y.append(buildingProp.t_eckSt??tzen[i])
            else:
                y = [0]

        ymax = 1.1*max(y)
        y = [y]
        ylabel = 'QS - Abmessungen [cm]'

            
        # X-Achse:
        if self.gui.comboBox_xAxis_261.currentText() == 'Soll - und Ist - Werte der Elemente':
            if buildingProp.tragwerk == 'Outrigger':
                xticklabels = ['']
                for i in range (0, len(buildingProp.posOut)):
                    j = i + 1
                    xticklabels.append('t_soll_{}'.format(j))
                    xticklabels.append('t_ist_Rand_{}'.format(j))
                    xticklabels.append('t_ist_Eck_{}'.format(j))
            else: 
                xticklabels = [0]


      
        # Legende:
        xlabel = None
        annotations = 'individual'
        annotations_position = 'center'
        title = 'Vergleich der QS- Abmessung nach Steifigkeitsverh??ltnis alpha (Soll) und der endg??ltigen Abmessung (Ist)'
        legend = ['QS - Abmessung']

        breiteMax = len(buildingProp.posOut)*3 - 0.5
        xlim=[-0.5, breiteMax]
        ylim=[0,ymax]

        # Speichern:
        dir_fileName = "alpha_Outrigger Vergleich"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_261.plotBarChart(y, xlabel=xlabel, ylabel=ylabel, title=title, legend=legend, xticklabels=xticklabels,
                                       annotations=annotations, annotations_position=annotations_position, xlim=xlim, ylim=ylim, dir_fileName=dir_fileName, vLines=None, savePlt=savePlt, saveTex=saveTex)

    def plotAusnutzung(self, buildingProp, materialProp):

        y=[]
        if self.gui.comboBox_yAxis_232.currentText() == 'Geb??udeh??he':
            if self.gui.comboBox_xAxis_271.currentText() == 'Belt':
                for i in buildingProp.H_outrigger:
                    y.append(i)
                    y.append(i-buildingProp.h_geschoss)
            else:
                #y.append(buildingProp.h_total)
                for i in range (0,len(buildingProp.I)):
                    a=max(buildingProp.h_total-i*buildingProp.n_abschnitt*buildingProp.h_geschoss,0)
                    y.append(a)
                    a=max(buildingProp.h_total-(i+1)*buildingProp.n_abschnitt*buildingProp.h_geschoss,0)
                    y.append(a)
            ylabel = 'H??he [m]'

            
        # X-Achse:
        if self.gui.comboBox_xAxis_271.currentText() == 'Kern':
            x = []
            for i in buildingProp.sigma_kern:
                x.append(i/materialProp.f)
                x.append(i/materialProp.f)
            
            xlabel = 'Ausnutzungsgrad Kern'

        if self.gui.comboBox_xAxis_271.currentText() == 'Innenst??tze':
            x = []
            for i in buildingProp.sigma_innenSt??tze:
                x.append(i/materialProp.f)
                x.append(i/materialProp.f)
            
            xlabel = 'Ausnutzungsgrad Innenst??tze'

        if self.gui.comboBox_xAxis_271.currentText() == 'Eckst??tze':
            x = []
            for i in buildingProp.sigma_eckSt??tze:
                x.append(i/materialProp.f)
                x.append(i/materialProp.f)
            
            xlabel = 'Ausnutzungsgrad Eckst??tze'
            
        if self.gui.comboBox_xAxis_271.currentText() == 'Randst??tze':
            x = []
            for i in buildingProp.sigma_randSt??tze:
                x.append(i/materialProp.f)
                x.append(i/materialProp.f)
            
            xlabel = 'Ausnutzungsgrad Randst??tze'


        if self.gui.comboBox_xAxis_271.currentText() == 'Belt':
            x = []
            for i in buildingProp.sigma_belt:
                x.append(i/materialProp.f)
                x.append(i/materialProp.f)
            
            xlabel = 'Ausnutzungsgrad Belt'

        title = 'Ausnutzungsgrad ??ber die H??he'

        title=None
        #mpl='plotStyle_legendeRechts'
        maxWert=1.0
        minWert=0.8*min(x)
        bereich=maxWert-minWert
        xlim=[minWert-0.1*bereich,maxWert+0.1*bereich]

        ylim = [0, buildingProp.h_total ]

        legend = ['Ausnutzungsgrad nach GzT (Achtung ohne GzG)']

        # Speichern als pdf_tex und als png:
        # Change to current file location
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Querschnitssausnutzung"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_271.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, xlim=xlim, ylim = ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)

    ####    Plots Multiberechnung:

    def plotResourceAnalysis(self,buildingProp,materialProp):

        # Y-Achse:
        y4=[]
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_parameter.currentText()!='Keiner':
            if self.gui.comboBox_yAxis_312.currentText() == 'Masse':
                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append(buildingProp.multi_G_decken_opt[-(i+1)]/1000)
                    y3.append(buildingProp.multi_G_innenSt??tzen_opt[-(i+1)]/1000)
                    y2.append(buildingProp.multi_G_au??enSt??tzen_opt[-(i+1)]/1000)
                    y1.append(buildingProp.multi_G_aussteifung_opt[-(i+1)]/1000)
            
                ylabel = 'Masse [kt]'
            
            if self.gui.comboBox_yAxis_312.currentText() == 'Treibhauspotenzial GWP':

                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append((buildingProp.multi_G_decken_opt[-(i+1)]*materialProp.GWP_decke)/1000)
                    y3.append((buildingProp.multi_G_innenSt??tzen_opt[-(i+1)]*materialProp.GWP_tragwerk)/1000)
                    y2.append((buildingProp.multi_G_au??enSt??tzen_opt[-(i+1)]*materialProp.GWP_tragwerk)/1000)
                    y1.append((buildingProp.multi_G_aussteifung_opt[-(i+1)]*materialProp.GWP_tragwerk)/1000)
            
                ylabel = 'Treibhauspotenzial [t CO2-??qui.]'

            if self.gui.comboBox_yAxis_312.currentText() == 'Prim??renergiebedarf PET':

                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append((buildingProp.multi_G_decken_opt[-(i+1)]*materialProp.PET_decke)/1000)
                    y3.append((buildingProp.multi_G_innenSt??tzen_opt[-(i+1)]*materialProp.PET_tragwerk)/1000)
                    y2.append((buildingProp.multi_G_au??enSt??tzen_opt[-(i+1)]*materialProp.PET_tragwerk)/1000)
                    y1.append((buildingProp.multi_G_aussteifung_opt[-(i+1)]*materialProp.PET_tragwerk)/1000)

                ylabel = 'Prim??renergiebedarf PET [GJ]'

        else:
            if self.gui.comboBox_yAxis_312.currentText() == 'Masse':
                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append(buildingProp.multi_G_decken[-(i+1)]/1000)
                    y3.append(buildingProp.multi_G_innenSt??tzen[-(i+1)]/1000)
                    y2.append(buildingProp.multi_G_au??enSt??tzen[-(i+1)]/1000)
                    y1.append(buildingProp.multi_G_aussteifung[-(i+1)]/1000)
            
                ylabel = 'Masse [kt]'
            
            if self.gui.comboBox_yAxis_312.currentText() == 'Treibhauspotenzial GWP':

                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append((buildingProp.multi_G_decken[-(i+1)]*materialProp.GWP_decke)/1000)
                    y3.append((buildingProp.multi_G_innenSt??tzen[-(i+1)]*materialProp.GWP_tragwerk)/1000)
                    y2.append((buildingProp.multi_G_au??enSt??tzen[-(i+1)]*materialProp.GWP_tragwerk)/1000)
                    y1.append((buildingProp.multi_G_aussteifung[-(i+1)]*materialProp.GWP_tragwerk)/1000)
            
                ylabel = 'Treibhauspotenzial [t CO2-??qui.]'

            if self.gui.comboBox_yAxis_312.currentText() == 'Prim??renergiebedarf PET':

                for i in range(0,len(buildingProp.multi_G_decken)):
                    y4.append((buildingProp.multi_G_decken[-(i+1)]*materialProp.PET_decke)/1000)
                    y3.append((buildingProp.multi_G_innenSt??tzen[-(i+1)]*materialProp.PET_tragwerk)/1000)
                    y2.append((buildingProp.multi_G_au??enSt??tzen[-(i+1)]*materialProp.PET_tragwerk)/1000)
                    y1.append((buildingProp.multi_G_aussteifung[-(i+1)]*materialProp.PET_tragwerk)/1000)

                ylabel = 'Prim??renergiebedarf PET [GJ]'


        y=[y1,y2,y3,y4]

        # X-Achse:
        x1=[]
        if self.gui.comboBox_xAxis_311.currentText() == 'Geb??udeh??he':
            for i in range (0,len(buildingProp.multi_h_total)):
                x1.append(str(buildingProp.multi_h_total[-(i+1)]))

            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_311.currentText() == 'Schlankheit':
            for i in range (0,len(buildingProp.multi_h_total)):
                x1.append(str(buildingProp.multi_schlankheit[-(i+1)]))

            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_311.currentText() == 'Anzahl Stockwerke':
            for i in range (0,len(buildingProp.multi_h_total)):
                x1.append(str(buildingProp.multi_n[-(i+1)]))

            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        xticks = numpy.arange(len(y[0]))
        xticklabels=x1
        annotations=None#'individual'
        annotations_position='center'

        xlim=[-0.5,(len(y[0])-0.5)]
        
        l1='Aussteifung'
        if buildingProp.tragwerk == 'Kerntragwerk':
            l1='Kern'
        if buildingProp.tragwerk == 'Braced Tube':
            l1='Diagonalen + Querstreben'
        if buildingProp.tragwerk == 'Rahmentragwerk':
            l1='Riegel'
        
        #legend=['Decken','Innenst??tzen','Au??enst??tzen', l1]
        legend=[l1,'Au??enst??tzen','Innenst??tzen','Decken',]
        title =None #'Resourcenverbrauch nach Bauteilen'

        mpl='plotStyle_barChart_small'
       
        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Bauteilmengenverteilung nach H??he"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_310.plotBarChart(y, xlabel=xlabel, ylabel=ylabel, title=title, xticks=xticks, xticklabels=xticklabels, legend=legend, barChart='stacked', mpl=mpl,
                annotations=annotations, annotations_position=annotations_position, xlim=xlim, dir_fileName=dir_fileName, vLines=None, savePlt=savePlt, saveTex=saveTex)
        

    def plotResourceAnalysis_m2(self,buildingProp,materialProp):

        # Y-Achse:
        y4=[]
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_322.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multi_G_decken)):
                y4.append(buildingProp.multi_G_decken[i]/buildingProp.multi_A[i])
                y3.append(buildingProp.multi_G_innenSt??tzen[i]/buildingProp.multi_A[i])
                y2.append(buildingProp.multi_G_au??enSt??tzen[i]/buildingProp.multi_A[i])
                y1.append(buildingProp.multi_G_aussteifung[i]/buildingProp.multi_A[i])
            
            ylabel = 'Masse [t/m??]'
            
        if self.gui.comboBox_yAxis_322.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y4.append(buildingProp.multi_G_decken[i]*materialProp.GWP_decke/buildingProp.multi_A[i])
                y3.append(buildingProp.multi_G_innenSt??tzen[i]*materialProp.GWP_tragwerk/buildingProp.multi_A[i])
                y2.append(buildingProp.multi_G_au??enSt??tzen[i]*materialProp.GWP_tragwerk/buildingProp.multi_A[i])
                y1.append(buildingProp.multi_G_aussteifung[i]*materialProp.GWP_tragwerk/buildingProp.multi_A[i])
            
            ylabel = 'Treibhauspotenzial [kg CO2-??qui./m??]'

        if self.gui.comboBox_yAxis_322.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y4.append(buildingProp.multi_G_decken[i]*materialProp.PET_decke/buildingProp.multi_A[i])
                y3.append(buildingProp.multi_G_innenSt??tzen[i]*materialProp.PET_tragwerk/buildingProp.multi_A[i])
                y2.append(buildingProp.multi_G_au??enSt??tzen[i]*materialProp.PET_tragwerk/buildingProp.multi_A[i])
                y1.append(buildingProp.multi_G_aussteifung[i]*materialProp.PET_tragwerk/buildingProp.multi_A[i])
            
            ylabel = 'Prim??renergiebedarf PET [MJ/m??]'

        # X-Achse:
        if self.gui.comboBox_xAxis_321.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_321.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_321.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'

        # Legende:
        l1='Aussteifung'
        if buildingProp.tragwerk == 'Kerntragwerk':
            l1='Kern'
        if buildingProp.tragwerk == 'Braced Tube':
            l1='Diagonalen'
        if buildingProp.tragwerk == 'Rahmentragwerk':
            l1='Riegel'

        legend = [l1,'Au??enst??tzen','Innenst??tzen','Decken']
        title =None     #'Resourcenverbrauch nach Bauteilen pro Geschossfl??che'
       
        x=[x1,x1,x1,x1]
        y=[y1,y2,y3,y4]

        maxWert=max(max(y1),max(y2),max(y3),max(y4))
        minWert=min(min(y1),min(y2),min(y3),min(y4))
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert]

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Bauteilmengen nach H??he pro m2"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_320.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)


    def plotPFHAnalysis(self,buildingProp,materialProp):

        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_332.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append(buildingProp.multi_G_decken[i])
                y2.append(buildingProp.multi_G_decken[i]+buildingProp.multi_G_totalOhnePFH[i])
                y1.append(buildingProp.multi_G_decken[i]+buildingProp.multi_G_total[i])
            
            ylabel = 'Masse [t]'
            
        if self.gui.comboBox_yAxis_332.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke)/1000)
                y2.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke+buildingProp.multi_G_totalOhnePFH[i]*materialProp.GWP_tragwerk)/1000)
                y1.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke+buildingProp.multi_G_total[i]*materialProp.GWP_tragwerk)/1000)
            
            ylabel = 'Treibhauspotenzial [t CO2-??qui.]'

        if self.gui.comboBox_yAxis_332.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke)/1000)
                y2.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke+buildingProp.multi_G_totalOhnePFH[i]*materialProp.PET_tragwerk)/1000)
                y1.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke+buildingProp.multi_G_total[i]*materialProp.PET_tragwerk)/1000)
            
            ylabel = 'Prim??renergiebedarf PET [GJ]'

        # X-Achse:
        if self.gui.comboBox_xAxis_331.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_331.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_331.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
         
        # Legende:
        legend = ['Tragwerk inkl. Horizontallastabtrag','Tragwerk ohne Horizontallastabtrag','Decken']
        title =None     #'Kummulierter Resourcenverbrauch nach Lastabtrag'
       
        x=[x1,x1,x1]
        y=[y1,y2,y3]
        
        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Ressourcen nach Lastabtrag und H??he"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_330.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)


    def plotPFHAnalysis_m2(self,buildingProp,materialProp):

        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_342.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append(buildingProp.multi_G_decken[i]/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken[i]+buildingProp.multi_G_totalOhnePFH[i])/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken[i]+buildingProp.multi_G_total[i])/buildingProp.multi_A[i])
                #y4.append(1)
            
            ylabel = 'Masse [t/m??]'
            
        if self.gui.comboBox_yAxis_342.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append(buildingProp.multi_G_decken[i]*materialProp.GWP_decke/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke+buildingProp.multi_G_totalOhnePFH[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke+buildingProp.multi_G_total[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Treibhauspotenzial [kg CO2-??qui./m??]'

        if self.gui.comboBox_yAxis_342.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multi_G_decken)):
                y3.append(buildingProp.multi_G_decken[i]*materialProp.PET_decke/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke+buildingProp.multi_G_totalOhnePFH[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke+buildingProp.multi_G_total[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Prim??renergiebedarf PET [MJ/m??]'

        # X-Achse:
        if self.gui.comboBox_xAxis_341.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_341.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_341.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        legend = ['Tragwerk inkl. Horizontallastabtrag','Tragwerk ohne Horizontallastabtrag','Decken']
        title =None     #'Kummulierter Resourcenverbrauch nach Lastabtrag pro Geschossfl??che'
        
        maxWert=max(max(y1),max(y2),max(y3))
        minWert=min(min(y1),min(y2),min(y3))
        bereich=maxWert-minWert
        ylim=[max(0,minWert-0.1*bereich),maxWert]
        #ylim=[0.3,5.3]

        mpl='plotStyle_plot2D_small'

        if buildingProp.i_wechselNW == 'none' or buildingProp.i_wechselNW == 0:
            vLines=None
            vTexts=None

        else:
            vLines = [(x1[buildingProp.i_wechselNW]+x1[buildingProp.i_wechselNW-1])/2]
            vTexts = ['Wechsel ma??g. NW']

       
        x=[x1,x1,x1]
        y=[y1,y2,y3]

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Ressourcen nach Lastabtrag und H??he pro m2"
        saveTex = False
        savePlt = True
        savePkl = False
        
        self.gui.MplWidget_340.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex, savePkl=savePkl)


    def plotEigenFrequency(self,buildingProp):

        # Y-Achse:       
        if self.gui.comboBox_yAxis_352.currentText() == 'Eigenfrequenz f':
            y=buildingProp.multi_eigenFrequenz
                       
            ylabel = 'f [Hz]'
            legend=['Eigenfrequenz']
            
        if self.gui.comboBox_yAxis_352.currentText() == 'Eigenkreisfrequenz Omega':
            y=[element*2*numpy.pi for element in buildingProp.multi_eigenFrequenz]
                       
            ylabel = 'Omega [rad/s]'
            legend=['Eigenkreisfrequenz']

        if self.gui.comboBox_yAxis_352.currentText() == 'Periodendauer T':
            y=[1/element for element in buildingProp.multi_eigenFrequenz]
                       
            ylabel = 'T [s]'
            legend=['Periodendauer']

        # X-Achse:
        if self.gui.comboBox_xAxis_351.currentText() == 'Geb??udeh??he':
            x=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_351.currentText() == 'Schlankheit':
            x=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_351.currentText() == 'Anzahl Stockwerke':
            x=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        title = None #'Kerntragwerk'       #'Eigenwerte nach Geb??udeh??he'
        ylim=[]

        maxWert=max(y)
        minWert=min(y)
        bereich=maxWert-minWert
        ylim=[max(minWert-0.1*bereich,0),maxWert]
        #ylim=[.05,0.45]

        mpl='plotStyle_plot2D_small'

        #mpl='plotStyle_legendeRechts'

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Eigenfrequenz nach H??he"  
        saveTex = False
        savePlt = True
        savePkl = False
        
        self.gui.MplWidget_350.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex,  savePkl=savePkl)
    
    def plotMoment(self,buildingProp):

        # Y-Achse:       
        if self.gui.comboBox_yAxis_362.currentText() == 'Maximales Biegemoment':
            y=buildingProp.multi_moment
                       
            ylabel = 'M [MNm]'
            legend=['Maximales Biegemoment']
            

        # X-Achse:
        if self.gui.comboBox_xAxis_361.currentText() == 'Geb??udeh??he':
            x=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_361.currentText() == 'Schlankheit':
            x=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_361.currentText() == 'Anzahl Stockwerke':
            x=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        title = None #'Eigenwerte nach Geb??udeh??he'
        ylim=[]

        maxWert=max(y)
        minWert=min(y)
        bereich=maxWert-minWert
        ylim=[max(minWert-0.1*bereich,0),maxWert]
        #ylim=[.05,0.45]

        mpl='plotStyle_plot2D_small'

        #mpl='plotStyle_legendeRechts'

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName =  "Biegemoment nach H??he"
        saveTex = False
        savePlt = True
        savePkl = False
        
        self.gui.MplWidget_360.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex,  savePkl=savePkl)
    

    ### Parameteranalyse:

    def checkIfPosOutVariing(self, buildingProp):
        if self.gui.comboBox_parameter.currentText() == 'Outriggeranordnung ??ber die H??he':
            buildingProp.varPosOut = True
        else: 
            buildingProp.varPosOut = False

    def checkIfPosOutVariing_2D(self, buildingProp):
        if self.gui.comboBox_parameter_2.currentText() == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
            buildingProp.varPosOut_2D = True
        else: 
            buildingProp.varPosOut_2D = False
    
            
    def plotDeformation(self,buildingProp,materialProp):
        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]

        for i in range(0,len(buildingProp.multiPar_w)):
                y3.append(buildingProp.multiPar_w_GA[i]*1000)
                y2.append(buildingProp.multiPar_w_EI[i]*1000)
                y1.append(buildingProp.multiPar_w[i]*1000)
            
        ylabel = 'Verformung [mm]'
        # X-Achse:
        x1=buildingProp.multi_p

        if buildingProp.parameter=='Staffelung n_abschnitt':
            xlabel='n_{abschnitt}'
        if buildingProp.parameter=='Deckengewicht g_k1':
            xlabel='g_{k1}'
        if buildingProp.parameter=='Windlast q_b':
            xlabel='q_b'
        if buildingProp.parameter=='Geb??udeform psi_r':
            xlabel='\psi_r'
        if buildingProp.parameter == 'Anzahl Stiele pro Raster (Framed Tube) n_stiele/b_raster':
            xlabel='n_{stiele}/b_{raster}'        
        if buildingProp.parameter == 'Verformungsverh??ltnis w_EI/w_GA':
            xlabel='w_{EI}/w_{GA}'
        if buildingProp.parameter == 'Mindestwirksamkeit Framed Tube tube_wirksam':
            xlabel='Zielwert tube_{wirksam}'
        if buildingProp.parameter == 'Anzahl Outrigger':
            xlabel = 'n_{Outrigger}'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
            xlabel = 'alpha_outrigger'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Beta (Outrigger)':
            xlabel = 'beta_outrigger'

        if buildingProp.parameter == 'Outriggeranordnung ??ber die H??he':
            xlabel = 'Outriggerposition: xi*H'

        # Legende:
        legend = ['w','w_EI', 'w_GA']

        title =None

        mpl='plotStyle_plot2D'

        maxWert=max(max(y1),max(y2),max(y3))
        minWert=min(min(y1),min(y2),min(y3))
        x=[x1,x1,x1]
        y=[y1,y2,y3]
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert+0.1*bereich]
        ylim=[min(y1),maxWert]

        vLines=None
        vTexts=None
        

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Verformung nach Parameter"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_412.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                    legend=legend, ylim=ylim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)

    def plotAlphaIst(self,buildingProp,materialProp):

        y2=[]
        y1=[]

        if self.gui.comboBox_yAxis_517.currentText() == 'alpha_ist Outrigger 1':
            position = 0
        elif self.gui.comboBox_yAxis_517.currentText() == 'alpha_ist Outrigger 2':
            position = 1
        elif self.gui.comboBox_yAxis_517.currentText() == 'alpha_ist Outrigger 3':
            position = 2
        elif self.gui.comboBox_yAxis_517.currentText() == 'alpha_ist Outrigger 4':
            position = 3

        if position < buildingProp.n_outrigger:
            for i in range(0,len(buildingProp.multi_p)):
                        y2.append(buildingProp.multiPar_alpha_ist_eck[i][position])
                        y1.append(buildingProp.multiPar_alpha_ist_rand[i][position])

            x1=buildingProp.multi_p
        
        else: 
            x1 = []
        
        ylabel = 'Wirkliches Steifigkeitsverh??ltnis zwischen Kern und St??tzen'

    
        xlabel = 'alpha_outrigger'

        # Legende:
        legend = ['alpha_ist Randst??tze','alpha_ist Eckst??tze']

        title =None

        mpl='plotStyle_plot2D'

        maxWert=max(max(y1),max(y2))
        minWert=min(min(y1),min(y2))
        x=[x1,x1]
        y=[y1,y2]
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert+0.1*bereich]
        #ylim=[min(y1),maxWert[0]]

        vLines=None
        vTexts=None
        

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Alpha_ist der Outriggerst??tzen"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_512.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                    legend=legend, ylim=ylim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)



    def plotParameterAnalysis(self,buildingProp,materialProp):

        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_412.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multiPar_G_decken)):
                y3.append(buildingProp.multiPar_G_decken[i])
                y2.append(buildingProp.multiPar_G_decken[i]+buildingProp.multiPar_G_totalOhnePFH[i])
                y1.append(buildingProp.multiPar_G_decken[i]+buildingProp.multiPar_G_total[i])
            
            ylabel = 'Masse [t]'
            
        if self.gui.comboBox_yAxis_412.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multiPar_G_decken)):
                y3.append((buildingProp.multiPar_G_decken[i]*materialProp.GWP_decke)/1000)
                y2.append((buildingProp.multiPar_G_decken[i]*materialProp.GWP_decke+buildingProp.multiPar_G_totalOhnePFH[i]*materialProp.GWP_tragwerk)/1000)
                y1.append((buildingProp.multiPar_G_decken[i]*materialProp.GWP_decke+buildingProp.multiPar_G_total[i]*materialProp.GWP_tragwerk)/1000)
            
            ylabel = 'Treibhauspotenzial [t CO2-??qui.]'

        if self.gui.comboBox_yAxis_412.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multiPar_G_decken)):
                y3.append((buildingProp.multiPar_G_decken[i]*materialProp.PET_decke)/1000)
                y2.append((buildingProp.multiPar_G_decken[i]*materialProp.PET_decke+buildingProp.multiPar_G_totalOhnePFH[i]*materialProp.PET_tragwerk)/1000)
                y1.append((buildingProp.multiPar_G_decken[i]*materialProp.PET_decke+buildingProp.multiPar_G_total[i]*materialProp.PET_tragwerk)/1000)
            
            ylabel = 'Prim??renergiebedarf PET [GJ]'

        if self.gui.comboBox_yAxis_412.currentText() == 'Verformung':

            for i in range(0,len(buildingProp.multiPar_w)):
                y3.append(buildingProp.multiPar_w_GA[i]*1000)
                y2.append(buildingProp.multiPar_w_EI[i]*1000)
                y1.append(buildingProp.multiPar_w[i]*1000)
            
            ylabel = 'Verformung [mm]'

        x1 = buildingProp.multi_p

        if buildingProp.parameter=='Staffelung n_abschnitt':
            xlabel='n_{abschnitt}'
        if buildingProp.parameter=='Deckengewicht g_k1':
            xlabel='g_{k1}'
        if buildingProp.parameter=='Windlast q_b':
            xlabel='q_b'
        if buildingProp.parameter=='Geb??udeform psi_r':
            xlabel='\psi_r'
        if buildingProp.parameter == 'Anzahl Stiele pro Raster (Framed Tube) n_stiele/b_raster':
            xlabel='n_{stiele}/b_{raster}'        
        if buildingProp.parameter == 'Verformungsverh??ltnis w_EI/w_GA':
            xlabel='w_{EI}/w_{GA}'
        if buildingProp.parameter == 'Mindestwirksamkeit Framed Tube tube_wirksam':
            xlabel='Zielwert tube_{wirksam}'
        if buildingProp.parameter == 'Anzahl Outrigger':
            xlabel = 'n_{Outrigger}'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
            xlabel = 'alpha'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Beta (Outrigger)':
            xlabel = 'beta'
        if buildingProp.parameter == 'Outriggeranordnung ??ber die H??he':
            xlabel = 'Outriggerposition: xi*H'
        
        
        # Legende:
        if self.gui.comboBox_yAxis_412.currentText() == 'Verformung':
            legend = ['w','w_EI', 'w_GA']
        else:
            legend = ['Tragwerk inkl. Horizontallastabtrag']#,'Tragwerk ohne Horizontallastabtrag','Decken']
        title =None     #'Parametereinfluss auf Kummulierten Resourcenverbrauch'
        
        mpl='plotStyle_plot2D'
        
        
        maxWert=max(max(y1),max(y2),max(y3))
        minWert=min(min(y1),min(y2),min(y3))
        x=[x1]#,x1,x1]
        y=[y1]#,y2,y3]
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert+0.1*bereich]
        ylim=[min(y1),maxWert]
            
        vLines=None
        vTexts=None
        

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Materialmenge nach Parameter"
        saveTex = False
        savePlt = True
        savePkl = False
        
        self.gui.MplWidget_410.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                    legend=legend, ylim=ylim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex, savePkl=savePkl)


    def plotOptimalParameter(self,buildingProp,materialProp):

        # Y-Achse:
        y1=buildingProp.multi_p_opt
        y2=buildingProp.multi_p_opt_min

        y=[y1,y2]
        
        if buildingProp.parameter=='Staffelung n_abschnitt':
            ylabel='n_{abschnitt}'
        if buildingProp.parameter=='Deckengewicht g_k1':
            ylabel='g_{k1}'
        if buildingProp.parameter=='Windlast q_b':
            ylabel='q_b'
        if buildingProp.parameter=='Geb??udeform psi_r':
            ylabel='\psi_r'
        if buildingProp.parameter == 'Anzahl Stiele pro Raster (Framed Tube) n_stiele/b_raster':
            ylabel='n_{stiele}/b_{raster}'        
        if buildingProp.parameter == 'Verformungsverh??ltnis w_EI/w_GA':
            ylabel='w_{EI}/w_{GA}'
        if buildingProp.parameter == 'Mindestwirksamkeit Framed Tube tube_wirksam':
            ylabel='Zielwert tube_{wirksam}'
        if buildingProp.parameter == 'Anzahl Outrigger':
            ylabel = 'n_{Outrigger}'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
            ylabel = 'alpha'
        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Beta (Outrigger)':
            ylabel = 'beta'
        if buildingProp.parameter == 'Outriggeranordnung ??ber die H??he':
            ylabel = 'Outriggerposition: xi*H'


        
        # X-Achse:
        if self.gui.comboBox_xAxis_421.currentText() == 'Geb??udeh??he':
            x=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_421.currentText() == 'Schlankheit':
            x=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_421.currentText() == 'Anzahl Stockwerke':
            x=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        x=[x,x]

        # Legende:
        legend = ['Obergrenze Optimum','Untergrenze Optimum']
        title =None     #'Parametereinfluss auf Kummulierten Resourcenverbrauch'
        
        maxWert=max(y1)
        minWert=min(y2)
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert+0.1*bereich]

        mpl='plotStyle_plot2D_small'

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Optimaler Parameter nach H??he"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_420.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)
    
    
    def plotPFHAnalysisOptimalParameter_m2(self,buildingProp,materialProp):

        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_432.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append(buildingProp.multi_G_decken_opt[i]/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken_opt[i]+buildingProp.multi_G_totalOhnePFH_opt[i])/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]+buildingProp.multi_G_total_opt[i])/buildingProp.multi_A[i])
            
            ylabel = 'Masse [t/m??]'
            
        if self.gui.comboBox_yAxis_432.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append(buildingProp.multi_G_decken_opt[i]*materialProp.GWP_decke/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken_opt[i]*materialProp.GWP_decke+buildingProp.multi_G_totalOhnePFH_opt[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]*materialProp.GWP_decke+buildingProp.multi_G_total_opt[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Treibhauspotenzial [kg CO2-??qui./m??]'

        if self.gui.comboBox_yAxis_432.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append(buildingProp.multi_G_decken_opt[i]*materialProp.PET_decke/buildingProp.multi_A[i])
                y2.append((buildingProp.multi_G_decken_opt[i]*materialProp.PET_decke+buildingProp.multi_G_totalOhnePFH_opt[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]*materialProp.PET_decke+buildingProp.multi_G_total_opt[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Prim??renergiebedarf PET [MJ/m??]'

        # X-Achse:
        if self.gui.comboBox_xAxis_431.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_431.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_431.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        legend = ['Tragwerk inkl. Horizontallastabtrag','Tragwerk ohne Horizontallastabtrag','Decken']
        title =None     #'Optimierter Kummulierter Resourcenverbrauch nach Lastabtrag pro Geschossfl??che'
        
        maxWert=max(max(y1),max(y2),max(y3))
        minWert=min(min(y1),min(y2),min(y3))
        bereich=maxWert-minWert
        ylim=[minWert-0.1*bereich,maxWert]
        #ylim=[0.60,2.05]
        mpl='plotStyle_plot2D_small'
       
        x=[x1,x1,x1]
        y=[y1,y2,y3]

        if buildingProp.i_wechselNW == 'none' or buildingProp.i_wechselNW ==0:
            vLines=None
            vTexts=None

        else:
            vLines = [(x1[buildingProp.i_wechselNW]+x1[buildingProp.i_wechselNW-1])/2]
            vTexts = ['Wechsel ma??g. NW']

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName ="Resourcen nach Lastabtrag und H??he mit optimalem Parameter pro m2"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_430.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)
    

    def plotPFHAnalysisParameterComparisson_m2(self,buildingProp,materialProp):

        # Y-Achse:
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_442.currentText() == 'Masse':
            
            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y2.append((buildingProp.multi_G_decken[i]+buildingProp.multi_G_total[i])/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]+buildingProp.multi_G_total_opt[i])/buildingProp.multi_A[i])
            
            ylabel = 'Masse [t/m??]'
            
        if self.gui.comboBox_yAxis_442.currentText() == 'Treibhauspotenzial GWP':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y2.append((buildingProp.multi_G_decken[i]*materialProp.GWP_decke+buildingProp.multi_G_total[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]*materialProp.GWP_decke+buildingProp.multi_G_total_opt[i]*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Treibhauspotenzial [kg CO2-??qui./m??]'

        if self.gui.comboBox_yAxis_442.currentText() == 'Prim??renergiebedarf PET':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y2.append((buildingProp.multi_G_decken[i]*materialProp.PET_decke+buildingProp.multi_G_total[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
                y1.append((buildingProp.multi_G_decken_opt[i]*materialProp.PET_decke+buildingProp.multi_G_total_opt[i]*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Prim??renergiebedarf PET [MJ/m??]'

        # X-Achse:
        if self.gui.comboBox_xAxis_441.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_441.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_441.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        legend = ['Optimierter Verlauf','Verlauf mit festem Parameter']
        title =None     #'Optimierter Gesamtresourcenverbrauch pro Geschossfl??che'
        
        maxWert=max(max(y1),max(y2))
        minWert=min(min(y1),min(y2))
        bereich=maxWert-minWert
        ylim=[max(0,minWert-0.1*bereich),maxWert]

        mpl='plotStyle_plot2D'

        x=[x1,x1]
        y=[y1,y2]

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Vergleich Resourcenverbrauch nach H??he mit und ohne optimalem Parameter"
        saveTex = False
        savePlt = True
        savePkl = False
        
        self.gui.MplWidget_440.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title, mpl=mpl,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex, savePkl=savePkl)


    def plotInfluenceOptimalParameter_m2(self,buildingProp,materialProp):

        # Y-Achse:
        y3=[]
        y2=[]
        y1=[]
        
        if self.gui.comboBox_yAxis_452.currentText() == 'Massendifferenz':
            
            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])/buildingProp.multi_A[i])
                y2.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])+(buildingProp.multi_G_totalOhnePFH[i]-buildingProp.multi_G_totalOhnePFH_opt[i]))/buildingProp.multi_A[i])
                y1.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])+(buildingProp.multi_G_total[i]-buildingProp.multi_G_total_opt[i]))/buildingProp.multi_A[i])
            
            ylabel = 'Differenz Masse [t/m??]'
            
        if self.gui.comboBox_yAxis_452.currentText() == 'Treibhauspotenzialdifferenz':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.GWP_decke/buildingProp.multi_A[i])
                y2.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.GWP_decke+(buildingProp.multi_G_totalOhnePFH[i]-buildingProp.multi_G_totalOhnePFH_opt[i])*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
                y1.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.GWP_decke+(buildingProp.multi_G_total[i]-buildingProp.multi_G_total_opt[i])*materialProp.GWP_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Differenz Treibhauspotenzial [kg CO2-??qui./m??]'

        if self.gui.comboBox_yAxis_452.currentText() == 'Prim??renergiebedarfdifferenz':

            for i in range(0,len(buildingProp.multi_G_decken_opt)):
                y3.append((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.PET_decke/buildingProp.multi_A[i])
                y2.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.PET_decke+(buildingProp.multi_G_totalOhnePFH[i]-buildingProp.multi_G_totalOhnePFH_opt[i])*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
                y1.append(((buildingProp.multi_G_decken[i]-buildingProp.multi_G_decken_opt[i])*materialProp.PET_decke+(buildingProp.multi_G_total[i]-buildingProp.multi_G_total_opt[i])*materialProp.PET_tragwerk)/buildingProp.multi_A[i])
            
            ylabel = 'Differenz Prim??renergiebedarf PET [MJ/m??]'

        # X-Achse:
        if self.gui.comboBox_xAxis_451.currentText() == 'Geb??udeh??he':
            x1=buildingProp.multi_h_total
            xlabel = 'H??he [m]'

        if self.gui.comboBox_xAxis_451.currentText() == 'Schlankheit':
            x1=buildingProp.multi_schlankheit
            xlabel = 'Schlankheit [h/b]'

        if self.gui.comboBox_xAxis_451.currentText() == 'Anzahl Stockwerke':
            x1=buildingProp.multi_n
            xlabel = 'Anzahl Stockwerke [-]'
        
        # Legende:
        legend = ['Tragwerk inkl. Horizontallastabtrag','Tragwerk ohne Horizontallastabtrag','Decken']
        title =None     #'Einfluss der Parameteroptimierung auf Resourcenverbrauch nach Lastabtrag pro Geschossfl??che'
        
        maxWert=max(max(y1),max(y2),max(y3))
        minWert=min(min(y1),min(y2),min(y3))
        bereich=maxWert-minWert
        ylim=[max(0,minWert-0.1*bereich),maxWert]
        #ylim=[0.3,5.3]
       
        x=[x1,x1,x1]
        y=[y1,y2,y3]

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Einsparungen nach H??he bei Parameteroptimierung"
        saveTex = False
        savePlt = True
        
        self.gui.MplWidget_450.plot2D(x, y, xlabel=xlabel, ylabel=ylabel, title=title,
                                      legend=legend, ylim=ylim, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)
    
    def plot3DParameter(self,buildingProp,materialProp):

        # if self.gui.comboBox_darstellung_473.currentText() == 'Scatter':

        #     # Gr????e z array
        #     #laenge = len(buildingProp.multi_p_1)

        #     # Z-Achse:
        #     #z=np.zeros((laenge,laenge))
        #     z = []
            
        #     if self.gui.comboBox_zAxis_472.currentText() == 'Masse':

        #         for i in range(0,len(buildingProp.multiPar_2D_G_decken)):
        #             z.append(buildingProp.multiPar_2D_G_decken[i]+buildingProp.multiPar_2D_G_total[i])
                
        #         #for zeile in range(0,laenge):
        #             #for spalte in range(0,laenge):
        #                 #z[i][j]=buildingProp.multiPar_2D_G_decken[i][j]+buildingProp.multiPar_2D_G_total[i][j]
        #                 #z_i = buildingProp.multiPar_2D_G_decken[zeile][spalte]+buildingProp.multiPar_2D_G_total[zeile][spalte]
        #                 #z.append(z_i)
                
        #         zlabel = 'Gesamtmasse [t]'

        #     #if self.gui.comboBox_zAxis_472.currentText() == 'Verformung':

        #         #for i in range(0,len(buildingProp.multiPar_w)):
        #             #z.append(buildingProp.multiPar_w[i]*1000)
                
        #         #zlabel = 'Verformung [mm]'
                

        #     # Horizontale Achsen:
        #     x=buildingProp.multi_p_1
        #     y=buildingProp.multi_p_2

        #     if buildingProp.parameter_2D == 'Steifigkeitsverh??ltnis Alpha und Beta (Outrigger)':
        #         xlabel = 'alpha'
        #         ylabel = 'beta'
        #     if buildingProp.parameter_2D == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
        #         xlabel = 'Outriggerposition 1 = xi_1*H'
        #         ylabel = 'Outriggerposition 2 = xi_2*H'
            
            
        #     # Legende:
        #     #if self.gui.comboBox_zAxis_472.currentText() == 'Verformung':
        #         #legend = ['w']
        #     #else:
        #     legend = ['Tragwerk inkl. Horizontallastabtrag']
        #     title =None 
            
        #     mpl='plotStyle_plot2D'

        #     #maxWert = z[0][0]
        #     #minWert = z[0][0]

            
        #     #for i in range(0, laenge):
        #         #for j in range(0,laenge):
        #             #if z[i][j]>maxWert:
        #                 #maxWert = z[i][j]
        #             #if z[i][j]< minWert:
        #                 #minWert = z[i][j]
                        
        #     maxWert_z=max(z)
        #     minWert_z=min(z)
        #     bereich_z=maxWert_z-minWert_z
        #     zlim=[minWert_z-0.1*bereich_z,maxWert_z+0.1*bereich_z]
        #     #zlim=[minWert_z,maxWert_z]

        #     maxWert_x=max(x)
        #     minWert_x=min(x)
        #     bereich_x=maxWert_x-minWert_x
        #     xlim = [minWert_x-0.1*bereich_x,maxWert_x+0.1*bereich_x]

        #     maxWert_y=max(y)
        #     minWert_y=min(y)
        #     bereich_y=maxWert_y-minWert_y
        #     ylim = [minWert_y-0.1*bereich_y,maxWert_y+0.1*bereich_y]
            
            
        #     vLines=None
        #     vTexts=None
            

        #     # Speichern:
        #     os.chdir(os.path.dirname(sys.argv[0]))
        #     dir_fileName = "Materialmenge nach Parameter"
        #     saveTex = False
        #     savePlt = True
            
        #     #self.gui.Widget_470.plot3DScatter(x, y, z, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, title=title, mpl=mpl,
        #                                 #legend=legend, xlim = xlim, ylim = ylim, zlim=zlim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)
        
        #if self.gui.comboBox_darstellung_473.currentText() == 'Oberfl??che':

        # Gr????e z array
        laenge = int(np.sqrt(len(buildingProp.multi_p_1)))

        # Z-Achse:
        z=np.zeros((laenge,laenge))
        
        
        #if self.gui.comboBox_zAxis_472.currentText() == 'Masse':
        j=0                
        for zeile in range(0,laenge):
            for spalte in range(0,laenge):
                z[zeile][spalte] = buildingProp.multiPar_2D_G_decken[j]+buildingProp.multiPar_2D_G_total[j]
                j += 1
        
        zlabel = 'Gesamtmasse [t]'

        # Minimum bestimmen
        minWert = buildingProp.multiPar_2D_G_decken[0]+buildingProp.multiPar_2D_G_total[0]
        x_min = buildingProp.multi_p_2[0]
        y_min = buildingProp.multi_p_1[0]
        for k in range(1, len(buildingProp.multiPar_2D_G_total)):
            masse = buildingProp.multiPar_2D_G_decken[k]+buildingProp.multiPar_2D_G_total[k]
            if  masse < minWert:
                minWert = masse
                x_min = buildingProp.multi_p_2[k]
                y_min = buildingProp.multi_p_1[k]
        TP = [round(x_min,2), round(y_min,2), round(minWert,2)]
        # ??bergebe Koordinaten
        self.gui.label_476.setText(str(TP))
            

        # Horizontale Achsen:
        x = []
        y = []
        for i in range(0,laenge):
            y.append(buildingProp.multi_p_1[laenge*i])
            x.append(buildingProp.multi_p_2[i])


        if buildingProp.parameter_2D == 'Steifigkeitsverh??ltnis Alpha und Beta (Outrigger)':
            ylabel = 'alpha'
            xlabel = 'beta'
        if buildingProp.parameter_2D == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
            xlabel = 'Outriggerposition 1 = xi_1*H'
            ylabel = 'Outriggerposition 2 = xi_2*H'
        
        
        legend = None#['Tragwerk inkl. Horizontallastabtrag']
        title =None 
        
        mpl='_3D'

        maxWert_z = z[0][0]
        minWert_z = z[0][0]

        
        for i in range(0, laenge):
            for j in range(0,laenge):
                if z[i][j]>maxWert_z:
                    maxWert_z = z[i][j]
                if z[i][j]< minWert_z:
                    minWert_z = z[i][j]
                    

        bereich_z=maxWert_z-minWert_z
        zlim=[minWert_z-0.1*bereich_z,maxWert_z+0.1*bereich_z]
        #zlim=[minWert_z,maxWert_z]

        maxWert_x=max(x)
        minWert_x=min(x)
        bereich_x=maxWert_x-minWert_x
        xlim = [minWert_x-0.1*bereich_x,maxWert_x+0.1*bereich_x]

        maxWert_y=max(y)
        minWert_y=min(y)
        bereich_y=maxWert_y-minWert_y
        ylim = [minWert_y-0.1*bereich_y,maxWert_y+0.1*bereich_y]
        
        
        vLines=None
        vTexts=None
        

        # Speichern:
        os.chdir(os.path.dirname(sys.argv[0]))
        dir_fileName = "Materialmenge nach Parameter 3D"
        saveTex = False
        savePlt = True
        
        self.gui.Widget_470.plot3DSurface(x, y, z, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, title=title, mpl=mpl,
                                        legend=legend, xlim = xlim, ylim = ylim, zlim=zlim, vLines=vLines, vTexts=vTexts, dir_fileName=dir_fileName, savePlt=savePlt, saveTex=saveTex)


    #------------------------------------------------------------------------------------------
    # Berechnungen:

    def mainCalculation(self,buildingProp,loads,materialProp,DataProp):
        # Ruft die dem Tragwerk entsprechende Datei auf

        # Berechnen der Schnittgr????en innerhalb der Klasse loads - Momente sind danch Teil der Klasse
        building.loads.calculateInternalForces(loads, buildingProp)

        if buildingProp.tragwerk == 'Kerntragwerk':
            str_core.design(buildingProp, loads, materialProp,DataProp)

        if buildingProp.tragwerk == 'Rahmentragwerk':
            str_frame.design(buildingProp, loads, materialProp,DataProp)

        if buildingProp.tragwerk == 'Braced Tube':
            str_bracedtube.design(buildingProp, loads, materialProp,DataProp)

        if buildingProp.tragwerk == 'Framed Tube':
            str_framedTube.design(buildingProp,loads,materialProp,DataProp)

        if buildingProp.tragwerk == 'Outrigger':
            str_outrigger.design(buildingProp,loads,materialProp,DataProp)
        
        building.buildingProp.calcFloorWeight(buildingProp,materialProp)


    def pushButton_einzelberechnung(self,buildingProp,materialProp,loads,DataProp):
        self.gui.progressBar.setValue(10)

        buildingProp.varPosOut_2D = False

        self.submit_buildingProp(buildingProp)
        self.submit_loads(loads, buildingProp)
        self.submit_materialProp(materialProp)
        self.gui.progressBar.setValue(20)

        self.mainCalculation(buildingProp,loads,materialProp,DataProp)
        self.gui.progressBar.setValue(30)

        self.showProfiles(buildingProp)
        self.gui.progressBar.setValue(40)

        self.plotMassDistribution(buildingProp)
        self.gui.progressBar.setValue(50)

        self.plotStiffnessDistribution(buildingProp,materialProp)
        self.gui.progressBar.setValue(60)
        
        self.plotInternalForcesCurve(buildingProp,loads)
        self.gui.progressBar.setValue(70)

        self.plotDeformationCurve(buildingProp)
        self.gui.progressBar.setValue(80)

        if buildingProp.tragwerk == 'Outrigger':
            self.plotalphaOutrigger(buildingProp)
            self.gui.progressBar.setValue(90)
            if buildingProp.Nachweis_GZT == True:
                self.plotAusnutzung(buildingProp, materialProp)
        self.gui.progressBar.setValue(100)


    def pushButton_multiberechnung(self,buildingProp,materialProp,loads,DataProp):
        self.gui.progressBar.setValue(10)

        self.gui.progressBar.setValue(20)

        n_min=self.gui.spinBox_n_min.value()
        n_max=self.gui.spinBox_n_max.value()
        varianten_n=self.gui.spinBox_varianten_n.value()
        delta_n=(n_max-n_min)/(varianten_n-1)       #Delta zum Erh??hen zwischen Iterationsschritten
        
        #n_0=buildingProp.n
        verh??ltnis_hw=self.gui.spinBox_w_max.value() #Maximale Verformung Verh??ltnis zu H??he

        buildingProp.multi_n=[]
        buildingProp.multi_h_total=[]
        buildingProp.multi_schlankheit=[]
        buildingProp.multi_A=[]

        buildingProp.multi_G_aussteifung=[]
        buildingProp.multi_G_au??enSt??tzen=[]
        buildingProp.multi_G_innenSt??tzen=[]
        buildingProp.multi_G_total=[]
        buildingProp.multi_G_decken=[]
        buildingProp.multi_G_totalOhnePFH=[]
        buildingProp.multi_eigenFrequenz=[]
        buildingProp.multi_moment = []

        buildingProp.multi_G_aussteifung_opt=[]
        buildingProp.multi_G_au??enSt??tzen_opt=[]
        buildingProp.multi_G_innenSt??tzen_opt=[]
        buildingProp.multi_G_total_opt=[]
        buildingProp.multi_G_decken_opt=[]
        buildingProp.multi_G_totalOhnePFH_opt=[]

        buildingProp.multi_p_opt=[]
        buildingProp.multi_p_opt_min=[]

        buildingProp.i_wechselNW='none'
        buildingProp.NW_ma??gebend_alt='Verformung'

        if self.gui.comboBox_parameter.currentText() == 'Outriggeranordnung ??ber die H??he':
            buildingProp.xi_True = False
        
        for i in range (0,varianten_n):
            self.submit_buildingProp(buildingProp)
            self.submit_loads(loads, buildingProp)
            self.submit_materialProp(materialProp)

            # Ver??ndern von n und den zugeh??rigen Variablen
            buildingProp.n=(int(n_max-i*delta_n))
            buildingProp.n_abschnitt = int(min(self.gui.spinBox_n_abschnitt.value(),buildingProp.n))
            buildingProp.multi_n.append(buildingProp.n)                                 # Liste der untersuchten Stockwersanzahlen
            buildingProp.h_total=buildingProp.h_geschoss*buildingProp.n
            buildingProp.multi_h_total.append(buildingProp.h_total)
            buildingProp.schlankheit=buildingProp.h_total/buildingProp.b_total
            buildingProp.multi_schlankheit.append(buildingProp.schlankheit)
            buildingProp.multi_A.append((buildingProp.b_raster*4)**2*buildingProp.n)    # Geschossfl??che*Geschossanzahl=Gesamtgeb??udefl??che

            buildingProp.x=int(buildingProp.n/buildingProp.n_abschnitt)                 # Anzahl der Abschnitte bei aktueller Geb??udeh??he
            loads.w_max=buildingProp.h_total/verh??ltnis_hw                              # max. Verformung

            print('-------',buildingProp.n,'-------')

            self.mainCalculation(buildingProp,loads,materialProp,DataProp)

            buildingProp.multi_G_aussteifung.append(buildingProp.G_aussteifung[-1])     # ??bernahme des Werts an Fu??punkt des Geb??udes (Summe des Gesamten), g=Ma??e
            buildingProp.multi_G_au??enSt??tzen.append(buildingProp.G_au??enSt??tzen[-1])
            buildingProp.multi_G_innenSt??tzen.append(buildingProp.G_innenSt??tzen[-1])
            buildingProp.multi_G_total.append(buildingProp.G_total[-1])
            buildingProp.multi_G_decken.append(buildingProp.G_decken)
            buildingProp.multi_G_totalOhnePFH.append(buildingProp.G_totalOhnePFH[-1])

            if buildingProp.NW_ma??gebend=='Tragf??higkeit' and buildingProp.NW_ma??gebend_alt=='Verformung':
                buildingProp.i_wechselNW=i

            buildingProp.NW_ma??gebend_alt=buildingProp.NW_ma??gebend

            buildingProp.multi_eigenFrequenz.append(buildingProp.eigenFrequenz)

            if buildingProp.tragwerk == 'Outrigger':
                momentmax=0
                for i in range(0,len(buildingProp.moment)):
                    if buildingProp.moment[i]>momentmax:
                        momentmax = buildingProp.moment[i]
                buildingProp.multi_moment.append(momentmax)
            else:
                buildingProp.multi_moment.append(loads.M[-1])

            if self.gui.comboBox_parameter.currentText()!='Keiner':     # Parametervariation
                self.multiberechnungParameter(buildingProp,materialProp,loads, DataProp)

                index=buildingProp.multiPar_G_total.index(min(buildingProp.multiPar_G_total))

                buildingProp.multi_G_aussteifung_opt.append(buildingProp.multiPar_G_aussteifung[index])
                buildingProp.multi_G_au??enSt??tzen_opt.append(buildingProp.multiPar_G_au??enSt??tzen[index])
                buildingProp.multi_G_innenSt??tzen_opt.append(buildingProp.multiPar_G_innenSt??tzen[index])
                buildingProp.multi_G_total_opt.append(buildingProp.multiPar_G_total[index])
                buildingProp.multi_G_decken_opt.append(buildingProp.multiPar_G_decken[index])
                buildingProp.multi_G_totalOhnePFH_opt.append(buildingProp.multiPar_G_totalOhnePFH[index])
                buildingProp.multi_eigenFrequenz[-1]=buildingProp.multiPar_eigenFrequenz[index]

                indexMax=int(self.gui.spinBox_varianten_p.value())-1

                index2=index

                while index2 < indexMax and buildingProp.multiPar_G_total[index2+1] == buildingProp.multiPar_G_total[index2]:
                    index2=index2+1
                    
                buildingProp.multi_p_opt.append(buildingProp.multi_p[index])
                buildingProp.multi_p_opt_min.append(buildingProp.multi_p[index2])
                

            self.gui.progressBar.setValue(20+30*int(i/varianten_n))

        self.plotResourceAnalysis(buildingProp,materialProp)
        self.gui.progressBar.setValue(55)
        self.plotResourceAnalysis_m2(buildingProp,materialProp)
        self.gui.progressBar.setValue(60)
        self.plotPFHAnalysis(buildingProp,materialProp)
        self.gui.progressBar.setValue(65)
        self.plotPFHAnalysis_m2(buildingProp,materialProp)
        self.gui.progressBar.setValue(70)
        self.plotEigenFrequency(buildingProp)
        self.gui.progressBar.setValue(75)
        self.plotMoment(buildingProp)
        self.gui.progressBar.setValue(80)
        
        if self.gui.comboBox_parameter.currentText()!='Keiner':
            self.plotOptimalParameter(buildingProp,materialProp)
            self.plotPFHAnalysisOptimalParameter_m2(buildingProp,materialProp)
            self.plotPFHAnalysisParameterComparisson_m2(buildingProp,materialProp)
            self.plotInfluenceOptimalParameter_m2(buildingProp,materialProp)

        self.gui.progressBar.setValue(100)


    def pushButton_multiberechnungParameter(self,buildingProp,materialProp,loads,DataProp):  #Parametervariation
        self.gui.progressBar.setValue(10)

        self.submit_buildingProp(buildingProp)
        self.submit_loads(loads, buildingProp)
        self.submit_materialProp(materialProp)
        self.gui.progressBar.setValue(20)

        self.multiberechnungParameter(buildingProp,materialProp,loads,DataProp)

        self.plotParameterAnalysis(buildingProp, materialProp)
        self.gui.progressBar.setValue(80)

        if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
            self.plotAlphaIst(buildingProp,materialProp)
            self.gui.progressBar.setValue(90)

        self.gui.progressBar.setValue(100)

    def pushButton_multiberechnungParameter_2D(self,buildingProp,materialProp,loads,DataProp):  # Parametervariation 2D
        self.gui.progressBar.setValue(10)

        self.checkIfPosOutVariing_2D(buildingProp)

        self.submit_buildingProp(buildingProp)
        self.submit_loads(loads, buildingProp)
        self.submit_materialProp(materialProp)
        self.gui.progressBar.setValue(20)

        self.multiberechnungParameter2D(buildingProp,materialProp,loads,DataProp)

        self.plot3DParameter(buildingProp,materialProp)
        self.gui.progressBar.setValue(100)


    def multiberechnungParameter(self,buildingProp,materialProp,loads,DataProp):
        p_min=self.gui.spinBox_p_min.value()
        p_max=self.gui.spinBox_p_max.value()
        varianten_p=int(self.gui.spinBox_varianten_p.value())
        delta_p=(p_max-p_min)/(varianten_p-1)

        buildingProp.parameter=self.gui.comboBox_parameter.currentText()

        buildingProp.multi_p = []

        buildingProp.multiPar_G_aussteifung=[]
        buildingProp.multiPar_G_au??enSt??tzen=[]
        buildingProp.multiPar_G_innenSt??tzen=[]
        buildingProp.multiPar_G_total=[]
        buildingProp.multiPar_G_decken=[]
        buildingProp.multiPar_G_totalOhnePFH=[]
        buildingProp.multiPar_eigenFrequenz=[]
        buildingProp.multiPar_t_rand_1 = []
        buildingProp.multiPar_t_eck_1 = []
        buildingProp.multiPar_t_rand_2 = []
        buildingProp.multiPar_t_eck_2 = []

        buildingProp.multiPar_w_EI = []
        buildingProp.multiPar_w_GA = []
        buildingProp.multiPar_w = []

        buildingProp.multiPar_alpha_ist_rand = []
        buildingProp.multiPar_alpha_ist_eck = []


        if buildingProp.parameter == 'Outriggeranordnung ??ber die H??he':
            if p_max > 1.0:
                print(str('Die gew??hlte Outriggeranordnung ist nicht verf??gbar. Bitte Zahlen zwischen 0 und 1 angeben!'))
                quit()

            if buildingProp.n_outrigger != 1:       # Esrtmal nur f??r einen Outrigger Optimierung
                print(str('Die gew??hlte Outriggeranzahl ist nicht m??glich. Bitte nur 1  Outriggerangeben!'))
                quit()
        
        for i in range (0,varianten_p):
            # Ver??ndern von n und den zugeh??rigen Variablen
            if buildingProp.parameter == 'Staffelung n_abschnitt':
                buildingProp.p=int(p_max-i*delta_p)
                buildingProp.n_abschnitt=min(buildingProp.p,buildingProp.n)
                buildingProp.x=int(buildingProp.n/buildingProp.p)

            elif buildingProp.parameter == 'Anzahl Stiele pro Raster (Framed Tube) n_stiele/b_raster':
                buildingProp.p=int(p_max-i*delta_p)
                buildingProp.n_stiele=buildingProp.p

            elif buildingProp.parameter == 'Anzahl Outrigger':
                buildingProp.p = int(p_max-i*delta_p)
                buildingProp.n_outrigger = buildingProp.p
                
            else:
                buildingProp.p=p_max-i*delta_p

            if buildingProp.parameter == 'Deckengewicht g_k1':
                g_k1_alt=self.gui.spinBox_gk1.value()
                loads.g_k1=buildingProp.p
                #self.gui.spinBox_gk1.setValue(buildingProp.p)       # Dadurch wird g_d neu berechnet
                
                materialProp.G_decke=100*buildingProp.p

                loads.gd = loads.gd+loads.gamma_g*(loads.g_k1-g_k1_alt)
                print('gd=',loads.gd)

            if buildingProp.parameter == 'Windlast q_b':  
                #loads.wk = self.gui.spinBox_wk.value()
                loads.wk=loads.wk/self.gui.spinBox_qb_k.value()*buildingProp.p

            if buildingProp.parameter == 'Geb??udeform psi_r':
                #self.gui.spinBox_psi_r.setValue(buildingProp.p)
                loads.wk=loads.wk/self.gui.spinBox_psi_r.value()*buildingProp.p

            if buildingProp.parameter == 'Verformungsverh??ltnis w_EI/w_GA':
                loads.w_verh??ltnis=buildingProp.p

            if buildingProp.parameter == 'Mindestwirksamkeit Framed Tube tube_wirksam':
                buildingProp.tube_wirksam_min=buildingProp.p

            if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
                buildingProp.alpha_outrigger = buildingProp.p

            if buildingProp.parameter == 'Steifigkeitsverh??ltnis Beta (Outrigger)':
                buildingProp.beta_outrigger = buildingProp.p

            if buildingProp.parameter == 'Outriggeranordnung ??ber die H??he':
                buildingProp.varPosOut = True
                buildingProp.xi_True = True
                buildingProp.xi = buildingProp.p

            
            buildingProp.multi_p.append(buildingProp.p)
            
            print('-------',buildingProp.p,'-------')

            self.mainCalculation(buildingProp,loads,materialProp,DataProp)

            if buildingProp.parameter == 'Steifigkeitsverh??ltnis Alpha (Outrigger)':
                buildingProp.multiPar_alpha_ist_rand.append(buildingProp.alpha_ist_rand)
                buildingProp.multiPar_alpha_ist_eck.append(buildingProp.alpha_ist_eck)

            buildingProp.multiPar_G_aussteifung.append(buildingProp.G_aussteifung[-1])
            buildingProp.multiPar_G_au??enSt??tzen.append(buildingProp.G_au??enSt??tzen[-1])
            buildingProp.multiPar_G_innenSt??tzen.append(buildingProp.G_innenSt??tzen[-1])
            buildingProp.multiPar_G_total.append(buildingProp.G_total[-1])
            buildingProp.multiPar_G_decken.append(buildingProp.G_decken)
            buildingProp.multiPar_G_totalOhnePFH.append(buildingProp.G_totalOhnePFH[-1])
            buildingProp.multiPar_eigenFrequenz.append(buildingProp.eigenFrequenz)
            buildingProp.multiPar_w_EI.append(buildingProp.w_EI[0])
            buildingProp.multiPar_w_GA.append(buildingProp.w_GA[0])
            buildingProp.multiPar_w.append(buildingProp.w)
            if buildingProp.tragwerk == 'Outrigger':
                for n, posOut in enumerate(buildingProp.posOut_abschnitt):
                    if n == 0:
                        buildingProp.multiPar_t_rand_1.append(buildingProp.t_randSt??tzen[posOut])
                        buildingProp.multiPar_t_eck_1.append(buildingProp.t_eckSt??tzen[posOut])
                    elif n == 1:
                        buildingProp.multiPar_t_rand_2.append(buildingProp.t_randSt??tzen[posOut])
                        buildingProp.multiPar_t_eck_2.append(buildingProp.t_eckSt??tzen[posOut])

    def multiberechnungParameter2D(self,buildingProp,materialProp,loads,DataProp):
        p_min_1=self.gui.spinBox_p_min_1.value()
        p_max_1=self.gui.spinBox_p_max_1.value()
        p_min_2=self.gui.spinBox_p_min_2.value()
        p_max_2=self.gui.spinBox_p_max_2.value()

        varianten_p_2D=int(self.gui.spinBox_varianten_p_2.value())
        delta_p_1=(p_max_1-p_min_1)/(varianten_p_2D-1)
        delta_p_2=(p_max_2-p_min_2)/(varianten_p_2D-1)

        buildingProp.parameter_2D=self.gui.comboBox_parameter_2.currentText()
        

        if buildingProp.tragwerk != 'Outrigger':
            print(str('Das gew??hlte Aussteifungssystem kann hier nicht analysiert werden. Bitte Outrigger ausw??hlen'))
            quit()

        if buildingProp.parameter_2D == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
            buildingProp.varPosOut_2D = True
            if p_max_1 > 1.0 or p_max_2 >1.0:
                print(str('Die gew??hlte Outriggeranordnung ist nicht verf??gbar. Bitte Zahlen zwischen 0 und 1 angeben!'))
                quit()

            if buildingProp.n_outrigger == 1:
                print(str('Die gew??hlte Outriggeranzahl ist nicht m??glich. Bitte nur 2 Outrigger angeben!'))
                quit()


        buildingProp.multi_p_1=[]
        buildingProp.multi_p_2=[]

        buildingProp.multiPar_2D_G_aussteifung = []
        buildingProp.multiPar_2D_G_au??enSt??tzen = []
        buildingProp.multiPar_2D_G_innenSt??tzen = []
        buildingProp.multiPar_2D_G_total = []
        buildingProp.multiPar_2D_G_decken = []
        buildingProp.multiPar_2D_G_totalOhnePFH = []

        buildingProp.multiPar_2D_w_EI = []
        buildingProp.multiPar_2D_w_GA = []
        buildingProp.multiPar_2D_w = []


        for i in range (0,varianten_p_2D):
            # Ver??ndern von n und den zugeh??rigen Variablen
            
            buildingProp.p1=p_max_1-i*delta_p_1

            for j in range(0,varianten_p_2D):

                buildingProp.p2=p_max_2-j*delta_p_2
                

                if buildingProp.parameter_2D == 'Steifigkeitsverh??ltnis Alpha und Beta (Outrigger)':
                    buildingProp.alpha_outrigger = buildingProp.p1
                    buildingProp.beta_outrigger = buildingProp.p2

                if buildingProp.parameter_2D == 'Outriggeranordnung ??ber die H??he f??r zwei Outrigger':
                    buildingProp.xi_1 = buildingProp.p1
                    buildingProp.xi_2 = buildingProp.p2

            
                print('-------','p1:',buildingProp.p1,'p2:',buildingProp.p2,'-------')

                self.mainCalculation(buildingProp,loads,materialProp,DataProp)

                buildingProp.multi_p_1.append(buildingProp.p1)
                buildingProp.multi_p_2.append(buildingProp.p2)

                buildingProp.multiPar_2D_G_aussteifung.append(buildingProp.G_aussteifung[-1])
                buildingProp.multiPar_2D_G_au??enSt??tzen.append(buildingProp.G_au??enSt??tzen[-1])
                buildingProp.multiPar_2D_G_innenSt??tzen.append(buildingProp.G_innenSt??tzen[-1])
                buildingProp.multiPar_2D_G_total.append(buildingProp.G_total[-1])
                buildingProp.multiPar_2D_G_decken.append(buildingProp.G_decken)
                buildingProp.multiPar_2D_G_totalOhnePFH.append(buildingProp.G_totalOhnePFH[-1])
            


    # ??bergabe der Werte in Objekte der Klassen buildingProp, loads und material:

    def submit_buildingProp(self,buildingProp):
        '''??bergibt Geb??udewerte aus GUI
        '''
        buildingProp.n = self.gui.spinBox_n.value()
        buildingProp.n_abschnitt = int(self.gui.spinBox_n_abschnitt.value())
        buildingProp.h_geschoss=self.gui.spinBox_h_geschoss.value()
        buildingProp.h_total=self.gui.spinBox_h_total.value()
        buildingProp.b_raster=self.gui.spinBox_b_raster.value()
        buildingProp.b_total=self.gui.spinBox_b_total.value()
        buildingProp.schlankheit=self.gui.spinBox_lambda.value()
        buildingProp.tragwerk=self.gui.comboBox_tragwerk.currentText()
        buildingProp.n_stiele=self.gui.spinBox_n_stiele.value()
        buildingProp.tube_wirksam_min=self.gui.spinBox_tube_wirksam.value()
        buildingProp.n_outrigger=self.gui.spinBox_n_outrigger.value()
        buildingProp.alpha_outrigger=self.gui.spinBox_alpha_outrigger.value()
        buildingProp.beta_outrigger=self.gui.spinBox_beta_outrigger.value()

        if buildingProp.n_abschnitt > buildingProp.n:     # n_abschnitt kann maximal so gro?? sein wie n
            buildingProp.n_abschnitt = buildingProp.n
        
        buildingProp.x=int(buildingProp.n/buildingProp.n_abschnitt)
        
        buildingProp.GA=[]
        buildingProp.I=[]


    def submit_loads(self,loads,buildingProp):
        '''??bergibt Lasten aus GUI'''

        loads.gd = self.gui.spinBox_gd.value()
        loads.qd = self.gui.spinBox_qd.value()
        loads.gk1 = self.gui.spinBox_gk1.value()
        loads.gk2 = self.gui.spinBox_gk2.value()
        loads.qk = self.gui.spinBox_qk.value()
        loads.gd_fassade = self.gui.spinBox_gd_fassade.value()
        loads.gk_fassade = self.gui.spinBox_gk_fassade.value()
        loads.gamma_g=self.gui.spinBox_gamma_g.value()
        loads.wk = self.gui.spinBox_wk.value()
        loads.gamma_w = self.gui.spinBox_gamma_w.value()
        loads.GK = self.gui.spinBox_GK.value()
        #if self.gui.checkBox_170_02.isChecked():
            #loads.GK = int(self.gui.comboBox_170_23.currentText())
            #loads.alpha_v = self.gui.doubleSpinBox_170_33.value()
            #loads.v_bk = self.gui.doubleSpinBox_170_43.value()
            #loads.D = self.gui.doubleSpinBox_170_53.value()
            #loads.dynamicAnalysis = True
        loads.Psi_q = self.gui.spinBox_Psi_q.value()
        loads.Psi_w = self.gui.spinBox_Psi_w.value()
        loads.gamma_gdyn = self.gui.spinBox_gamma_g_dyn.value()
        loads.gamma_qdyn = self.gui.spinBox_gamma_q_dyn.value()
        loads.alpha_a = self.gui.comboBox_alpha_a.currentText()
        loads.alpha_n = self.gui.comboBox_alpha_n.currentText()
            
        loads.w_max = buildingProp.h_total/self.gui.spinBox_w_max.value()
        loads.maxTeta_i = buildingProp.h_geschoss/self.gui.spinBox_maxTeta_i.value()
        loads.w_verh??ltnis = self.gui.spinBox_verhaeltnisVerformungen.value()


    def submit_materialProp(self,materialProp):
        '''??bergibt Materialkennwerte'''
        #Eingangsparameter Material und Profil

        materialProp.f = 1000*self.gui.spinBox_f.value() # N/mm?? -> kN/m??
        materialProp.E = 1000*self.gui.spinBox_E.value() # N/mm?? -> kN/m??
        materialProp.G = 1000*self.gui.spinBox_G.value() # N/mm?? -> kN/m??      
        materialProp.gamma = self.gui.spinBox_gamma.value()
        materialProp.t_min = self.gui.spinBox_tmin.value()
        materialProp.delta_t = self.gui.spinBox_delta_t.value() #Schrittweite Querschnittserh??hung
        materialProp.verh??ltnis_td = 1/self.gui.spinBox_verhaeltnis_td.value()
        materialProp.minderung_A = self.gui.spinBox_minderung_A.value()/100    # [%] -> [-]
        materialProp.minderung_I = self.gui.spinBox_minderung_I.value()/100    # [%] -> [-]
        materialProp.G_decke = self.gui.spinBox_G_decke.value()
        materialProp.GWP_decke = self.gui.spinBox_GWP_decke.value()
        materialProp.PET_decke = self.gui.spinBox_PET_decke.value()
        materialProp.GWP_tragwerk = self.gui.spinBox_GWP_tragwerk.value()
        materialProp.PET_tragwerk = self.gui.spinBox_PET_tragwerk.value()

        #print (materialProp.gamma)

#----------------------------------------------------------------------------------------------------------
def start():
    """Starting the application
    """
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()
