# ------------------------------------------------------------------------------
# Description:  Finite Elemente Berechnungen
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-27      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke

# Co-Author:    Carola Grupp
# Created:      2021-11-02  
# Projekt:      MAHS+ - MA Carola Grupp
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
import numpy as np
# For Static solver (feastruct)
from feastruct.pre.material import Material
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
import feastruct.fea.elements.frame2d as frame2d
import feastruct.fea.fea as fea
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings

from WindData import data

import pyLEK.plotters.plot2D as plt

# ------------------------------------------------------------------------------
# Abbreviations
# ------------------------------------------------------------------------------
# p...  load
# r...  response
# ms... model scale
# fs... full scale 
# L...  lift
# D...  drag
# M...  moment
# F...  force
# H...  (at) height of building
# sp..  sample
# fq...  frequency
# dn... direction
# ------------------------------------------------------------------------------


class feModel:
    """Class containing finite element analyis of the building
    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    """
    def __init__(self, buildingProp,loads,materialProp):
        """Inits the class, setting up calculation model
        :param buildingProp: full scale properties of the building
        :type buildingProp: :class:`~modelProp.buildProp`
        """

        # Setting up calculation model
        # ------------
        # preprocessor
        # ---------
        # constants & lists
        self.L = buildingProp.h_total                  # length of the beam [m]
        self.n=buildingProp.n                          # number of nodes between support and tip [-] gewählt: pro Etage ein Knoten
        if loads.windData == True:
            DataProp = data.DataProp()
            #self.z = DataProp.z_lev                     # coordinates for datapoints [m]

        # everything starts with the analysis object
        self.analysis = FrameAnalysis2D()


        ### NODES
        # ---------------
        # nodes are objects
        self.nodes = []
        for i in range(0,self.n+1):                     # 0 for tip
            node = self.analysis.create_node(coords=[0, self.L-i*buildingProp.h_geschoss])
            self.nodes.append(node)
            

          
        # materials and sections are objects
        # self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue, colour='w')

        # boundary conditions are objects
        self.freedom_case = cases.FreedomCase()
        
        # dof=0: Bewegung in x, dof=1: Bewegung in y, dof=2: Bewegung in z, dof=3: Rotation um x, dof=4 Rotation um y, dof=5 Rotation um z
        # Erlaubt nur Bewegung der Knoten in x und z Richtung und Rotation um y-Achse (-> zweidimensional)

        # Der support Knoten (letzter Knoten, daher -1) ist eingespannt und hat daher gar keine Freiheiten, keine Federsteifigkeit
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=0)      #add_nodal_support in fea>cases
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=1)
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=5)
        

        if buildingProp.tragwerk=='Outrigger':
            for j,i in enumerate(buildingProp.posOut):
                self.freedom_case.add_nodal_spring(node=self.nodes[i], val=buildingProp.K[j], dof=5)
                #self.freedom_case.add_nodal_support(node=self.nodes[i], val=0, dof=5, stiff=buildingProp.K[j])
                


        # Anzahl Abschnitte mit unterschiedlichen Steifigkeiten
        self.x=buildingProp.x

        # Adding loads
        self.load_case = cases.LoadCase()

        if loads.windData == False:
            for i in range(0,self.n+1):
                self.load_case.add_nodal_load(node=self.nodes[i], val=loads.F_p[i] , dof=0) # i+1 (support, tip)

        elif loads.windData == True:
            # Loop over all nodes
            for j in range(len(loads.F_p_j)):
                # Search for node by z coordinate
                F_p = np.mean(loads.F_p_j[j])        #[in KN]
                self.load_case.add_nodal_load(node=self.nodes[j], val=F_p , dof=0) 
                # Check applied loading
                # print('at z = '+ str(z_j) )
                # print('F_p  = '+ str(F_p) )

        # an analysis case relates a support case to a load case
        self.analysis_case = cases.AnalysisCase(freedom_case=self.freedom_case, load_case=self.load_case)

        ### BEAMS
        # ----------
        self.beams = []

        for i in range (0,self.x):
            ## SECTION
            self.section = Section(ixx=buildingProp.I[i])

            ## MATERIAL
            self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue[i], colour='w')

            # self.section = Section(area=1, ixx=1000)

            for j in range (i*buildingProp.n_abschnitt, (i+1)*buildingProp.n_abschnitt):
                # Add beams                               
                beam = self.analysis.create_element(el_type='EB2-2D', nodes=[self.nodes[j], self.nodes[j+1]], material=self.material, section=self.section)
                # print("Beam von j=" + str(j) + "bis j+1=" + str(j+1))
                self.beams.append(beam)

        if buildingProp.x*buildingProp.n_abschnitt<buildingProp.n:      # unterster Abschnitt, wenn Abschnitte nicht aufgehen
            self.section= Section(area=1, ixx=buildingProp.I[-1])
            z=buildingProp.x*buildingProp.n_abschnitt
            self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue[-1], colour='w')

            for j in range(0,buildingProp.n-buildingProp.x*buildingProp.n_abschnitt):   

                beam = self.analysis.create_element(el_type='EB2-2D', nodes=[self.nodes[z+j], self.nodes[z+j+1]], material=self.material, section=self.section)
                #print("Beam von j=" + str(z+j) + "bis j+1=" + str(z+j+1))
                self.beams.append(beam)

    def calcBendingMoment(self):
         # you can easily change the solver settings
        settings = SolverSettings()
        settings.linear_static.time_info = False

        # the linear static solver is an object and acts on the analysis object
        solver = LinearStatic(analysis=self.analysis, analysis_cases=[self.analysis_case], solver_settings=settings)
        solver.solve()

        #beamBernoulli = frame2d.EulerBernoulli2D_2N(nodes = self.nodes, material = self.material , section = self.section)
        #fea_int = fea.FiniteElement(nodes = self.nodes, material = self.material, efs = beamBernoulli.efs)
        
        bm = []
        #for i in range (0, len(self.nodes)):
        for el in solver.analysis.elements:
            bmd = el.get_bmd(2, self.analysis_case)
            bm_i1 = bmd[1][0]
            bm_i2 = bmd[1][1]
            #bm_i = max(abs(bmd[1][1]), abs(bmd[1][0]))
            #bm_i = bmd[1][1]
            bm.append(bm_i1)
            bm.append(bm_i2)
                 
        return bm



    def calcreactions(self):    #Fkt s.bcs
         # you can easily change the solver settings
        settings = SolverSettings()
        settings.linear_static.time_info = False

        # the linear static solver is an object and acts on the analysis object
        solver = LinearStatic(analysis=self.analysis, analysis_cases=[self.analysis_case], solver_settings=settings)
        solver.solve()

        reactions = []
        for support in self.analysis_case.freedom_case.items:
            reaction = support.get_reaction(analysis_case=self.analysis_case)
            reactions.append(reaction)
       
        return reactions 


    def calcStaticWindloadDeflection(self):        
        
        # you can easily change the solver settings
        settings = SolverSettings()
        settings.linear_static.time_info = False

        # the linear static solver is an object and acts on the analysis object
        solver = LinearStatic(analysis=self.analysis, analysis_cases=[self.analysis_case], solver_settings=settings)
        solver.solve()

        # ----
        # post
        # ----
        # there are plenty of post processing options!
        # self.analysis.post.plot_geom(analysis_case=analysis_case)
        # analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=1e2)
        # analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
        # self.analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
        # analysis.post.plot_reactions(analysis_case=analysis_case)
        
        # Support reactions, to check bending moment for validation
        for support in self.analysis_case.freedom_case.items:
            if support.dof in [5]:
                reaction = support.get_reaction(analysis_case=self.analysis_case)

        #print (reaction)
        # read out deformation at top 
        displacements = []
        for i in range (0, len(self.nodes)):    #len = Anzahl Stockwerke
            displacements.append(self.nodes[i].get_displacements(self.analysis_case)[0])
            #delta_tip = self.nodes[0].get_displacements(self.analysis_case)[0]
       
        return displacements #delta_tip   


    def calcInterstoryDrift(self,buildingProp):
        
        settings = SolverSettings()
        settings.linear_static.time_info = False

        # the linear static solver is an object and acts on the analysis object
        solver = LinearStatic(analysis=self.analysis, analysis_cases=[self.analysis_case], solver_settings=settings)
        solver.solve()
       
        # Support reactions, to check bending moment for validation
        
        # read out deformation at top 
        maxTeta_i=0
        Teta_i_EI=[]
        Teta_i_GA=[]

        for i in range (0,len(self.nodes)-1): # -1 damit für den untersten Knoten Schleife nicht ausgeführt wird (Kein Knoten unterhalb)
            
            #Verschiebung oben des Geschosses in x-Richtung (deswegen 0)
            w_EI_1 = self.nodes[i].get_displacements(self.analysis_case)[0]
            w_GA_1 = buildingProp.w_GA[i]
            
            #Verschiebung unten des Geschosses
            w_EI_2 = self.nodes[i+1].get_displacements(self.analysis_case)[0]
            w_GA_2 = buildingProp.w_GA[i+1]

            Teta_i_EI.append(w_EI_1-w_EI_2)
            Teta_i_GA.append(w_GA_1-w_GA_2)

            if maxTeta_i < Teta_i_EI[-1]+Teta_i_GA[-1]: #Verschiebung an aktuellem Geschoss (letzter Eintrag) 
                maxTeta_i=Teta_i_EI[-1]+Teta_i_GA[-1]
                maxTeta_i_EI=Teta_i_EI[-1]
                maxTeta_i_GA=Teta_i_GA[-1]

        buildingProp.Teta_i_EI=Teta_i_EI
        buildingProp.Teta_i_GA=Teta_i_GA

        buildingProp.Teta_i_EI.append(0)    # Hinzufügen einer Null ans Ende für das unterste Geschoss
        buildingProp.Teta_i_GA.append(0)

        w_EI_max=self.nodes[0].get_displacements(self.analysis_case)[0] #Verschiebung ganz oben

        return maxTeta_i, w_EI_max


    def getEigenfrequency(self):
        # ----------------
        # preprocessor
        # ----------------
        # an analysis case relates a support case to a load case
        analysis_case = cases.AnalysisCase(freedom_case=self.freedom_case, load_case=cases.LoadCase())

        # ----------------
        # frequency solver
        # ----------------
        settings = SolverSettings()
        settings.natural_frequency.time_info = True
        settings.natural_frequency.num_modes = 1

        solver = NaturalFrequency(
            analysis=self.analysis, analysis_cases=[analysis_case], solver_settings=settings)

        # Manual solver, see feastruct/solvers/naturalfrequency.py, in order
        # to extract mass/stiffness-matrix and eigenvectors      
        # assign the global degree of freedom numbers
        solver.assign_dofs()

        # Get the global stiffness / mass matrix
        (K, Kg) = solver.assemble_stiff_matrix()
        self.M = solver.assemble_mass_matrix()

        (self.f_ext, self.f_eq) = solver.assemble_fext(analysis_case=analysis_case)

        # apply spring condition
        self.K_mod = solver.apply_spring(K=K, analysis_case=analysis_case)

        # apply the boundary conditions
        self.K_mod = solver.remove_constrained_dofs(K=K, analysis_case=analysis_case)
        self.M_mod = solver.remove_constrained_dofs(K=self.M, analysis_case=analysis_case)

        # Solve for the eigenvalues
        (self.fq_e, self.v) = solver.solve_eigenvalue(A=self.K_mod, M=self.M_mod, eigen_settings=settings.natural_frequency)

        # compute natural frequencies in Hz
        self.fq_e = np.sqrt(self.fq_e[0]) / 2 / np.pi

        # Normalize Eigenvector acc. to Boggs 1991, 4.3.1, p.109
        self.v = self.v / self.v[0] * self.L

        # Store stiffness matrix as np array
        self.K_mod = self.K_mod.toarray()
        self.M_mod = self.M_mod.toarray()

        return self.fq_e

    def calcGeneralizedQuantitites(self):
        # Get generalized quantities
        self.K_gen = np.dot(np.dot(self.v.T, self.K_mod), self.v)[0][0]
        self.M_gen = np.dot(np.dot(self.v.T, self.M_mod), self.v)[0][0] 

        # To check, compute 
        # f_k = np.sqrt(K_gen/M_gen) / 2 / np.pi  d
        # print(f/f_k)
        # K_gen = 3 * E * I /(L^3) / 

   
def calcStress(element,buildingProp,loads,materialProp):   
    N_max= element.Ng_darüberliegend+(element.A_einzug*(loads.gd+loads.qd*element.alpha)+element.l_Fassade*loads.gd_Fassade)*buildingProp.n_abschnitt*element.i+element.A*materialProp.gamma*buildingProp.n_abschnitt*buildingProp.h_geschoss*1.35
    N_kombi= element.Ng_darüberliegend+(element.A_einzug*(loads.gd+loads.qd*element.alpha*loads.Psi_q)+element.l_Fassade*loads.gd_Fassade)*buildingProp.n_abschnitt*element.i+element.A*materialProp.gamma*buildingProp.n_abschnitt*buildingProp.h_geschoss*1.35
    
    sigma_Nmax=1.5*loads.Psi_w*loads.M[element.i-1]/element.W+N_max/element.A      #Abdeckung beider Kombinationen aus den veränderlichen Lasten durch Wind und Verkehr
    sigma_Mmax=1.5*loads.M[element.i-1]/element.W+N_kombi/element.A

    sigma=max(sigma_Nmax,sigma_Mmax)


    return sigma


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   