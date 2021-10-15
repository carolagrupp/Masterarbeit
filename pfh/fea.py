# ------------------------------------------------------------------------------
# Description:  Finite Elemente Berechnungen
#
# ------------------------------------------------------------------------------
# Author:       st169687@stud.uni-suttgart.de
# Created:      2021-04-27      (YYYY-MM-DD)
# Projekt:      Premium for Height - MA Christian Engelke
# ------------------------------------------------------------------------------
# Sources:
# ------------------------------------------------------------------------------
# Imports:      
import numpy as np
# For Static solver (feastruct)
from feastruct.pre.material import Material
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings

import plotters.plot2D as plt

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
        self.n=buildingProp.n     # number of nodes between support and tip [-] gewählt: pro Etage ein Knoten
        
        # everything starts with the analysis object
        self.analysis = FrameAnalysis2D()

        # nodes are objects
        self.nodes = []
        for i in range(0,self.n+1): # 0 for tip
            node = self.analysis.create_node(coords=[0, self.L-i*buildingProp.h_geschoss])
            self.nodes.append(node)
            
            # print(node.y) # Koordinaten - Kontrolle

        # bs/ce: Test: Anzahl Beams = Anzahl nodes - 1, dass passt
        #print("Anzahl nodes = " + str(len(self.nodes)))
          
        # materials and sections are objects
        # self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue, colour='w')

        # boundary conditions are objects
        self.freedom_case = cases.FreedomCase()
        
        # dof=0: Bewegung in x, dof=1: Bewegung in y, dof=2: Bewegung in z, dof=3: Rotation um x, dof=4 Rotation um y, dof=5 Rotation um z
        # Erlaubt nur Bewegung der Knoten in x und z Richtung und Rotation um y-Achse (-> zweidimensional)

        # Der support Knoten (letzter Knoten, daher -1) ist eingespannt und hat daher gar keine Freiheiten
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=0)
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=1)
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=5)

        # Anzahl Abschnitte mit unterschiedlichen Steifigkeiten
        self.x=buildingProp.x

        # Adding loads
        self.load_case = cases.LoadCase()

        for i in range(0,self.n+1):
            self.load_case.add_nodal_load(node=self.nodes[i], val=loads.F_p[i] , dof=0) # i+1 (support, tip)

        # an analysis case relates a support case to a load case
        self.analysis_case = cases.AnalysisCase(freedom_case=self.freedom_case, load_case=self.load_case)

        self.beams = []

        for i in range (0,self.x):
  
            self.section = Section(ixx=buildingProp.I[i])
            self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue[i], colour='w')

            # self.section = Section(area=1, ixx=1000)

            for j in range (i*buildingProp.n_abschnitt, (i+1)*buildingProp.n_abschnitt):
                # bs/ce: 
                # Auch hier passt was noch nicht so ganz
                # Beams werden von j bis j+1 erstellt. Das passt, aber sowie ich das sehe 
                # gibt es nodes zwischen denen keine beams erstellt werden
                               
                beam = self.analysis.create_element(el_type='EB2-2D', nodes=[self.nodes[j], self.nodes[j+1]], material=self.material, section=self.section)
                # print("Beam von j=" + str(j) + "bis j+1=" + str(j+1))
                self.beams.append(beam)

        if buildingProp.x*buildingProp.n_abschnitt<buildingProp.n:      # unterster Abschnitt
            self.section= Section(area=1, ixx=buildingProp.I[-1])
            z=buildingProp.x*buildingProp.n_abschnitt
            self.material = Material("Dummy", materialProp.E, 0.2, buildingProp.mue[-1], colour='w')

            for j in range(0,buildingProp.n-buildingProp.x*buildingProp.n_abschnitt):   

                beam = self.analysis.create_element(el_type='EB2-2D', nodes=[self.nodes[z+j], self.nodes[z+j+1]], material=self.material, section=self.section)
                #print("Beam von j=" + str(z+j) + "bis j+1=" + str(z+j+1))
                self.beams.append(beam)


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
        for i in range (0, len(self.nodes)):
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

            w_EI_1 = self.nodes[i].get_displacements(self.analysis_case)[0]
            w_GA_1 = buildingProp.w_GA[i]
            
            w_EI_2 = self.nodes[i+1].get_displacements(self.analysis_case)[0]
            w_GA_2 = buildingProp.w_GA[i+1]

            Teta_i_EI.append(w_EI_1-w_EI_2)
            Teta_i_GA.append(w_GA_1-w_GA_2)

            if maxTeta_i < Teta_i_EI[-1]+Teta_i_GA[-1]:
                maxTeta_i=Teta_i_EI[-1]+Teta_i_GA[-1]
                maxTeta_i_EI=Teta_i_EI[-1]
                maxTeta_i_GA=Teta_i_GA[-1]

        buildingProp.Teta_i_EI=Teta_i_EI
        buildingProp.Teta_i_GA=Teta_i_GA

        buildingProp.Teta_i_EI.append(0)    # Hinzufügen einer Null ans Ende für das unterste Geschoss
        buildingProp.Teta_i_GA.append(0)

        w_EI_max=self.nodes[0].get_displacements(self.analysis_case)[0]

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