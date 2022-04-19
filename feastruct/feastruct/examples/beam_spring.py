from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.naturalfrequency import NaturalFrequency

from feastruct.solvers.feasolve import SolverSettings

import numpy as np

# ------------
# preprocessor
# ------------

# constants & lists
length = 5  # length of the beam
num_nodes = 6  # number of nodes to use
nodes = []  # list holding the node objects
elements = []  # list holding the element objects

# create 2d frame analysis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=1, ixx=1)

# create nodes
for i in range(num_nodes):
    nodes.append(analysis.create_node(coords=[length / (num_nodes - 1) * i]))

# create beam elements
for i in range(num_nodes - 1):
    elements.append(analysis.create_element(
        el_type='EB2-2D',
        nodes=[nodes[i], nodes[i+1]],
        material=steel,
        section=section
    ))

# add supports
freedom_case = cases.FreedomCase()
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
freedom_case.add_nodal_support(node=nodes[-1], val=0, dof=1)

# add springs
freedom_case.add_nodal_spring(node=nodes[0], val=50000000, dof=5) # rotational spring, not plotted


# add loads
load_case = cases.LoadCase()

load_case.add_nodal_load(node=nodes[0], val=100, dof=5)
load_case.add_nodal_load(node=nodes[2], val=-1000, dof=1)

# add analysis case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

# ------
# solver
# ------

settings = SolverSettings()
settings.linear_static.time_info = True

LinearStatic(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# ----
# post
# ----

# Support reactions, to check for validation
for support in analysis_case.freedom_case.items:
    reaction = support.get_reaction(analysis_case=analysis_case)
    print(reaction)

# Nodal displacement, to check for validation
for node in nodes:
    vec = node.get_displacements(analysis_case)
    print("u: " + str(vec[0]) + "; w: " + str(vec[1]) + "; phi: " + str(vec[5]))

analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=25)
# analysis.post.plot_frame_forces(analysis_case=analysis_case, axial=True)
# analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
# analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)

# --------------------------------
# natural frequency
# --------------------------------

# ----------------
# preprocessor
# ----------------
# an analysis case relates a support case to a load case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=cases.LoadCase())

# ----------------
# frequency solver
# ----------------
settings = SolverSettings()
settings.natural_frequency.time_info = True
settings.natural_frequency.num_modes = 1

solver = NaturalFrequency(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings)

# Manual solver, see feastruct/solvers/naturalfrequency.py, in order
# to extract mass/stiffness-matrix and eigenvectors      
# assign the global degree of freedom numbers
solver.assign_dofs()

# Get the global stiffness / mass matrix
(K, Kg) = solver.assemble_stiff_matrix()
M = solver.assemble_mass_matrix()

(f_ext, f_eq) = solver.assemble_fext(analysis_case=analysis_case)

# !apply spring condition
K_mod = solver.apply_spring(K=K, analysis_case=analysis_case)

# apply the boundary conditions
K_mod = solver.remove_constrained_dofs(K=K_mod, analysis_case=analysis_case)
M_mod = solver.remove_constrained_dofs(K=M, analysis_case=analysis_case)

# Solve for the eigenvalues
(fq_e, v) = solver.solve_eigenvalue(A=K_mod, M=M_mod, eigen_settings=settings.natural_frequency)

# compute natural frequencies in Hz
fq_e = np.sqrt(fq_e[0]) / 2 / np.pi

# In the case, the dominant eigenform is in longitudinal direction
# To see the influence of the spring, add a spring for DOF[0]
print("f: " + str(fq_e) + "Hz")