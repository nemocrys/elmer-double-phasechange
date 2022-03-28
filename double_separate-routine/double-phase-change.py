from objectgmsh import Model, Shape
import gmsh
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver, run_elmerf90
from pyelmer.post import scan_logfile

def adjust_sif_for_phasechange2(boundary_name, sif="./case.sif"):
    """Function to adjust variable names etc to the new routine

    Args:
        boundary_name (str): name of the respective boundary
        sif (str, optional): Sif file path. Defaults to "./case.sif".
    """
    with open(sif) as f:
        data = f.readlines()

    in_bc_section = False
    for i in range(len(data)):
        if data[i] == f"! {boundary_name}\n":
            in_bc_section = True
        if in_bc_section and data[i] == "End\n":
            in_bc_section = False
        if in_bc_section:
            if "Equals PhaseSurface" in data[i]:
                data[i] = data[i].replace("Equals PhaseSurface", "Equals PhaseSurface2")
            if 'Real Procedure "SteadyPhaseChange" "MeltingHeat"' in data[i]:
                data[i] = data[i].replace('Real Procedure "SteadyPhaseChange" "MeltingHeat"', 'Real Procedure "SteadyPhaseChange2.so" "MeltingHeat"')
            if "Phase Change = Logical True" in data[i]:
                data[i] = data[i].replace("Phase Change = Logical True", "Phase Change 2 = Logical True")
    with open("./case.sif", "w") as f:
        f.writelines(data)

####################
# geometry modeling
####################
occ = gmsh.model.occ

model = Model()
crystal = Shape(model, 2, "crystal", [occ.add_rectangle(0, 0, 0, 1, 1)])
crystal.mesh_size = 0.1
melt = Shape(model, 2, "melt", [occ.add_rectangle(0, 1, 0, 1, 1)])
melt.mesh_size = 0.1
feed = Shape(model, 2, "feed", [occ.add_rectangle(0, 2, 0, 1, 1)])
feed.mesh_size = 0.1

crystal.set_interface(melt)
melt.set_interface(feed)

if_crystal_melt = Shape(model, 1, "if_crystal_melt", melt.get_interface(crystal))
if_feed_melt = Shape(model, 1, "if_feed_melt", melt.get_interface(feed))
bnd_melt = Shape(model, 1, "bnd_melt", [melt.right_boundary])
bnd_crys = Shape(model, 1, "bnd_crys", [crystal.bottom_boundary])
bnd_feed = Shape(model, 1, "bnd_feed", [feed.top_boundary])
bnd_others = Shape(model, 1, "bnd_others", model.symmetry_axis + [feed.left_boundary, feed.right_boundary, crystal.left_boundary, crystal.right_boundary])

model.synchronize()
model.make_physical()
model.set_const_mesh_sizes()
model.generate_mesh()
model.write_msh("case.msh")
# model.show()

####################
# simulation setup
####################

sim = elmer.load_simulation("axi-symmetric_steady", "elmer-config.yml")
solver_heat = elmer.load_solver("HeatSolver", sim, "elmer-config.yml")
solver_phase1 = elmer.load_solver("SteadyPhaseChange", sim, "elmer-config.yml")
solver_phase2 = elmer.load_solver("SteadyPhaseChange2", sim, "elmer-config.yml")
solver_mesh = elmer.load_solver("MeshUpdate", sim, "elmer-config.yml")
solver_out = elmer.load_solver("ResultOutputSolver", sim, "elmer-config.yml")
eqn_main = elmer.Equation(
    sim,
    "eqn_main",
    [
        solver_heat,
        solver_mesh
    ]
)
eqn_phase_1 = elmer.Equation(sim, "eqn_phase_1", [solver_phase1])
eqn_phase_2 = elmer.Equation(sim, "eqn_phase_2", [solver_phase2])

tin_l = elmer.load_material("tin_liquid", sim, "elmer-config.yml")
tin_s = elmer.load_material("tin_solid", sim, "elmer-config.yml")

crystal = elmer.Body(sim, "crystal", [crystal.ph_id])
crystal.material = tin_s
crystal.equation = eqn_main
melt = elmer.Body(sim, "melt", [melt.ph_id])
melt.material = tin_l
melt.equation = eqn_main
feed = elmer.Body(sim, "feed", [feed.ph_id])
feed.material = tin_s
feed.equation = eqn_main

bnd_melt = elmer.Boundary(sim, "bnd_melt", [bnd_melt.ph_id])
bnd_feed = elmer.Boundary(sim, "bnd_feed", [bnd_feed.ph_id])
bnd_crys = elmer.Boundary(sim, "bnd_crys", [bnd_crys.ph_id])

bnd_crys.fixed_temperature = 300
bnd_crys.mesh_update = [0, 0]
bnd_melt.fixed_heatflux = 1.5e4
bnd_melt.mesh_update = [0, None]
bnd_feed.fixed_temperature = 300
bnd_feed.mesh_update = [0, 0]

bnd_others = elmer.Boundary(sim, "bnd_others", [bnd_others.ph_id])
bnd_others.mesh_update = [0, None]

# phase change
phase_change_feed = elmer.Body(sim, "phase_change_feed", [if_feed_melt.ph_id])
phase_change_feed.material = tin_s
phase_change_feed.equation = eqn_phase_1
if_feed_melt = elmer.Boundary(sim, "if_feed_melt", [if_feed_melt.ph_id])
if_feed_melt.phase_change_steady = True
if_feed_melt.material = tin_s
if_feed_melt.normal_target_body = feed
if_feed_melt.phase_change_vel = 1e-6
if_feed_melt.phase_change_body = phase_change_feed

phase_change_crys = elmer.Body(sim, "phase_change_crys", [if_crystal_melt.ph_id])
phase_change_crys.material = tin_s
phase_change_crys.equation = eqn_phase_2
if_crystal_melt = elmer.Boundary(sim, "if_crystal_melt", [if_crystal_melt.ph_id])
if_crystal_melt.phase_change_steady = True
if_crystal_melt.material = tin_s
if_crystal_melt.normal_target_body = crystal
if_crystal_melt.phase_change_vel = 1e-6
if_crystal_melt.phase_change_body = phase_change_crys

sim.write_sif("./")
adjust_sif_for_phasechange2("if_crystal_melt")

####################
# simulation run
####################

run_elmerf90("SteadyPhaseChange2.F90", "SteadyPhaseChange2.so")
run_elmer_grid("./", "case.msh")
run_elmer_solver("./")

err, warn, stats = scan_logfile("./")
print("Errors:", err)
print("Warnings:", warn)
print("Statistic:", stats)
