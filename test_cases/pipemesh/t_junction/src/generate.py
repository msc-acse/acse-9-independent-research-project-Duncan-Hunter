from pipemesh import pipes
from pipemesh.icferst import auto_mpml

network = pipes.Network(0.1, 0.3, [0, 0, -1], 0.1)
network.add_t_junction([1, 0, 1], 0.1)
network.generate(filename="pipe", binary=True, write_info=True, write_xml=True, run_gui=False)

entry_phys_ids = network.get_inlet_outlet_phys_ids()
cyl_phys_ids = network.get_cyl_phys_ids()
inlets = entry_phys_ids[:2]
outlets = entry_phys_ids[2:]

vel = 0.02
inlet_velocities = network.get_velocities_vel_mag(inlets, vel)

options = auto_mpml.AutoMPML()
options.set_all(sim_name="junction_flow_test_case",
                msh_file="src/pipe",
                dump_ids=entry_phys_ids,
                density=1000,
                viscosity=1e-3,
                inlet_phys_ids=inlets, inlet_velocities=inlet_velocities,
                outlet_phys_ids=outlets,
                cyl_phys_ids=cyl_phys_ids,
                max_no_nodes=10000,
                min_mesh_size=0.08,
                finish_time=1.0,
                t_adapt_delay=0.5
                )
options.write_mpml("../3d_pipe_FEM")