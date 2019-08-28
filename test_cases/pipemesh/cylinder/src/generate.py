from pipemesh import pipes
from pipemesh.icferst import auto_mpml


network = pipes.Network(1.2, 0.4, [0, 0, -1], 0.1)

network.generate("pipe", binary=True, write_info=True, write_xml=False)

entry_phys_ids = network.get_inlet_outlet_phys_ids()
cyl_phys_ids = network.get_cyl_phys_ids()
inlets = entry_phys_ids[:1]
outlets = entry_phys_ids[1:]

vel = 0.02
inlet_velocities = network.get_velocities_vel_mag(inlets, vel)

options = auto_mpml.AutoMPML()
options.set_all(sim_name="3d_pipe_cylinder_test_case",
                msh_file="src/pipe",
                dump_ids=entry_phys_ids,
                inlet_phys_ids=inlets, inlet_velocities=inlet_velocities,
                density=1000,
                viscosity=1e-3,
                outlet_phys_ids=outlets,
                cyl_phys_ids=cyl_phys_ids,
                max_no_nodes=10000,
                min_mesh_size=0.08,
                finish_time=1.0,
                t_adapt_delay=0.5
                )
options.write_mpml("../3d_pipe_FEM")