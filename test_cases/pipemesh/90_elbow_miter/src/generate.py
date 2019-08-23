from pipemesh import pipes

network = pipes.Network(0.5, 0.3, [0, 0, -1], 0.1)
network.add_mitered([0, 1, 0], 0.06)
network.add_cylinder(0.5, 0.1)

network.generate(filename="pipe", binary=True, write_info=True, write_xml=True, run_gui=False)