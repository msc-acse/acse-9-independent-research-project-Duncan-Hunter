from pipemesh import pipes

network = pipes.Network(0.5, 0.3, [0, 0, -1], 0.1)
network.add_curve([0, -1, 0], 0.5, 0.1)
network.add_cylinder(0.5, 0.1)

network.generate(filename="pipe", binary=False, write_info=True, write_xml=True)