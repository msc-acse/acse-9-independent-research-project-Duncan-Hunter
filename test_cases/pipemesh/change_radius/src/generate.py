from pipemesh import pipes

network = pipes.Network(0.1, 0.2, [0, 0, -1], 0.1)
network.add_change_radius(1, 0.4, 0.9, 0.1)
network.add_cylinder(0.1, 0.1)

network.generate("pipe", binary=True, write_info=True, write_xml=True)