from pipemesh import pipes

network = pipes.Network(1.2, 0.4, [0, 0, -1], 0.1)

network.generate("pipe", binary=True, write_info=True, write_xml=True)