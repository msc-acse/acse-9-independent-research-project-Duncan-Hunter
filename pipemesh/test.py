from pipemesh import pipes

network = pipes.Network(0.1, 0.4, [1, 0, 0], 0.2)

network.generate(run_gui=True)