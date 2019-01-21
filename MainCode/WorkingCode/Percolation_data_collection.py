import Atrium as AC

def CMP2D_timestep_perc(Atrium):
    """A single timestep"""
    Atrium.SinusRhythm()
    Atrium.Relaxing()
    Atrium.Conduct()
    while len(Atrium.states[0]) != 0:
        Atrium.Relaxing()
        Atrium.Conduct()
        Atrium.t += 1
        for i in Atrium.states[0]:
            if i % Atrium.size == Atrium.size-1:
                return 1

    return 0


Atrium = AC.Atrium()