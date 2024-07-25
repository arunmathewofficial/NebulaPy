import NebulaPy.src as nebula


line_emission = nebula.emissionline('o_4')
temperature = 1.e+9
ne = 1.e+4
line_emission.Luminosity_1DGrid(0.0, temperature, ne, 0.0, 0.0)


