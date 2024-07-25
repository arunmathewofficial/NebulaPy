import NebulaPy.src as nebula


line_emission = nebula.emissionline('o_4', verbose=True)

line = 'OIV25' # Once can either put string or float
temperature = [1.e+3, 1.e+4]
ne = [0.1, 10]
line_emission.lineluminosity_1Dgrid(line=line, temperature=temperature, ne=ne, ns=0.0, dV=0.0)
print(line_emission.LineLuminosity)



