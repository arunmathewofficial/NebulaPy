import NebulaPy.src as nebula


pion = nebula.pion_silos(verbose=True)

silo_dir = '/home/tony/Desktop/NebulaPy/tests/wind-wind-jm'

# batch the silo file according to the time instant
silo_instant_set = pion.batch_silos(silo_dir, 'e7_WRwind_d1l5n256_v0750')

# get_chemistry will extract all chemistry information from
# the silo file into chemistry container
pion.get_chemistry(silo_instant_set[0])

# display chemistry container
print(pion.chemistry_container)

pion.get_cell_ne()





# get the number density of any species in the silo file









'''
for i in range(len(files[0])):
  datafile=[]
  for v in range(0,lev):
    datafile.append(files[v][i])
  print(i,datafile[0])
'''










'''
line_emission = nebula.emissionline('o_4', verbose=True)

line = 'OIV25' # Once can either put string or float
temperature = [1.e+3, 1.e+4]
ne = [0.1, 10]
line_emission.lineluminosity_1Dgrid(line=line, temperature=temperature, ne=ne, ns=0.0, dV=0.0)
print(line_emission.LineLuminosity)
'''


