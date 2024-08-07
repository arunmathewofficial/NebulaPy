import NebulaPy.src as nebula
import sys
#sys.path.insert(0, "/Users/tony/.local/silo/lib")
#import Silo

# bundle up the silo files
silo = nebula.pion()
silo_dir = '/Users/tony/Desktop/NebulaPy/tests/wind-wind-jm'
files = silo.collate_files(dir=silo_dir, filebase='e7_WRwind_d1l5n256_v0750')
silo.get_all(files[0][0])






'''
line_emission = nebula.emissionline('o_4', verbose=True)

line = 'OIV25' # Once can either put string or float
temperature = [1.e+3, 1.e+4]
ne = [0.1, 10]
line_emission.lineluminosity_1Dgrid(line=line, temperature=temperature, ne=ne, ns=0.0, dV=0.0)
print(line_emission.LineLuminosity)
'''


