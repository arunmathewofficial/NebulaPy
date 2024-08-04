import NebulaPy.src as nebula

EnergyBins = [[7.902470e+00, 1.126030e+01], [1.126030e+01, 1.359840e+01],
                 [1.359840e+01, 1.453410e+01], [1.453410e+01, 1.619920e+01],
                 [1.619920e+01, 2.156450e+01], [2.156450e+01, 2.438310e+01],
                 [2.438310e+01, 2.960130e+01], [2.960130e+01, 3.065100e+01],
                 [3.065100e+01, 3.512110e+01], [3.512110e+01, 4.096300e+01],
                 [4.096300e+01, 4.514180e+01], [4.514180e+01, 4.788780e+01],
                 [4.788780e+01, 5.441780e+01], [5.441780e+01, 6.342330e+01],
                 [6.342330e+01, 7.700000e+01]]

plot_dir = '/Users/tony/Desktop/NebulaPy/tests'
pion_formar_path = '/Users/tony/Desktop/NebulaPy/tests'

atlas_sed = nebula.sed(EnergyBins, verbose=True, plot=plot_dir, pion=pion_formar_path)
#atlas_sed = nebula.sed(EnergyBins, verbose=True)

atlas_sed.CastelliKuruczAtlas(0.0, 4.5)
#print(atlas_sed.container)