import NebulaPy.src as nebula

EnergyBins = [[7.902470e+00,1.359840e+01],
              [1.359840e+01, 2.438310e+01],
              [2.438310e+01, 5.441780e+01],
              [5.441780e+01, 7.700000e+01]]

plot_dir = '/mnt/local/jm/code/arun/NebulaPy/problems/atlas_n04'
pion_format = '/mnt/local/jm/code/arun/NebulaPy/problems/atlas_n04'

atlas_sed = nebula.sed(
    database='/mnt/local/jm/code/arun/NebulaPy/NebulaPy-DB',
    energy_bins=EnergyBins,
    verbose=True,
    plot=None,
    pion=pion_format
)

#atlas_sed = nebula.sed(EnergyBins, verbose=True)
atlas_sed.CastelliKuruczAtlas(metallicity=0.0, gravity=4.5)
atlas_sed.CastelliKuruczAtlas(metallicity=0.0, gravity=4.0)
atlas_sed.CastelliKuruczAtlas(metallicity=0.0, gravity=3.5)

bb_sed = nebula.sed(
    database='/mnt/local/jm/code/arun/NebulaPy/NebulaPy-DB',
    energy_bins=EnergyBins,
    verbose=True,
    plot=None,
    pion=pion_format
)
bb_sed.Blackbody()



