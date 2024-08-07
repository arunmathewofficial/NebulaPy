import glob
import sys
#sys.path.insert(0, "/Users/tony/.local/silo/lib/")
#import Silo
from pypion.ReadData import ReadData




class pion():
    '''
    This class is not an alternative to Pypion; rather, it is a
    bundle of methods useful for creating synthetic emission
    maps from the Silo file.
    '''

    '''
    def __init__(self, single, multiple, verbose):
        self.single = single # for single instant
        self.multiple = multiple # for multiple instants
    '''

    ######################################################################################
    # collate
    ######################################################################################
    def collate_files(self, dir, filebase):
        '''
        find files and put them in order of levels.
        :param path:
        :param filebase:
        :return: level_batch_files
        '''
        files = []
        level = 0
        for i in range(0, 20):
            seek = dir + "/" + filebase + "_level" + str(i).zfill(2) + "_0000.*.silo"
            level_files = sorted(glob.glob(seek))
            if len(level_files) < 1:
                break
            else:
                files.append(level_files)
                level = level + 1

        if len(files) < 1:
            print(f"empty directory error: no files in {dir}")
            quit()

        return files


    ######################################################################################
    # ?????
    ######################################################################################
    def get_all(self, file):
        dataio = ReadData(file)
        print(dataio)















    '''
    def get_basic_data(self):

        basic = self.get_1Darray("Density")
        return {'sim_time': basic['sim_time'],
                'min_extents': basic['min_extents'],
                'max_extents': basic['max_extents']}
    '''

