



class xray:

    def __init__(self, min_xray_energy, max_xray_energy, elements, verbose=True,):

        self.min_energy = min_xray_energy
        self.max_energy = max_xray_energy
        self.elements = elements
        self.verbose = verbose
        self.setup()

    def setup(self):
        print(self.elements)


