import numpy as np

class Room(object):
    
    def __init__(self, dim, F_abs=None, abs_coeff=None):
        self.x_val = dim[0]
        self.y_val = dim[1]
        self.z_val = dim[2]
        self.room_size = np.array(dim)
        self.freq_dep_absorption = {}
        if F_abs is None:
            self.freq_dep_absorption['F_abs'] = np.array([125, 250, 500, 1000, 2000, 4000, 8000])
        else:
            self.freq_dep_absorption['F_abs'] = np.array(F_abs)
        if abs_coeff is None:
            self.__set_absorption()
        else:
            if isinstance(abs_coeff, float) or isinstance(abs_coeff, int):
                self.__set_absorption(abs_val=abs_coeff)
            else:
                self.freq_dep_absorption['Ax1'] = np.array(abs_coeff[0])
                self.freq_dep_absorption['Ax2'] = np.array(abs_coeff[1])
                self.freq_dep_absorption['Ay1'] = np.array(abs_coeff[2])
                self.freq_dep_absorption['Ay2'] = np.array(abs_coeff[3])
                self.freq_dep_absorption['Az1'] = np.array(abs_coeff[4])
                self.freq_dep_absorption['Az2'] = np.array(abs_coeff[5])

    def __set_absorption(self, abs_val=0.671):
        self.freq_dep_absorption['Ax1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ax2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ay1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Ay2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Az1'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
        self.freq_dep_absorption['Az2'] = np.array([abs_val] * len(self.freq_dep_absorption['F_abs']))
