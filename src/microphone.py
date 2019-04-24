import os

class Microphone(object):
    '''
        Deal with a single microphone
    '''
    def __init__(self, pos, id_val,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional', micro_config_path = './micro_config'):
        self.x_pos = pos[0]
        self.y_pos = pos[1]
        self.z_pos = pos[2]
        self.pos = pos
        self._id = str(id_val)
        self.orientation = orientation
        self.direction = os.path.join(micro_config_path, direction)
