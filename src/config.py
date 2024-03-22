from .microphone import Microphone
from .room import Room

class Config(object):
    '''
    Interface to read config files and put it to the right objects
    '''
    def __init__(self, config_file):
        self._file = config_file
        self.config = {}
        with open(config_file) as fid:
            for line in fid:
                line = line.strip()
                if line.startswith('%') or line == '':
                # This is a comment. Ignore
                    continue
                temp = line.split()
                try:
                    self.config[temp[0]] = [float(temp_) for temp_ in temp[1:]]
                except:
                    self.config[temp[0]] = [temp_ for temp_ in temp[1:]]
        self.config['Fs'] = int(self.config['Fs'][0])
        dict_keys = self.config.keys()
        self.sp_keys = [ke for  ke in dict_keys if ke.startswith('sp')]
        self.sd_keys = [ke for  ke in dict_keys if ke.startswith('sd')]
        self.so_keys = [ke for  ke in dict_keys if ke.startswith('so')]
        self.__verify_config()


    def __verify_config(self):
        assert 'room_size' in self.config, 'room_size not found in config'
        assert 'F_abs' in self.config, 'F_abs not found in config'
        assert 'Ax1' in self.config, 'Ax1 not found in config'
        assert 'Ax2' in self.config, 'Ax2 not found in config'
        assert 'Ay1' in self.config, 'Ay1 not found in config'
        assert 'Ay2' in self.config, 'Ay2 not found in config'
        assert 'Az1' in self.config, 'Az1 not found in config'
        assert 'Az2' in self.config, 'Az2 not found in config'
        assert 'sp1' in self.config, 'sp1 not found in config'
        assert 'sd1' in self.config, 'sd1 not found in config'
        assert 'so1' in self.config, 'so1 not found in config'
        assert len(self.sp_keys) == len(self.sd_keys) == len(self.so_keys), \
            'sp, sd and so are not of same length'

    def create_room_et_mic_objects(self):
        room_size = [float(_) for _ in self.config['room_size']]
        F_abs = [float(_) for _ in self.config['F_abs']]
        Ax1 = [float(_) for _ in self.config['Ax1']]
        Ax2 = [float(_) for _ in self.config['Ax2']]
        Ay1 = [float(_) for _ in self.config['Ay1']]
        Ay2 = [float(_) for _ in self.config['Ay2']]
        Az1 = [float(_) for _ in self.config['Az1']]
        Az2 = [float(_) for _ in self.config['Az2']]
        room = Room(room_size, F_abs, [Ax1, Ax2, Ay1, Ay2, Az1, Az2])
        mics = []
        for mic_idx in range(len(self.sp_keys)):
            mic_idx += 1
            _xp, _yp, _zp = self.config['sp'+str(mic_idx)]
            orientation = self.config['so'+str(mic_idx)]
            direction = self.config['sd'+str(mic_idx)][0].replace("'",'')
            mics.append(Microphone([_xp, _yp, _zp], mic_idx,\
                                  orientation = orientation, direction = direction))
        return[self.config['Fs'], room, mics]
