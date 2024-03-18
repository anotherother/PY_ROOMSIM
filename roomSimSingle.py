import argparse
import numpy as np

from src.microphone import Microphone
from src.room import Room
from src.room_sim import RoomSim
from src.rt60_tools import rt60_to_absorption, get_rt60

def do_everything(room_dim, mic_positions, source_pos, rt60):
    absorption = rt60_to_absorption(room_dim, rt60)
    room = Room(room_dim, abs_coeff=absorption)
    mics = []
    for idx, mic in enumerate(mic_positions):
        temp_mic = Microphone(mic, idx,  \
            orientation=[0.0, 0.0, 0.0], direction='omnidirectional')
        mics.append(temp_mic)
    sim_rir = RoomSim(16000, room, mics, RT60=rt60)
    rir = sim_rir.create_rir(source_pos)
    return rir


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Config file')
    parser.add_argument('source_pos_x', help='Source x pos')
    parser.add_argument('source_pos_y', help='Source y pos')
    parser.add_argument('source_pos_z', help='Source z pos')
    parser.add_argument('out_file', help='File to write the RIR')
    args = parser.parse_args()
    source_pos = [float(args.source_pos_x), \
                    float(args.source_pos_y),\
                    float(args.source_pos_z)]
    sim_rir = RoomSim.init_from_config_file(args.config)
    rir = sim_rir.create_rir(source_pos)
    np.savetxt(args.out_file, rir)
