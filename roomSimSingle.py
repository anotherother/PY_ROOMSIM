import argparse
import numpy as np

from src.microphone import Microphone
from src.room import Room
from src.room_sim import RoomSim

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

def get_rt60(F_abs, room_size, A):
    '''
    Get RT 60 given the room characteristics
    '''
    m_air = 6.875e-4*(F_abs.T/1000)**(1.7)
    # attenuation factors for one metre travelled in air
    room_size = np.array(room_size)
    atten_air = np.exp(-0.5*m_air).T
    Lx = room_size[0]
    Ly = room_size[1]
    Lz = room_size[2]
    #Volume of room m^3
    V_room=Lx*Ly*Lz
    area_xz=Lx*Lz
    area_yz=Ly*Lz
    area_xy=Lx*Ly
    total_area = 2*(area_xz+area_yz+area_xy)# Total area of shoebox room surfaces
    # Effective absorbing area of room surfaces at each frequency
    Se=area_yz*(A[0]+A[1])+area_xz*(A[2]+A[3])+area_xy*(A[5]+A[4])
    a_bar=Se/total_area # Mean absorption of each room surface
    # Norris-Eyring estimate adjusted for air absorption
    RT60=0.1611*V_room/(4*m_air.T*V_room-total_area*np.log(1-a_bar))
    return RT60

def rt60_to_absorption(room_obj_dim, rt60):
    '''
    Norris-Eyring formula %%
     Converts a given reverberation time into a single absorption coefficient for all surfaces 
    '''
    room_vol = np.prod(room_obj_dim)
    area_xz=room_obj_dim[0] * room_obj_dim[2]
    area_yz=room_obj_dim[1] * room_obj_dim[2]
    area_xy=room_obj_dim[0] * room_obj_dim[1]
    total_area =2*(area_xz+area_yz+area_xy) # Total area of shoebox room surfaces
    absorption = 1-np.exp(-0.1611*room_vol/(total_area*rt60))
    return absorption


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
