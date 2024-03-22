import numpy as np

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