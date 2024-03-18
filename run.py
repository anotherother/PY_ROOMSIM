import soundfile as sf
from src.overlap_add_filter import OverlapAddFilter
import roomSimSingle
from utils.audio_loader import AudioLoader

def main():
    rt60 = 0.0 # in seconds
    room_dim = [3.5, 4, 2.5] # in meters
    mic_pos1 = [2, 2, 1] # in  meters
    mic_pos2 = [2, 2, 1.5] # in  meters
    source_pos = [2, 2, 1.2] # in  meters
    sampling_rate = 16000

    absorption = roomSimSingle.rt60_to_absorption(room_dim, rt60)

    room = roomSimSingle.Room(room_dim, abs_coeff=absorption)
    mic1 = roomSimSingle.Microphone(mic_pos1, 1, orientation=[0.0, 0.0, 0.0],
                                    direction='omnidirectional', micro_config_path='./micro_config')

    mic2 = roomSimSingle.Microphone(mic_pos2, 2, orientation=[0.0, 0.0, 0.0],
                                    direction='cardioid', micro_config_path='./micro_config')
    mics = [mic1, mic2]
    sim_rir = roomSimSingle.RoomSim(sampling_rate, room, mics, RT60=rt60)
    rir = sim_rir.create_rir(source_pos)

    al = AudioLoader(sampling_rate)
    fs, data = al.load_audio('./wav/4.wav')
    print(data.shape)
    data_rev = OverlapAddFilter(rir[:,0], data)
    sf.write('./wav/data_rev4.wav', data_rev, fs)

if __name__ == '__main__':
    main()