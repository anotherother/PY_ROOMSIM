import scipy.signal as scipy_sig
import numpy as np

from src.acoustics import get_rt60
from src.config import Config

_TEXT_CACHE = {}


def _load_text(path):
    cached = _TEXT_CACHE.get(path)
    if cached is None:
        cached = np.loadtxt(path)
        _TEXT_CACHE[path] = cached
    return cached


class RoomSim(object):
    '''
    Class to handle RIR creation:
        Input
        -----
        room_config : Roomconfig object
    '''

    def __init__(self, fs, room, mics, RT60=None):
        self._do_init(fs, room, mics, RT60)
        self.verify_positions()

    def verify_positions(self):
        '''
        Method to verify if all the microphones are inside the room
        '''
        for mic in self.mics:
            assert mic.x_pos < self.room.x_val,\
                    mic._id+' x position is outside the room'
            assert mic.y_pos < self.room.y_val,\
                    mic._id+' y position is outside the room'
            assert mic.z_pos < self.room.z_val,\
                    mic._id+' z position is outside the room'


    @classmethod
    def init_from_config_file(cls, room_config_file, RT60=None):
        '''
        constructor to read config file and initialize an instance
        '''
        config = Config(room_config_file)
        sample_rate, room, mics = config.create_room_et_mic_objects()
        obj = cls(sample_rate, room, mics, RT60)
        return obj

    def _do_init(self, fs, room, mics, RT60):
        self.sampling_rate = fs
        self.room = room
        self.mics = mics
        mic_count = 0
        for mic in self.mics:
            mic_count += 1
            mic._id = str(mic_count)
        self.channels = len(mics)
        self.room_size = room.room_size
        self.F_abs = room.freq_dep_absorption['F_abs']
        Ax1 = room.freq_dep_absorption['Ax1']
        Ax2 = room.freq_dep_absorption['Ax2']
        Ay1 = room.freq_dep_absorption['Ay1']
        Ay2 = room.freq_dep_absorption['Ay2']
        Az1 = room.freq_dep_absorption['Az1']
        Az2 = room.freq_dep_absorption['Az2']
        self.A = np.array([Ax1, Ax2, Ay1, Ay2, Az1, Az2])
        self.A = self.A[:, self.F_abs<=self.sampling_rate/2.0]
        self.F_abs = self.F_abs[self.F_abs<=self.sampling_rate/2.0]
        if self.F_abs[0] != 0:
            self.A = np.vstack((self.A.T[0], self.A.T)).T
            self.F_abs = np.hstack((0, self.F_abs))
        if self.F_abs[-1] != self.sampling_rate/2.0:
            self.A = np.vstack((self.A.T, self.A.T[-1]))
            self.F_abs = np.hstack((self.F_abs, self.sampling_rate/2.0))

        self.tm_sensor = np.zeros((self.channels, 3, 3))
        self.sensor_xyz = np.zeros((self.channels, 3))
        self.sensor_off = np.zeros((self.channels, 3))
        for idx, mic in enumerate(self.mics):
            self.sensor_xyz[idx, :] = mic.pos
            self.sensor_off[idx, :] = mic.orientation
            self.tm_sensor[idx, :, :] = self.__create_tm(\
                                    self.__create_psi_theta_phi(mic.orientation))
        if RT60 is None:
            self.RT60 = get_rt60(self.F_abs, self.room_size, self.A)
        else:
            self.RT60 = np.array([RT60] * len(self.F_abs))



    def create_rir(self, source_xyz, source_off=None, source_dir=None):
        '''
        Create the RIR
        source_xyz : list containing xyz position of the source
        source_off: 3 x 1 list representing the source orientation (azimuth,
        elevation, roll)
        source_dir: source directivity np txt file of dimension 181 x 361
        '''
        source_xyz = np.array(source_xyz)
        if source_dir is None:
            source_dir = np.ones((181, 361))
        else:
            source_dir = _load_text(source_dir)
        if source_off is None:
            source_off = np.zeros(source_xyz.shape)

        tm_source = self.__create_tm(self.__create_psi_theta_phi(source_off))
        sampling_period = 1.0 / self.sampling_rate
        nyquist = self.sampling_rate / 2.0
        Fs_c = self.sampling_rate / 343.0
        H_length = np.floor(np.max(self.RT60) * self.sampling_rate)
        range_ = H_length / Fs_c
        Lx = self.room_size[0]
        Ly = self.room_size[1]
        Lz = self.room_size[2]
        order_x = np.ceil(range_ / (2 * Lx))
        order_y = np.ceil(range_ / (2 * Ly))
        order_z = np.ceil(range_ / (2 * Lz))
        delay_s = Fs_c * np.sqrt(np.sum((source_xyz.T - self.sensor_xyz) ** 2, axis=1))
        H_length = int(np.max((H_length, np.ceil(np.max(np.max(delay_s))) + 200)))
        N_frac = 32
        Tw = N_frac * sampling_period
        Two_pi_Tw = 2 * np.pi / Tw
        t = np.arange(-Tw / 2, Tw / 2 + sampling_period, sampling_period)
        w = 2 * np.pi * 20
        r1 = np.exp(-w * sampling_period)
        r2 = np.exp(-w * sampling_period)
        b1 = -(1 + r2)
        b2 = np.copy(r2)
        a1 = 2 * r1 * np.cos(w * sampling_period)
        a2 = -r1 * r1
        HP_gain = (1 - b1 + b2) / (1 + a1 - a2)
        b_HP = [1, b1, b2] / HP_gain
        a_HP = [1, -a1, -a2]
        Two_Lx = 2 * self.room_size[0]
        Two_Ly = 2 * self.room_size[1]
        Two_Lz = 2 * self.room_size[2]
        isource_ident = np.array([
            [-1, -1, -1],
            [-1, -1, +1],
            [-1, +1, -1],
            [-1, +1, +1],
            [+1, -1, -1],
            [+1, -1, +1],
            [+1, +1, -1],
            [+1, +1, +1],
        ])
        surface_coeff = np.array([
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1],
        ])
        qq = surface_coeff[:, 0]
        jj = surface_coeff[:, 1]
        kk = surface_coeff[:, 2]
        F_abs_N = self.F_abs / nyquist
        N_refl = int(2 * np.round(nyquist / self.F_abs[1]))
        Half_I = int(N_refl / 2)
        Half_I_plusone = Half_I + 1
        window = 0.5 * (1 - np.cos(2 * np.pi * np.arange(0, N_refl + 1).T / N_refl))
        isource_xyz, refl = self._image_sources(
            source_xyz, order_x, order_y, order_z, H_length, Fs_c,
            Two_Lx, Two_Ly, Two_Lz, isource_ident, qq, jj, kk,
        )
        n_images = isource_xyz.shape[1]
        H = np.zeros((H_length, self.channels))
        m_air = 6.875e-4 * (self.F_abs / 1000) ** (1.7)
        atten_air = np.exp(-0.5 * m_air).T
        freq_grid = (1.0 / Half_I) * np.arange(Half_I + 1)
        for mic in self.mics:
            sensor_dir = _load_text(mic.direction + '.txt')
            sensor_no = int(mic._id) - 1
            if n_images:
                self._render_sensor(
                    H, sensor_no, isource_xyz, refl, atten_air, sensor_dir,
                    source_dir, tm_source, Fs_c, t, Two_pi_Tw, sampling_period,
                    N_frac, N_refl, Half_I_plusone, window, F_abs_N,
                    freq_grid, n_images,
                )
            H[:, sensor_no] = scipy_sig.lfilter(b_HP, a_HP, H[:, sensor_no])
        return H

    def _image_sources(self, source_xyz, order_x, order_y, order_z, H_length, Fs_c,
                       Two_Lx, Two_Ly, Two_Lz, isource_ident, qq, jj, kk):
        xx_yy_zz = np.array([
            isource_ident[:, 0] * source_xyz[0],
            isource_ident[:, 1] * source_xyz[1],
            isource_ident[:, 2] * source_xyz[2],
        ])
        B = np.sqrt(1 - self.A)
        bx1, bx2, by1, by2, bz1, bz2 = B
        coordinates = []
        reflections = []
        for n in np.arange(-order_x, order_x + 1, 1):
            bx2_abs_n = bx2 ** np.abs(n)
            two_n_lx = n * Two_Lx
            for l in np.arange(-order_y, order_y + 1, 1):
                bx2y2_abs_nl = bx2_abs_n * (by2 ** np.abs(l))
                two_l_ly = l * Two_Ly
                for m in np.arange(-order_z, order_z + 1, 1):
                    bx2y2z2_abs_nlm = bx2y2_abs_nl * (bz2 ** np.abs(m))
                    shift = np.array([two_n_lx, two_l_ly, m * Two_Lz])
                    for permu in np.arange(8):
                        xyz = shift - xx_yy_zz[:, permu]
                        delay = np.min(Fs_c * np.sqrt(np.sum((xyz - self.sensor_xyz) ** 2, axis=1)))
                        if delay <= H_length:
                            refl = (
                                bx1 ** np.abs(n - qq[permu])
                                * by1 ** np.abs(l - jj[permu])
                                * bz1 ** np.abs(m - kk[permu])
                                * bx2y2z2_abs_nlm
                            )
                            if np.sum(refl) < 1e-6:
                                continue
                            coordinates.append(xyz)
                            reflections.append(refl)
        if not coordinates:
            return np.zeros((3, 0)), np.zeros((len(self.F_abs), 0))
        return np.stack(coordinates, axis=1), np.stack(reflections, axis=1)

    def _render_sensor(self, H, sensor_no, isource_xyz, refl, atten_air, sensor_dir,
                       source_dir, tm_source, Fs_c, t, Two_pi_Tw, sampling_period,
                       N_frac, N_refl, Half_I_plusone, window, F_abs_N,
                       freq_grid, n_images):
        xyz = isource_xyz.T - self.sensor_xyz[sensor_no]
        dist = np.sqrt(np.sum(xyz ** 2, axis=1))
        b_refl = (refl.T / dist[:, None]) * (atten_air ** dist[:, None])
        half = np.empty((n_images, freq_grid.size))
        for idx in range(n_images):
            half[idx] = np.interp(freq_grid, F_abs_N, b_refl[idx])
        spectrum = np.concatenate((half, half[:, -2:0:-1]), axis=1)
        h_refl = np.real(np.fft.ifft(spectrum, n=N_refl, axis=1))
        h_refl = window * np.concatenate(
            (h_refl[:, Half_I_plusone - 1:N_refl], h_refl[:, :Half_I_plusone]),
            axis=1,
        )
        if n_images == 1:
            active = np.array([0])
        else:
            peaks = np.max(np.abs(h_refl[:, :Half_I_plusone]), axis=1)
            active = np.flatnonzero(peaks >= 1e-5)
        if active.size == 0:
            return

        delay = Fs_c * dist[active]
        rdelay = np.round(delay)
        frac = (delay - rdelay) * sampling_period
        t_Td = t[None, :] - frac[:, None]
        hsf = 0.5 * (1 + np.cos(Two_pi_Tw * t_Td)) * np.sinc(self.sampling_rate * t_Td)

        eps = np.finfo(float).eps
        xyz_active = xyz[active]
        xyz_source = xyz_active @ self.tm_sensor[sensor_no].T
        hyp = np.sqrt(xyz_source[:, 0] ** 2 + xyz_source[:, 1] ** 2)
        elevation = np.arctan(xyz_source[:, 2] / (hyp + eps))
        azimuth = np.arctan2(xyz_source[:, 1], xyz_source[:, 0])
        sensor_e = (np.round(elevation * 180 / np.pi) + 90).astype(int)
        sensor_a = (np.round(azimuth * 180 / np.pi) + 180).astype(int)

        xyz_sensor = -xyz_active @ tm_source.T
        hyp = np.sqrt(xyz_sensor[:, 0] ** 2 + xyz_sensor[:, 1] ** 2)
        elevation = np.arctan(xyz_sensor[:, 2] / (hyp + eps))
        azimuth = np.arctan2(xyz_sensor[:, 1], xyz_sensor[:, 0])
        source_e = (np.round(elevation * 180 / np.pi) + 90).astype(int)
        source_a = (np.round(azimuth * 180 / np.pi) + 180).astype(int)
        sensor_gain = sensor_dir[sensor_e, sensor_a]
        source_gain = source_dir[source_e, source_a]
        pad = np.zeros(N_frac)

        for row, image_index in enumerate(active):
            sig = np.concatenate((h_refl[image_index], pad))
            impulse = scipy_sig.lfilter(hsf[row], 1, sig)
            impulse = (impulse * sensor_gain[row]) * source_gain[row]
            len_h = len(impulse)
            adjust_delay = int(rdelay[row] - np.ceil(len_h / 2.0))
            start_index_Hp = max(adjust_delay + (adjust_delay >= 0), 0)
            stop_index_Hp = min(adjust_delay + len_h, H.shape[0])
            start_index_h = max(-adjust_delay, 0)
            stop_index_h = start_index_h + (stop_index_Hp - start_index_Hp)
            if stop_index_h < 0:
                continue
            H[start_index_Hp:stop_index_Hp, sensor_no] = (
                H[start_index_Hp:stop_index_Hp, sensor_no]
                + impulse[start_index_h:stop_index_h]
            )

    def __create_psi_theta_phi(self, source_off):
        c_psi = np.cos(np.pi/180*source_off[0])
        s_psi = np.sin(np.pi/180*source_off[0])
        c_theta = np.cos(-np.pi/180*source_off[1])
        s_theta = np.sin(-np.pi/180*source_off[1])
        c_phi = np.cos(np.pi/180*source_off[2])
        s_phi = np.sin(np.pi/180*source_off[2])
        return [c_psi, s_psi, c_theta, s_theta, c_phi, s_phi]

    def __create_tm(self, psi_theta_phi):
        c_psi, s_psi, c_theta, s_theta, c_phi, s_phi = psi_theta_phi
        tm_source = np.array([[c_theta*c_psi, \
                        c_theta*s_psi, \
                        -s_theta], \
                   [s_phi*s_theta*c_psi-c_phi*s_psi, \
                        s_phi*s_theta*s_psi+c_phi*c_psi, \
                        s_phi*c_theta], \
                   [c_phi*s_theta*c_psi+s_phi*s_psi, \
                        c_phi*s_theta*s_psi-s_phi*c_psi, \
                        c_phi*c_theta]])
        return tm_source