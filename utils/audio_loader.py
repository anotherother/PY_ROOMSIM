import soundfile as sf
import numpy as np
import samplerate

class AudioLoader:
    def __init__(self, standart_sr=16000):
        self.standart_sr = standart_sr

    def load_audio(self, path):
        """
        :param path:source path
        :return: samplerate, signal
        """
        signal, sr = sf.read(path,always_2d=True)
        if signal.ndim > 1:
            signal = np.mean(signal, axis=1)

        ratio = self.standart_sr / sr
        if sr > self.standart_sr:
            converter = 'sinc_best'
            reasample_signal = samplerate.resample(signal, ratio, converter)
        elif sr < self.standart_sr:
            assert False, "sample rate must be greater or equal then target_sr"
        else:
            reasample_signal = signal

        return int(sr * ratio), reasample_signal