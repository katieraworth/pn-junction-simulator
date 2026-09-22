import numpy as np
from constants import *

def built_in_voltage(NA, ND):
    """Built-in voltage (V) of a PN junction.
    NA = acceptor doping on the p-side (cm^-3)
    ND = donor doping on the n-side (cm^-3)
    """
    return Vt * np.log(NA * ND / ni**2)


if __name__ == "__main__":
    Vbi = built_in_voltage(1e17, 1e16)
    print(f"Built-in voltage: {Vbi:.3f} V")