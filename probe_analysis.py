from char_compare import Manager
from char_compare import data_slicer as ds
from char_compare import langmuir_models2 as lm
from Fluctuation_analysis import lp_spectrum as ls
import numpy as np
import matplotlib.pyplot as plt



def create_series(iprobe):
        iprobe.voltage *= 100
        iprobe.current *= -2
        vmcs = ds.VmaxClusterSlicer(iprobe, offset = -40, half_periods = 2)
        model = lm.WeinlichSheathExpansion
        iprobe.fit(model, vmcs)
        
        Te = iprobe.fit_result(model, vmcs, 'params', 'Te')
        ne = iprobe.fit_result(model, vmcs, 'params', 'ne')
        Vf = iprobe.fit_result(model, vmcs, 'params', 'Vf')
        time = iprobe.fit_result(model, vmcs, 'time')
        Vp = Vf + 3.0*(Te)
        
        return (time, ne, Te, Vf, Vp)