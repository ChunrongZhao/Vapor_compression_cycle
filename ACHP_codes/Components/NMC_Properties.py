# Cell_Properties.py
#
# Created: Oct. 2022, C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import json
from scipy.interpolate import LinearNDInterpolator
import os


# ----------------------------------------------------------------------
#  Battery cell properties
# ----------------------------------------------------------------------
class NMC_Properties:
    def __init__(self):
        self.tag = "LIB_NMC_18650"
        self.diameter = 0.018  # [m]
        self.height = 0.065  # [m]
        self.mass = 0.048  # [kg]
        self.surface_area = (np.pi * self.height * self.diameter) + (0.5 * np.pi * self.diameter ** 2)  # [m^2]
        self.volume = np.pi * (0.5 * self.diameter) ** 2 * self.height
        self.density = self.mass / self.volume  # [kg/m^3]
        self.electrode_area = 0.0342  # [m^2]

        self.max_voltage = 4.2  # [V]
        self.nominal_capacity = 3.55  # [Amp-Hrs]
        self.nominal_voltage = 3.6  # [V]
        self.charging_voltage = self.nominal_voltage  # [V]

        self.watt_hour_rating = self.nominal_capacity * self.nominal_voltage  # [Watt-hours]
        self.specific_energy = self.watt_hour_rating * 3600 / self.mass  # [J/kg]
        self.specific_power = self.specific_energy / self.nominal_capacity  # [W/kg]  need to be required
        self.resistance = 0.025  # [Ohms]

        self.specific_heat_capacity = 1108  # [J/kgK]
        self.radial_thermal_conductivity = 0.4  # [J/kgK]
        self.axial_thermal_conductivity = 32.2  # [J/kgK] # estimated

        return

    @staticmethod
    def compute_NMC_cell_state_variables(I, T, SOC):
        with open('NMC_Raw_Data.txt') as f:
            NMC_Raw_Data = f.read()

        battery_data = json.loads(NMC_Raw_Data)
        battery_data_voltage = battery_data['Voltage']

        # convert the data dictionary to the numpy array with the shape of (626, 4)
        i, j = 0, 0
        DATA = np.zeros((0, 4))
        for keys in battery_data_voltage.keys():
            x = battery_data_voltage[keys]
            for sub_keys in x.keys():
                # update current and temp here
                current, temp = i, j
                y = np.array(x[sub_keys])
                y[:, 0] = y[:, 0] / 3550.0
                z = np.zeros((y.shape[0], 2))
                z[:, 0] = current
                z[:, 1] = temp
                Z = np.concatenate((z, y), axis=1)
                DATA = np.concatenate((DATA, Z), axis=0)
                j += 10
            i += 2
            j = 0

        inputs = DATA[:, [0, 1, 2]]
        output = DATA[:, 3]
        Voltage_func = LinearNDInterpolator(inputs, output)
        xx = np.array([I, T, SOC])
        V_oc = Voltage_func(xx)

        return V_oc


if __name__ == "__main__":
    cwd = os.getcwd()
    files = os.listdir(cwd)
    print("Files in %r: %s" % (cwd, files))

    ospath = os.path.abspath('NMC_Raw_Data.txt')
    print(ospath)

    separator = os.path.sep
    x = os.path.dirname(ospath)
    rel_path = x + separator
    print(rel_path)
