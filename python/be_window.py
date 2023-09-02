import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization import quantity_support
quantity_support()

from roentgen.absorption import Material

thickness = 5 * u.mm

material_list = ['Be', 'Al']
energy = u.Quantity(np.arange(5, 1000, 0.01), 'keV')

for material in material_list:
    mat = Material(material, thickness)
    plt.plot(energy, mat.transmission(energy), label=mat.name)

plt.xlim(0, 100)
plt.ylim(0, 1)
plt.xlabel('Energy [' + str(energy.unit) + ']')
plt.ylabel('Transmission')
plt.legend(loc='lower right')
#plt.show()
plt.savefig("be_window_transmission.png")