#source /opt/homebrew/bin/thisroot.sh
import glob
import uproot
import matplotlib.pyplot as plt
import numpy as np

from roentgen.absorption import Material
import astropy.units as u

thickness = 3 * u.mm

energy = u.Quantity(np.arange(5, 150, 0.1), 'keV')
mat = Material("Be", thickness)
trans = mat.transmission(energy)


angle_list = []
flux_list = []

for filename in glob.glob('../data/angle_run/DetectorHists_angleX*.root'):
    print(f"Processing {filename}")
    file = uproot.open(filename)
    this_angle = int(filename.split('angleX')[1].split('.')[0])

    data = file['h2'].to_numpy()
    plt.figure()
    plt.plot(data[1][:-1], data[0], '+', label='geant4')
    intensity = data[0][-30:].mean()
    plt.plot(energy.value, trans * intensity, label='roentgen 3mm Be total')
    plt.legend(loc='lower right')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Intensity [counts]')
    plt.savefig(filename + '.pdf')
    plt.close()

    angle_list.append(this_angle)
    flux_list.append(data[0].sum())

angles = np.array(angle_list)
fluxes = np.array(flux_list)

# normalize to maximum which should be angle = 0
fluxes = fluxes / fluxes.max()

sort = np.argsort(angles)
angles = angles[sort]
fluxes = fluxes[sort]

plt.figure()
plt.plot(angles, fluxes, 'o')
plt.plot(angles, fluxes)
plt.xlabel('source angle [deg]')
plt.ylabel('source intensity [normalized]')
plt.savefig('../data/angle_run/aperture_vs_angles.pdf')