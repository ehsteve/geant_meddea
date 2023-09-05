import numpy as np
from datetime import datetime
import os
import subprocess
import shutil
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import os

# clean up
if os.path.exists('../data/angle_run'):
    shutil.rmtree('../data/angle_run/')
os.mkdir('../data/angle_run/')

#angle_arr = np.sort(np.concatenate([np.arange(-45, 45, 5), np.arange(-4, 5, 1)]))
#angle_arr = np.linspace(-1, 1, 3, dtype='int')
angle_arr1 = np.arange(-10, 11, 1)  # degrees
angle_arr2 = np.array([0, 10, 15, 20, 25, 30])
angle_arr = np.sort(np.concatenate([angle_arr1, angle_arr2]))

print(angle_arr)

tstart = datetime.now()

num_beam = 1000000

orig_beam_direction = [0, 0, -1]
rot_axis = [1, 0, 0]
rot_axis = rot_axis / norm(rot_axis)
theta = np.deg2rad(45)
beam_z_position = 10. # cm


for this_angle in angle_arr:
    print(f"Run angle = {this_angle}")
    # create new macro file
    shutil.copy('../build/run_angles_base.mac', '../build/this_angle_run.mac')
    theta = np.deg2rad(this_angle)
    rot_matrix = Rotation.from_rotvec(theta * rot_axis)
    new_beam_direction = rot_matrix.apply(orig_beam_direction)
    print(new_beam_direction)
    # reposition the beam to keep illumination centered on detector
    with open('../build/this_angle_run.mac', 'a') as fp:
        fp.write(f'/gps/direction {new_beam_direction[0]:.4f} {new_beam_direction[1]:.4f} {new_beam_direction[2]:.4f}\n')
        fp.write(f'/gps/pos/centre 0. {-beam_z_position * np.sin(theta):.2f} {beam_z_position:.2f} cm\n')
        fp.write(f'/run/beamOn {num_beam}')

    # run geant4 sim
    subprocess.run(["../build/meddea_sim", "../build/this_angle_run.mac"])
    # rename output root file and move it into new directory
    new_filename = f"../data/angle_run/DetectorHists_angleX{this_angle}.root"
    os.rename("DetectorHists.root", new_filename)

tend = datetime.now()
dt = tend - tstart
print(f"Total run time: {dt}")
print(f"Average time per run: {dt/len(angle_arr)}")