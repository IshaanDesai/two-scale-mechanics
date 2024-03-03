"""
Dummy meso solver which writes random strain values to and reads stresses from preCICE.
"""
from __future__ import division
import numpy as np
import precice

mesh_name = 'dummy-mesh'

# RVE_ID = Unique ID for each macro gauss point [INTEGER]
# MOD_ID = ID for the micromechanical model [INTEGER]
# RUC_SIZE = Problem size of the RUC model [INTEGER]
write_data_names = ['RVE_ID', 'MOD_ID', 'RUC_SIZE', 'strains1to3', 'strains4to6']
read_data_names = ['stresses1to3', 'stresses4to6']

num_vertices = 10  # Number of vertices

participant = precice.Participant('Meso-scale-dummy', '../precice-config-nasmat-scaling.xml', 0, 1)

assert (participant.requires_mesh_connectivity_for(mesh_name) is False)

vertices = np.zeros((num_vertices, participant.get_mesh_dimensions(mesh_name)))

read_data = dict()
for name in read_data_names:
    read_data[name] = np.zeros((num_vertices, participant.get_data_dimensions(mesh_name, name)))

write_data = dict()
for name in write_data_names:
    write_data[name] = np.zeros((num_vertices, participant.get_data_dimensions(mesh_name, name)))

for x in range(num_vertices):
    for y in range(participant.get_mesh_dimensions(mesh_name)):
        vertices[x, y] = x

for name in read_data_names:
    for x in range(num_vertices):
        for y in range(participant.get_data_dimensions(mesh_name, name)):
            read_data[name][x, y] = x

for name in write_data_names:
    for x in range(num_vertices):
        for y in range(participant.get_data_dimensions(mesh_name, name)):
            write_data[name][x, y] = x

vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

participant.initialize()

while participant.is_coupling_ongoing():
    if participant.requires_writing_checkpoint():
        print("DUMMY: Writing iteration checkpoint")

    dt = participant.get_max_time_step_size()

    for name in read_data_names:
        read_data[name] = participant.read_data(mesh_name, name, vertex_ids, dt)

    for name in write_data_names:
        participant.write_data(mesh_name, name, vertex_ids, write_data[name])

    print("DUMMY: Advancing in time")
    participant.advance(dt)

    if participant.requires_reading_checkpoint():
        print("DUMMY: Reading iteration checkpoint")

participant.finalize()
print("DUMMY: Closing meso solver dummy...")
