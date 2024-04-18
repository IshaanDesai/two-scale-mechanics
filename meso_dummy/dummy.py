"""
Dummy meso solver which writes random strain values to and reads stresses from preCICE.
"""
from __future__ import division
import numpy as np
import precice

mesh_name = 'dummy-mesh'
mod_id = 103
ruc_size = 64
num_vertices = 128  # Number of vertices

# Data names
write_data_names = ['rve_id', 'mod_id', 'ruc_size', 'strains1to3', 'strains4to6']
read_data_names = ['stresses1to3', 'stresses4to6','cmat1', 'cmat2', 'cmat3', 'cmat4', 'cmat5', 'cmat6', 'cmat7', 'conv']

# Creating precice participant
participant = precice.Participant('Meso-scale-dummy', '../precice-config-nasmat-scaling.xml', 0, 1)

# Initializing vertices
vertices = np.zeros((num_vertices, participant.get_mesh_dimensions(mesh_name)))
for x in range(num_vertices):
    vertices[x, 0] = float(x)/float(num_vertices)
vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

# Initializing write data
write_data = dict()
for name in write_data_names:
    write_data[name] = np.zeros((num_vertices, participant.get_data_dimensions(mesh_name, name)))

for x in range(num_vertices):
    write_data['rve_id'][x, 0] = x+1
    write_data['mod_id'][x, 0] = mod_id
    write_data['ruc_size'][x, 0] = ruc_size
    write_data['strains1to3'][x, :] = np.random.uniform(0, 1e-4, 3) # Random strain values
    write_data['strains4to6'][x, :] = np.random.uniform(0, 1e-6, 3) # Random strain values

# Initializing read data
read_data = dict()
for name in read_data_names:
    read_data[name] = np.zeros((num_vertices, participant.get_data_dimensions(mesh_name, name)))

participant.initialize()

print("Starting meso solver dummy...")

while participant.is_coupling_ongoing():
    dt = participant.get_max_time_step_size()

    for name in write_data_names:
        participant.write_data(mesh_name, name, vertex_ids, write_data[name])

    participant.advance(dt)

    if (participant.get_max_time_step_size()>0.0):
        for name in read_data_names:
            read_data[name] = participant.read_data(mesh_name, name, vertex_ids, dt)

        print("Read data at the end of the time step")
        #print(read_data['cmat6'])

participant.finalize()
print("Closing meso solver dummy...")
