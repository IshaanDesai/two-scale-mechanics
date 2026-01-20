1. [__Configuration__](Configuration.md)
2. [Data Representation](Data_Format.md)
3. [Simulation Types](Simulation.md)
4. [Source Doc](SourceDoc.md)
   - [config.py](SourceDoc-Config.md)
   - [coupling.py](SourceDoc-Coupling.md)
   - [fnx.py](SourceDoc-Fnx.md)
   - [main.py](SourceDoc-Main.md)
   - [mesh_utils.py](SourceDoc-MeshUtils.md)
   - [meshes.py](SourceDoc-Meshes.md)
   - [simulation.py](SourceDoc-Simulation.md)
   - [util.py](SourceDoc-Util.md)

The input file format is regular JSON and on the highest level must have the following form:

```json
{
  "mesh": MESH_BLOCK,
  "problem": PROBLEM_BLOCK,
  "simulation": SIMULATION_BLOCK
}
```

For examples see the pre-configured files in the [examples](../examples) folder.
In the following descriptions, if no default is provided, the entry is not optional and must be provided.

### Description of all top level entries

| Name        | type | Description                             | Default |
|-------------|------|-----------------------------------------|---------|
| output_path | str  | path to output VTX file location        | mf_out  |
| mesh        | obj  | Block describing mesh                   | -       |
| problem     | obj  | Block describing problem parameters     | -       |
| simulation  | obj  | Block describing what simulation to use | -       |

### Mesh Block

| Name          | Type  | Description                                                        | Default | Options                                                                               |
|---------------|-------|--------------------------------------------------------------------|---------|---------------------------------------------------------------------------------------|
| path          | str   | Path to mesh file or name of shipped example                       | -       | path, BAR, NOTCH                                                                      |
| bc_dc_use_tag | bool  | Locate Dirichlet Boundary using mesh tags (=2)                     | -       | -                                                                                     |
| bc_dc_locator | str   | If no tag is used, then location is determined by locator function | -       | plane_xy_low, plane_xz_low, plane_yz_low, plane_xy_high, plane_xz_high, plane_yz_high |
| bc_dc_value   | float | Dirichlet Boundary Value                                           | 0.0     | -                                                                                     |
| bc_nm_use_tag | bool  | Locate Neumann Boundary using mesh tags (=3)                       | -       | -                                                                                     |
| bc_nm_locator | str   | If no tag is used, then location is determined by locator function | -       | plane_xy_low, plane_xz_low, plane_yz_low, plane_xy_high, plane_xz_high, plane_yz_high |
| bc_nm_value   | float | Neumann Boundary Value                                             | 0.0     | -                                                                                     |
| bc_nm_dim     | int   | Dimension in which nm bound should be applied                      | 1       | -                                                                                     |

### Problem Block

| Name        | Type  | Description                          | Default      | Options                    |
|-------------|-------|--------------------------------------|--------------|----------------------------|
| lambda      | float | TODO                                 | 10.0         | -                          |
| mu          | float | TODO                                 | 5.0          | -                          |
| alpha       | float | TODO                                 | 100.0        | -                          |
| strain_type | str   | Select between small or large strain | small_strain | small_strain, large_strain |
| elem_degree | int   | Element Degree                       | 2            | -                          |

### Simulation Block

| Name             | Type      | Description                                                         | Default                      | Options                               |
|------------------|-----------|---------------------------------------------------------------------|------------------------------|---------------------------------------|
| type             | str       | What simulation method to use                                       | -                            | MesoSim, PseudoCoupledSim, CoupledSim |
| input            | str       | Path to input h5 file if PseudoCoupledSim is used                   | - (only if PseudoCoupledSim) | -                                     |
| micro_type       | str       | What micro solver is being used in the coupled case                 | - (only if CoupledSim)       | ADA, PYFANS, NASMAT                   |
| write_state      | str       | Path to which internal meso state should be written to in h5 format | None                         | -                                     |
| write_state_type | list(str) | What data should be written                                         | None                         | E (eps), S (sig), U (displacement)    |
