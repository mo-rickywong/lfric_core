.. -----------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

******************
Mesh Tools
******************
Programs primarily for mesh generation. Generated meshes are in NetCDF(.nc) files which aim to be compliant with UGRID (v1.0) conventions.
Additional mesh attributes are added within the UGRID convention to support informed application use cases not covered by the UGRID convention.

******************
Cubed-sphere
******************
Cubed-sphere mesh for applciations working with a global model.
Topology = Cube
Geometry = Spherical

******************
Planar mesh
******************
Planar mesh for applciations working with a regional model.
Topology = Grid
Geometry = * Curved surface (spherical region used for LAMs) implicitly used with lon-lat coordinate system
           * Flat surface (Flat region used for reasearch/idealised models)  implicitly used with cartesian coordinate system
