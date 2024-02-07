.. -----------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

******************
Mesh Tools
******************
Programs for the generation of 2D mesh topologies. Output aims to be formatted as UGRID (v1.0) compliant NetCDF(.nc) files. Development of these tools are driven by application requirements on the LFRic Core infrastructure.

Generated meshes are two-dimensional, using quadrilateral cells arranged for a given use case. Additional mesh attributes not covered by the UGRID convention are added via NetCDF attributes in order to support specific use cases arising from application requirements.

******************
Cubed-Sphere
******************
The Cubed-Sphere mesh topology's primary use case is for Global General Circulation Models (GCMs) as an alternative to the regular lon-lat grid. Cells in this mesh topology are connected as a gridded-cube as the name suggests, though the geographical location of the the nodes are such that the mesh geometry is spherical.


******************
Planar mesh
******************
The planar mesh topology's use case is for Regional GCMs. The planar mesh topology is connected as a regular grid with nodes located on a:

Curved plane
    A gridded regional domain on a spherical surface used for Limited Area Models (LAMs). This mesh topology implicitly uses a spherical coordinate system (longitude-latitude) for the node locations. Full periodicity at the domain boundaries is not supported.

Flat plane
    A gridded regional domain on a flat plane, generally used for idealised modelling cases. This mesh topology implicity uses a cartesian coordinate system (xyz) for node locations. Full-periodicity at the domain boundaries is permitted.

