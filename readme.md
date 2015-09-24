This small program can be used to compute an approximation of the 2-points and 3-points correlation functions of point sets. It does so by using a double tree walk algorithm whose accuracy can be adjusted via an angular opening criterion.

I) Installing
-------------

Edit the `Makefile` and change options if needed, then type `make` and the binary will be created in the bin directory.

II) Usage
---------

The program is run by giving the name of a binary or ascii file containing a point set distribution in the NDfield format (see below):

    npcf.exe distribution.ND

by default, the auto-correlation function of that point set will be computed. The following options are available:

    -dist2 <filename>
     Used to specify the name of a second point set. If this option is used, the inter-correlation function of the two sets will be computed.

    -npoints N
     Used to specify the number of points for the correlation function. So far, only the 2-point and 3-point correlation function are available.
    
    -decimate V
     If this option is used, only a fraction V of the points in the distribution will be considered (0<V<=1).
     
    -periodic
     Use to force periodic boundary conditons

    -periodicity P
     Use to set periodic boundary conditions along specific axes only. `P=100` would enable periodic boundaries along dimension 0 (X-axis) of a 3D space.
     
    -theta T
     Cell angular opening criterion. The tree algorithm works by computing the grouped contribution of points clusters seen under a small angular size. The value of `T` sets the maximum angular size of a cell below which the contribution of the points it contains are grouped. T=0 implies an exact correlation function.

    -weight <filename> / -weight2 <filename>
     Used to specifie the name of a NDfield file containing the weights of each point in dist1 and dist2. If no weight is given, it is set to 1.

    -LS 
     Specify this option to force using the Landy-Szalay estimator: (DD-2.DR+RR)/RR.
     See Landy & Szalay, 1993, Apj (http://adsabs.harvard.edu/full/1993ApJ...412...64L)

    -outdir <DIR>
     Specify a directory to output the results
    
    -nbins <N>
     Number of bins to use in the distance dimension
     
    -nangbins <N>
     Number of bins to use in the angular dimensions (only for 3-point correlation functions)

    -linearbins
     Use that to force linear size for the bins (they are logarithmic by default)

    -threadlevel <L=3>
     The inter- or auto- correlation functions of cells at level L are given to individual threads. This number should be high enough to allow good work balance but not too high so that the system is not saturated !

III) File formats
-----------------

**ASCII format:**

A typical ASCII file looks like this for of set of 10000 3D points within a box of origin (0,0,0) and size (1,1,1):

    ANDFIELD COORDS    # Header
    3 10000            # 3 dimensions, 10000 points
    BBOX 0 0 0 1 1 1   # Bounding box 
    X1 Y1 Z1           # Points coordinates
    X2 Y2 Z2
    ...
    ...
    ...
    X9999 Y9999 Z9999

nb: The weights associated the point sets should be given as 1 dimensional set

**BINARY format:**
The binary format is organized as follows. 


|field       |type         |size   |comment|
|:-----------|:-----------:|:-----:|:------|
|dummy       |int(4B)      |1      |for FORTRAN compatibility                  |
|tag         |char(1B)     |16     |identifies the file type. Value : "NDFIELD"|
|dummy       |int(4B)      |1	   |                                           |
|dummy       |int(4B)      |1	   |                                           |
|ndims       |int(4B)      |1      |number of dimensions of the embedding space|
|dims        |int(4B)      |20     |size of the grid in pixels along each dimension, or [ndims,nparticles] if data represents particle coordinates (i.e. fdims_index=1)|
|fdims_index |int(4B)      |1      |0 if data represents a regular grid, 1 if it represents coordinates of tracer particles (always 1 in our case !)|
|datatype    |int(4B)      |1      |type of data stored (see below, should agree with Makefile)|
|x0          |double(8B)   |20     |origin of bounding box (first ndims val. are meaningfull)
|delta       |double(8B)   |20     |size of bounding box (first ndims val. are meaningfull)|
|dummy_ext   |char(1B)     |160    |dummy data reserved for future extensions  |
|dummy       |int(4B)      |1	   |                                           |
|dummy       |int(4B)      |1	   |                                           |
|data        |sizeof(type) |N      |data itself (the coordinates of particles, as `dims` times the number of particles values)|
|dummy       |int(4B)      |1	   |                                           |

with the possible 'datatypes' values being (64 bits system):

|type       |size(B) |numeric_type  |code       |
|-----------|:------:|:------------:|-----------|
|ND_CHAR    | 1      | integer      |1 (=1<<0)  |
|ND_UCHAR   | 1      | integer      |2 (=1<<1)  |
|ND_SHORT   | 2      | integer      |4 (=1<<2)  |
|ND_USHORT  | 2      | integer      |8 (=1<<3)  |
|ND_INT     | 4      | integer      |16 (=1<<4) |
|ND_UINT    | 4      | integer      |32 (=1<<5) |
|ND_LONG    | 8      | integer      |64 (=1<<6) |
|ND_ULONG   | 8      | integer      |128 (=1<<7)|
|ND_FLOAT   | 4      | float        |256 (=1<<8)|
|ND_DOUBLE  | 8      | float        |512 (=1<<9)|

nb1: blocks are delimited by dummy variables indicating the size of the blocks for FORTRAN compatibility, but they are ignored in C.
nb2: The weights associated the point sets should be given as 1 dimensional set