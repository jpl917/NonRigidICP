## Non-Rigid ICP

An implementation of non-rigid ICP to register the scanned model to the unified template.

## Dependencies
The code are based on the following third parties.
- Eigen
- SuiteSparse
- Boost
- Opencv
- dlib
- Trimesh2


## Compliation:
Platforms: Linux

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make -j4`

## Usage
`./nonrigid`


## Input folder include:
1. **image/**		the folder contains all the images  
2. **P/**		        the folder contains the projection matrix for each camera  
3. **T/**			the folder contains the transformation matrix for each camera  
4. **photoscan.ply**       the generated dense mesh model   

## Output
1. **nonrigid_landmarks/**  the intermediate 2D face landmark detection results   
2. **nonrigid_results/**	the non-rigid ICP results  
