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

## Dependency Install Steps:
    ```
    sudo apt-get install libopencv-dev liblz4-dev libpcl-dev libdlib-dev libsuitesparse-dev libomp-dev cmake python python-pip
    ```


## Compliation:
Platforms: Linux

```
mkdir build
cd build
cmake ..
make -j4
```

## Usage
`./nonrigid`

