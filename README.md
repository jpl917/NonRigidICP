## Non-Rigid ICP

An implementation of non-rigid ICP to register the scanned model to the unified template.

## Dependencies
Platforms: Ubuntu

The code are based on the following third parties.
- Eigen
- SuiteSparse
- Boost
- Opencv
- dlib
- Trimesh2 (included in 3rdparty)

## Install Steps:
- 1. run blow command
```
sudo apt-get install libopencv-dev liblz4-dev libpcl-dev libdlib-dev libsuitesparse-dev libomp-dev cmake python python-pip
```

- 2.`sudo vim /usr/include/dlib/gui_core/gui_core_kernel_2.h`  and comment line 11 and 12.<br>
```
//#error "DLIB_NO_GUI_SUPPORT is defined so you can't use the GUI code.  Turn DLIB_NO_GUI_SUPPORT off if you want to use it."
//#error "Also make sure you have libx11-dev installed on your system"
```

- 3. recomplie trimesh
```
cd 3rdparty/trimesh2/
make clean && make
```

- 4. compliation:
```
mkdir build
cd build
cmake ..
make -j4
```

## Usage
`./nonrigid`

