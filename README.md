# Spin ED: A quick&dirty *parallel* Exact Diagonalization algorithm for interacting spins on a lattice. 

This is a beta version software.

A Lanczos based solver for interacting electrons on a lattice exploiting distributed memory MPI parallelization.  

### Dependencies

The code is written around the SciFortran library. Dependencies are:   

* gfortran > 9.0.0 **OR** ifort  > 13.0
* cmake > 3.0.0    
* [MPI](https://github.com/open-mpi/ompi)
* [SciFortran](https://github.com/QcmPlab/SciFortran)


### Installation
Installation is available using CMake.    

Clone the repo:

`git clone https://github.com/QcmPlab/spin_ed`

and from the just created directory make a standard out-of-source CMake compilation:

`mkdir build`  
 `cd build`  
`cmake ..`     
`make`     
`make install`   

The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/>` 
* `-DUSE_MPI=<yes>/no`  
* `-DVERBOSE=yes/<no> `  
* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  

### System loading: 

The library can be loaded into the operative system using one of the following, automatically generated, methods:    

1. environment module file 
2. homebrew `bash` script 
3. pkg-config file 


Method 2. has the advantage of making `uninstall` operation feasible. 
Information about loading the library are printed on the screen at the end of the installation process.

### Uninstall

The library is removed with the command:

`make uninstall`

from the same building directory as for the installation part. 



For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***COPYRIGHT & LICENSING***  
Copyright  (c), Adriano Amaricci.  
All rights reserved. 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

--



