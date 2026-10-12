# GEOS626_seis

A public GitHub repository within the organization
[uafgeoteach](https://github.com/uafgeoteach). Contains materials for GEOS 626 Applied Seismology, a class at the University of Alaska Fairbanks by [Carl Tape](https://sites.google.com/alaska.edu/carltape/) ([ctape@alaska.edu](mailto:ctape@alaska.edu))

Course webpage: [GEOS 626](https://sites.google.com/alaska.edu/carltape/home/teaching/aseis)  

The repository can be obtained from GitHub with this command:
```
git clone --depth=1 https://github.com/uafgeoteach/GEOS626_seis.git
```

Several students and seismologists have contributed toward improving these materials; please see history.txt for details.

### Setup
---
A `.yml` file (see setup/ folder) lists dependencies. This file, executed within conda or docker, enables a user to establish the software tools needed to execute the iPython notebooks. (A dockerfile is also provided in setup/)

### How to run using Conda
---

- install conda (miniconda) if not done already
- if you already have miniconda, you might want to update to the latest version
- navigate to the setup folder
  ```bash
  cd GEOS626_seis/setup
  ```
- setup the conda environment
  ```bash
  conda env create -f seismo.yml
  ```
- activate the conda environment once the setup is complete
  ```bash
  conda activate seismo
  ```
- navigate back to the root of repository and launch jupyter
  ```bash
  cd ..
  jupyter lab
  ```
- browse and run notebooks as desired
