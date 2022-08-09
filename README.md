# GEOS626_seis

A public GitHub repository within the organization
[uafgeoteach](https://github.com/uafgeoteach). Contains materials for GEOS 626 Applied Seismology, a class at the University of Alaska Fairbanks by [Carl Tape](https://sites.google.com/alaska.edu/carltape/) ([ctape@alaska.edu](mailto:ctape@alaska.edu))

Course webpage: [GEOS 626](https://sites.google.com/alaska.edu/carltape/home/teaching/aseis)  

The repository can be obtained from GitHub with this command:
```
git clone --depth=1 https://github.com/uafgeoteach/GEOS626_seis.git
```

### Setup
---
A `.yml` file (see setup/ folder) lists dependencies. This file, executed within conda or docker, enables a user to establish the software tools needed to execute the iPython notebooks. (A dockerfile is also provided in setup/)

### How to Run Using Docker
---

First, go into setup directory by using following command:

```bash
# move into setup directory
cd setup

# (optional): check Makefile exists 
ls | grep Makefile
```


Once you are in the setup directory, run following command to build and run docker container.

``` bash
make
```

For more details, see `README.md` located in setup directory.