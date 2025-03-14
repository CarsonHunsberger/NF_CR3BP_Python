# CR3BP Normal Form Code (Python)

This project contains all of the data an functions necessary for transforming between Cartesian CR3BP states, normal form states, and action-angle states for both the Birkhoff and resonant normal forms, at all three collinear libration points.

## Installation

To install, simply download this repository, extract it, and then copy the "NF_CR3BP" folder into your desired working directory.  
The normal form functions can then be loaded in by including the following line in any of your python scripts.

    from NF_CR3BP.NF_CR3BP import *

## Capabilities

1. Three different values of $\mu$
    - $\mu = 0.012154$ (Earth-Moon, used in some academic papers)
    - $\mu = 0.012150505856\dots$ (Also Earth-Moon)
    - $\mu = 3.0542\times 10^{-6}$ (Sun-Earth)
2. Two normal forms
    - Birkhoff
    - Resonant
3. Three libration points
    - $L_1$
    - $L_2$
    - $L_3$

All normal form approximations are of degree 11.

## Getting Started

Tutorial.ipybn has been included to provide an introduction to all of the functions included in this project.  
It is highly recommended to read through the comments and try running and fiddling with the examples.

For a higher-level view of the functionality offered, refer to the following diagram.

![Structure of normal form code](/NFStructure.png)


