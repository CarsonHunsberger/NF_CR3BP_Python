# CR3BP Normal Form Code (Python)

This project contains all of the data an functions necessary for transforming between Cartesian CR3BP states, normal form states, and action-angle states for both the Birkhoff and resonant normal forms, at all three collinear libration points.

## Installation

This code is designed to run "out of the box", and has been verified to work on MATLAB versions 2022b and later.

To install, simply download this repository and extract it to wherever you'd like.

To access the NF_CR3BP functions from other folders, right click on the "NF_CR3BP_MATLAB" folder while in MATLAB and select "Add to path > Selected folders".


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

Tutorial.m has been included to provide an introduction to all of the functions included in this project.  
It is highly recommended to read through the comments and try running and fiddling with the examples.

For a higher-level view of the functionality offered, refer to the following diagram.

![Structure of normal form code](/NFStructure.png)

All of the functions have explanatory comments, meaning that typing "help \<function name\>" in MATLAB's command line will provide additional details.


