# madeal: Masoud Ahmadi's deal.II codes

![Project Logo](logo.png)

## Project Overview

Utilizing **deal.II** finite element library, **madeal** is a collection of codes developed by Masoud Ahmadi as part of his PhD research at the University of Glasgow.  Initiated in 2022, this project focuses on deal.II's capabilities for analyzing large deformations of elastomeric composites reinforced by randomly dispersed fibres.

## Features

- **Large Elastic Deformation**: analysis large deformations of hyperelastic materials.
- **Incompressibility**: handles nearly-incompressible materials using mixed three-field formulations.
- **Plane-stress approximation**: facilitates nonlinear FEM under plane-stress conditions.
- **Canonical Boundary Conditions**: implements both affine and periodic boundary. conditions
- **2D RVE Generation**: utilizes pixel meshing technique to generate 2D representative volume elements (RVEs).
- **Benchmarking**: Includes several examples to validate and benchmark the models.

## Installation

Ensure the system has a C++ compiler (e.g., GCC, Clang) and CMake installed. 
To install **madeal**, follow these steps:

1. **Install deal.II**: Install deal.II: The deal.II library must be pre-installed on the system.
    Follow the installation guide available at [deal.II official website](https://www.dealii.org/).
2. **Clone the Repository**: Clone this repository into the deal.II directory.
```sh
git clone https://github.com/Masoud16ahm/madeal.git
```

## Templates

**Templates**  for various physical scenarios are provided. For example, to solve a small deformation elastic problem, use "template1". Modify this template and the inputs according to your specific problem requirements.


## Examples

Navigate to the **examples** folder to view and run benchmark examples demonstrating madeal's capabilities. To compile and run an example:
```sh
cmake .
make run
```

## License

**madeal** is distributed under the **MIT** License, which grants users specific rights to use, modify, and distribute the package both in original and modified form under certain conditions:
1. **Usage**: Users are allowed to use the package for any purpose, including commercial and academic.
2. **Modification**: Users may modify the source code, ensuring that modifications are documented.
3. **Distribution**: Redistribution of the original or modified software is permitted as long as the copies retain the original license terms and a notice of any modifications.
4. **Attribution**: Users must credit the original authors in any publications or software distributions that utilize **madeal**.

## Testing

Testing is implemented using **Google Test**. To run the tests:
```sh
cmake .
make
./gtest
```

## Contact

For inquiries and additional information, contact Masoud Ahmadi at:

- Emails: Masoud.Ahmadi@glasgow.ac.uk and Masoud_Ahmadi_pr@yahoo.com

