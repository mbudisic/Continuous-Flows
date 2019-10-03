# Continuous Flows

This is a MATLAB object-oriented framework for continuous dynamical
systems. It provides a common interface for computing trajectories of velocity
fields, 2D Hamiltonian systems (stream-function based fluid flows), etc. The
intent of the package is to enable quickly swapping out the analyzed flow
while testing various algorithms for analysis of dynamical systems. It is also
intended to be easily extendable, so that adding a new flow to the library
requires minimal boilerplate code.

## How to install

To install the library, `git clone` the repository to a folder and then add
the directory *containing* `+ContinuousFlows` to MATLAB path. (Or click on "Download" link and unzip the downloaded file which creates the folder).

For example, if you clone the repository to
`/Users/mbudisic/Packages/ContinuousFlows` such that directory
`/Users/mbudisic/Packages/ContinuousFlows/+ContinuousFlows` exists, you can
add the package to MATLAB's path by
`addpath('/Users/mbudisic/Packages/ContinuousFlows')`.

After installing, you can check that the package is installed by running `doc
ContinuousFlows` which should give you a summary of the flows available in the
package.

## Demonstration

Open and run `demo.m` for demonstration of basic use of the package.
