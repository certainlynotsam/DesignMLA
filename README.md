# DesignMLA.m
A MATLAB script to help fabricate a microlens array with PowerPhotonic. 

## Overview

**DesignMLA.m** is a MATLAB script developed by Sam Daly (sgd46@cam.ac.uk) at the University of Cambridge, designed to fabricate four microlens array (MLA) patterns that are compatible with [PowerPhotonic LightForge](https://www.powerphotonic.com/products/lightforge/). 

## Current features
- Generate four distinct MLA patterns based on optical configurations.
- Compatibility with the LightForge system from PowerPhotonic.
- Flexible design parameters for different fabrication requirements.
- Packaged with MicrolensRelayCalculator.xlsx to help in the development of MLA optical configurations.

## Known issues
- When producing the full 15x15 layout (lines 213--216) the tiling must be changed to accomodate different size BFPs. 

## Requirements
- MATLAB (version 2022b or later)

## Acknowledgments
The foundations of this script is based on the contributions of Boya Zhang and Olivia Dovernor.

## Contributing
Contributions are welcome! Please contact Sam Daly (sgd46@cam.ac.uk).
