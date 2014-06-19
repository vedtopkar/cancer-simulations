# Cancer Lattice-Based Simulations
Authored by Dr. Bartek Waclaw, modified by Ved Topkar

## Compilation
Using the GNU C++ compiler (recommended):
`g++ bit_simulation.cpp bit_main.cpp functions.cpp -w -O3 -o executable`

(OR) Using GCC:
`gcc bit_simulation.cpp bit_main.cpp functions.cpp -w -O3 -o executable`

After running this command, you should have a binary file called `executable`.

## Running
To run the `executable` file:
`./executable outfolder 3 1`
This will create a folder called `outfolder`, run `3` independent runs of the program with a random seed of `1`.

## Parameters
For now, you can use different parameters by altering values in the `params.h` file.