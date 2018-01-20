# Facilitator_scripts_for_reactive_flow_through_fractured_rock_with_PHREEQC

![Preview](https://numericalenvironmental.files.wordpress.com/2018/01/secondarycusulfides.png?w=1234&h=1048)

This repository contains python scripts for pre- and post-processing PHREEQC simulations for flow through a dual-porosity column, emulating reactive transport through fractured porous media. The conceptualization entails a one-dimensional column of mobile-water cells, where advective-dispersive solute transport is simulated, and adjacent immobile-water cells that host only diffusive transport. The post-processing routine is designed to extend PHREEQC dual-porosity model results into a three-dimensional block of fractured rock, with fractures of arbitrary orientation. A brief introduction to both scripts is provided at (https://numericalenvironmental.wordpress.com/2018/01/20/pre-and-post-processing-of-a-dual-porosity-phreeqc-reactive-transport-simulation-with-python/), along with some example visualizations. Both are written on Python 2.7, but with the from_future module included to provide compatibility with Python 3. Pandas is a required module for the post-processor.

The pre-processor module uses the following text input files:

* grid_params.txt - parameters to populate TRANSPORT and MIX keyword blocks in PHREEQC, including numbers and lengths of mobile- and immobile-water cells, number of shifts (transfer of a pore volume between adjacent mobile-water cells), time step size per shift, aqueous diffusion coefficient, column dispersivity, and angular departure from vertical (used to compute implied fracture hydraulic conductivity under a unit hydraulic gradient).

* eq_phases_top.txt and eq_phases_bottom.txt - equilibrium phases and initial masses (units = mol/kgw in porous media) to be distributed among cells across the top portion of the model, which includes both mobile-water and associated immobile-water cells, and those for the bottom portion. Cell numbers indicate which mobile-water cells belong to which group; immobile-water cells are inferred. Multiple zones, not just “top” and “bottom,” can also be defined.

* kinetics_top.txt - equivalent to the eq_phases_top.txt file, but for kinetically-constrained species (as defined under user-supplied RATES keyword blocks). Multiple kinetics zones can also be defined.

The various output files from the pre-processor contain keyword blocks that can be pasted into a single PHREEQC input file, together with separate user-supplied model information. A complete example input file is included. Familiarity with PHREEQC is obviously assumed.
The post-processor module reads the input from one or more PHREEQC SELECTED_OUTPUT tabular simulation result files, converts them to pandas dataframes, extrudes and rotates the geometry, in then merges the results into a three-dimensional lattice. The lattice, which includes chemical and mineralogical concentration data for each point, can then be imported into a suitable visualization package (e.g., Golden Software’s Voxler, for example; http://www.goldensoftware.com/products/voxler). Required text input files include:

* frac_setup.txt - grid information for each fracture, one set of parameters per row, including gridding scheme used in corresponding PHREEQC dual-porosity reactive transport runs, parameters for extruding the fracture into a third dimension along the y-axis (i.e., ny and dy), a datum and x-axis intercept to match the rock block grid, fracture strike and dip, and the number of times the grid is to be subdivided and linearly interpolated along the z-axis to provide better spatial resolution of fracture details (nz_divide).

* grid_setup.txt - discretization constraints for the 3-D rock block containing the fractures

* PHREEQC SELECTED_OUTPUT files (e.g., dt2.txt and dt10.txt, as provided among the example files) - raw output files from PHREEQC dual-porosity run, renamed but otherwise unaltered

* background.txt - initial solid-phase concentrations of mineral or kinetic phases, or user-defined summary species (units = mol/kgw in porous media) in the background rock matrix, or non-fracture portion of the lattice domain

* ref_init.txt - initial solid-phase concentrations of mineral or kinetic phases, or user-defined summary species (units = mol/kgw in porous media) within a fracture zone; used to compute changes in mass over the course of a simulation, should such quantities be desired for visualization

* agglomeration list file (e.g., sec_Cu_sulfides.txt) - a file listing mineral or kinetic phases that will be grouped together under a single heading (a column in a dataframe); the prefix of the agglomeration list file will be used for the name

I'd appreciate hearing back from you if you find either script useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

