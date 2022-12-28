# exp_thermal_conductivity

The MATLAB codes for reproducing the figures and tables in our paper, *Thermal conductivity reconstruction method with application in a face milling  operation*, by Everton Boos, Fermín S.V. Bazán and Vanda M. Luchesi. 

Go to:

- [Example_1](/Example_1), for Example 1 codes
- [Example_2](/Example_2), for Example 2 codes
- [Example_section_4](/Example_section_4), for Section 4 examples (synthetic and experimental data)

Notice the use of the function **get_l**, which is part of the Regularization Tools package developed by Dr. Per Christian Hansen: http://www2.compute.dtu.dk/~pcha/Regutools/

## How to use

Main codes:
- **exp_direct.m**: runs the direct problem solver
- **exp_inverse.m** (and its variants **exp_inverse_over.m**, **exp_inverse_experimental.m** and **exp_inverse_synthetic.m**): runs the inverse problem solver for each specific scenario, used to build data in tables and figures

To see data used in figures and tables:
- **repetition_results.m** (and its variants **repetition_results_experimental.m** and **repetition_results_synthetic.m**): processes data in the .mat files (which were previously collected in the inverse problems solvers) into figures and tables
- **plot_line_average.m**: shows line plots present in Example 1 and 2 from data in the .mat files
- **plot_domain.m**: plots the domain for MS1 and MS2 measurement strategies present in Figure 2

The remaining .m and .txt files are auxiliary. All functions contain descriptions and comments in order to ease the reader.

## Licence Information

This library is free software. It can be redistributed and/or modified under the terms of the MIT License.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright (c) 2022 by Everton Boos
