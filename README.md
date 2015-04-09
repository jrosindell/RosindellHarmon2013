# RosindellHarmon2013

Code of the article: 

 * Rosindell, J., Harmon, L. J. (2013), A unified model of species immigration, extinction and abundance on islands. Journal of Biogeography, 40: 1107–1118. doi: 10.1111/jbi.12064

I made it to compile and run under Qt Creator.

Below is the PDF of the original authors converted to Markdown:

## readme.pdf

This appendix includes the following files licensed under an MIT license.
Full license information is available in the files themselves. 
Please refer to the paper when using the software and e­‐mail any enquiries to james@rosindell.org 

### Dynamics_BI_1.cpp
Edit, compile and run this code for simulations with basic immigration.

### Dynamics_CI_1.cpp
Edit, compile and run this code for simulations with clustered immigration.

### Dynamics_PI_1.cpp
Edit, compile and run this code for simulations with protracted immigration.

### metacommunity_DLS3.txt
### metacommunity_LS4.txt
These files give the metacommunity species abundances.You can change these
for your own files or add your own files (editing the code appropriately).The
code will not run unless these files or replacements are available in the same
folder as the executable file. 

### RandomWrap.cpp 
This is needed to compile the code - it includes a basic random number
generator.It is recommended you read the comments and follow the
instructions to upgrade the random number generator to a better one. 
Unfortunately the generator we use was not written by us so be cannot provide a
copy under an MIT open source license.
