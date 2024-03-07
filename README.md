How to use genotype data to estimate the census size of a wild population.  Bill Sherwin EERC, BEES, UNSW-Sydney   W.Sherwin@unsw.edu.au

BACKGROUND
- First read Sherwin, WB 2024. DnaDot - Fixing Ecology and Evolution's blind spot, Population size. Ecological Indicators {****** complete ref when published}
- Pay especial attention to "4.1 Calculation Protocol for Wild Populations"
- A number of files are provided, for users with different interests and platforms.  

RSEARCH PROGAM in MATLAB: This is the program used in the research, plus the function that it calls to calculate expectations:
- Extensively commented pdfs, human-readable without a MATLAB licence: NcHyper231006RandPa1a2.pdf; HypExactExp230607.pdf  
- MATLAB code: NcHyper231006RandPa1a2.m; HypExactExp230607.m    
- If you do not have a MATLAB licence, see alternatives below.

POINT-AND-CLICK FOR USERS: This is probably the best option for most users; it is a point-and-click application which allows users to analyse their own genetic data in an EXCEL file, without having a MATLAB licence. You need four files in one folder. 
- Prepare your input of scenario choices and genotypic data in '.xls' files, following the style of the two input files from this GitHub ('Choices.xls' and 'Genotypes.xls').  
- "Choices.xls": Cell A1 is jj, the proportion of that sample to be included in each jacknife subsample; recommended values are 0.5 to 0.9; Cell A2 is minNtry, the minimum census size to hypothesise; eleven census sizes will be hypothesised, from minNtry to maxNtry=3*minNtry
-"Genotype.xls"  This contains only the genotype data for your sample, ie no row or column headings, no sex column, etc. Autosomal diploid data only, no missing data.  SNPs or Microsatellites are OK, or anything for which there are alleles with numerical codes eg '185' or '001' (with no quote marks); codes can be re-used in different loci, eg reference allele always '001', alternative allele always '002'. Each locus in two adjacent columns for the two diploid alleles.  Each individual on one row (order randomised - see how below). 
- An aside, how to randomise an EXCEL file:(1) Construct the file 'Genotypes.xls' as described in the previous point (2) Make a new column and fill every cell with     =RAND()      which will insert random numbers (3) Highlight all columns plus the new column (4) Go to Data>Sort (5) Choose the new column of random numbers as the basis for the sort (6) Then remove the random number column and your file is ready to use, so save it
- Put the two input files from this GitHub into the same folder as 'NcHyper240129.exe' and 'HypExactExp230607.exe', both also from this GitHub
- Double Click on 'NcHyper240129.exe'. 
- You should not get the error "Could not find MATLAB runtime", but if you do, then go to https://www.mathworks.com/products/compiler/matlab-runtime.html   and install the free version R2023b (23.2) into the same folder as the files from this GitHub; then try the above procedure again.
- A file called 'CensusOutput.csv' should appear in the same folder, it will be a comma-separated file (can be read by EXCEL)
- Results in this file are minNtry (the smallest hypothesised census size), maxNtry (the largest hypothesised census size), jj(proprtion of the total sample in each jackknife subsample), Average estimate of census size (Nest) over all loci, and standard deviation SD and standard error SE of Nest.
- NB if you run the sample input files, you will get zero SD and SE, because the sample dataset is so tiny; you should get nonzero SE and SD when you use your own data
- If you do multiple runs of the program it is recommended that you immediately rename each output file to tell you what was special about that run

OCTAVE FOR USERS: There is a translation of the same code into OCTAVE, which can be run without paying for a licence, and the OCTAVE console will allow users to read and modify the code if they wish.  
- Install GNU Octave (https://octave.org/) ; it is free.  In Windows, it will put a button called 'GNU Octave' on your desktop
- From this GitHub, put into a single folder on your computer: "HypExactExp240227Octave.m','NcHyper240229Octave.m','Choices.txt','Genotypes.txt'
- Modify the two input files to define your scenario and data, as explained above for the click-and-point app  
- Click on the 'GNU Octave' button to open it
- In the GNU Octave console, go to the directory containing the input files and the code
- Open 'NcHyper240229Octave.m', then click on the 'run' button on the console (green triangle, mid-top)
- In the same folder, you should get a new output file like the one from the click-and-point app.

COMMAND-LINE FOR USERS: Fourthly, for users who wish to incorporate the code into a larger pipeline that includes other software of their choice, there is a BASH script that calls and runs the OCTAVE code. 
- Install OCTAVE as described above
- You need five files in one folder. 
- From this GitHub, the BASH code 'nchyp-mcfarlane-github.sh'
- From this GitHub: "HypExactExp240227Octave.m','NcHyper240229Octave.m','Choices.txt','Genotypes.txt'
- Modify your input of scenario choices and genotypic data as for Octave above.  
- run the BASH script
- You should get '.txt' output as described for OCTAVE above

NEAR FUTURE: A molecular ecology organisation which has professional programmers and user support is building an R-version of the DnaDot program with greater flexibility.  The DnaDot program will soon be available as an R function ("gl.DnaDot") within the R package "dartR.popgen" which is available in CRAN and can be installed using the R command: install.packages("dartR.popgen"), and of course will have usage instructions after testing on users.
