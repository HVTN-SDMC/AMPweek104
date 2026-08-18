# Analysis of HIV-1 acquisition in the Antibody Mediated Prevention (AMP) trials through the Week 104 final study visit.

R code implementing data analyses that generated figures and tables in the manuscript

deCamp et al., Extended in vivo activity of VRC01 HIV-1 monoprophylaxis in the AMP 
trials evidenced by initial suppression and breakthrough resistance


## 1. System requirements and installation guide

Unless otherwise specified, the software was tested on macOS 26.5.2 running 
R version 4.5.1. R packages required for analyses are listed in the renv.lock 
file. To restore all the R packages specified in the renv.lock file, follow the 
instructions below:
  
  * Install required version of R.
* Clone or download the repository to your local machine and open it as an 
RStudio project.
* Install R package "renv".
* Run renv::restore(). It creates the folder renv/library/ and installs 
packages or links recorded in renv.lock file into the project library.

## 2. Data access

The datasets required to run these analyses are hosted on the Harvard Dataverse 
and are available at:

https://doi.org/10.7910/DVN/WEJPUF

To reproduce the analyses in this repository:

* Download all data files from the Dataverse repository at the DOI above.
* Place the downloaded files into the `data/` directory of this repository.
* File names should match those referenced in the analysis scripts (see `code/`).

Data are made available under the terms specified on the Dataverse page. Please 
refer to the Dataverse listing for the data use agreement, citation requirements, 
and any access restrictions that may apply.

## 3. Generating derived datasets

After downloading the raw data files into `data/` (see Section 2), run the 
following scripts to generate the derived datasets used by the figure and 
table scripts. Scripts can be run in any order, and each writes its output 
back into the `data/` directory:

* `code/mk_data1.R`
* `code/mk_data2.R`
* `code/mk_data3.R`
* `code/mk_data4.R`
* `code/mk_data5.R`
* `code/mk_data6.R`
* `code/mk_data7.R`

## 4. Running the analyses

Once the derived datasets have been generated (see Section 3), each figure 
and table in the manuscript can be reproduced by running its corresponding 
script in `code/`. Scripts are named according to the figure or table they 
produce:

* `figX.R` / `tableX.R` — main text figures and tables, where `X` is the 
  figure or table number (e.g. `fig4.R`, `table1a.R`).
* `figSX.R` / `tableSX.R` — supplemental figures and tables, using the same 
  numbering convention (e.g. `figS7.R`, `figS29_S30.R`).
* Suffixes such as `a`, `b`, `c` indicate a sub-panel of a figure or a subset 
  of a table (e.g. `fig2ab.R` produces panels 2a and 2b; `table1a.R` through 
  `table1h.R` each produce a portion of Table 1).

Each script can be run independently — there is no required order — and 
saves its output to `output/figures/` or `output/tables/`. Note that the 
output file name does not always match the script's file name; consult the 
script itself (e.g. the file name(s) passed to `ggsave()`, `write.csv()`, 
or similar output functions) to identify the corresponding output file.


