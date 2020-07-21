# Quality Control Dashboard for MaxQuant Proteomics Data

The Quality Control Dashboard for MaxQuant Proteomics Data is a R Shiny application, which takes the output text files (proteinGroups.txt, peptides.txt, msms.txt, summary.txt) from MaxQuant and a custom text file with the design of experiment information (DOE) for QC steps before further downstream analysis.

## Input files

### From MaxQuant output:

The files are located in the ~/combined/txt directory.

1. proteinGroups.txt (required)
2. peptides.txt (optional)
3. msms.txt (optional)
4. summary.txt (optional)

### Prepared by the user:

5. DOE (ex: doe.txt, optional but recommended)

 - The following 4 columns are required:        
   - Sample.id
   - Run.order
   - Type
   - Sample.type
 - Sample.id should match the data in MaxQuant.
 - Avoid naming the sample (Sample.id) starting with a number.
 - Only 'control' and 'sample' are admitted in the "Type" column.


## Run QC-MQ

### Method 1: Use QC-MQ online

QC-MQ is internally deployed at https://report.pri.bms.com/QC-MQ/ for online use.
Have the input files ready and follow the instructions from the dashboard to upload the files and at least input the **prefix** to run the QC.

### Method 2: Launch QC-MQ from R

#### Step 1: Install R and RStudio

- Please check CRAN (https://cran.r-project.org/) for the installation of R.

- Please check https://www.rstudio.com/ for the installation of RStudio.


#### Step 2: Install R Shiny package and other packages required by QC-MQ

- Once the R and RStudio are installed, start an R session using RStudio and run the following lines:

```r

load.lib <- c("shiny", "shinyBS", "shinydashboard", "shinyFiles", "colorspace", "tools", "shinyjs", "pvca", "Biobase", "ggthemes", "ggplot2", "scales", "shinyWidgets", "shinyalert", "reshape2", "stringr", "reshape2", "quantro", "plyr", "tidyr", "plotly", "dplyr")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)

```


#### Step 3: Download the source code to the working directory

- The zipped QC-MQ source code can be downloaded by clicking the **clone or download** button followed by clicking the **Download Zip** button from the page of this repository. 

- Move and unzip the downloaded file to the working directory.


#### Step 4: Start the app

- In the RStudio console, set the working directory at where the downloaded QC-MQ code is located by:

```r
setwd("~/path/workingdirectory")  # You have to put your own path instead of "~/path/workingdirectory"
```

- Launch the app by:

```r
library(shiny)
library(shinydashboard)
runApp(appDir = getwd())
```


#### Step 5: Run the QC

- Have the input files ready and follow the instructions from the dashboard to upload files and input **prefix** to run the QC for MaxQuant processed data


## Example files from MaxQuant output

Three different example datasets are included in the **QC-MQ/example_file_MQouput** folder. Users can upload the example files and set the prefix to test the app. (*note: msms.txt files are not included in this repository due to the size limitation, QC-MQ allows user running QC without uploading the msms.txt.*)

