rm(list=ls())
# ------- SET WORKING DIRECTORY ------- #
if (! "rstudioapi" %in% installed.packages() ) install.packages("rstudioapi")
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# ------- PACKAGES ------- #
.packages = c("roxygen2")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load package
load <- lapply(.packages, require, character.only=TRUE)
if (all(unlist(load))) print("packages loaded succesfully")

remove.packages("TIE")
# delete current compiled version
unlink(list.files(pattern = ".tar.gz"))
unlink("TIE", recursive = TRUE, force = TRUE)



package.skeleton(name="TIE",
                 code_files = c("R/estim_functions.R", 
                                "R/test_functions.R",
                                "R/aux.R"))

unlink("TIE/man", recursive = TRUE)
dir.create("TIE/man")
dir.create("TIE/inst")
# copy bib file:
file.copy("inst/REFERENCES.bib", "TIE/inst/REFERENCES.bib")

# delete default readme
unlink("TIE/Read-and-delete-me")
# delete default DESCRIPTION FILE
unlink("TIE/DESCRIPTION")

# copy package description
file.copy("DESCRIPTION", "TIE/DESCRIPTION")


roxygen2::roxygenise(file.path(getwd(),"/TIE"))

system("R CMD build  TIE")
install.packages("TIE_1.0.0.tar.gz", repos = NULL, type = "source")

