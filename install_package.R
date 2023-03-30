sparks_installed <- require("SPARKS")
if (!sparks_installed) {
  # install.packages("./", repos=NULL)
  devtools::install_github('Xinglab/sparks')
}
sparks_installed <- require("SPARKS")
if (!sparks_installed) {
  stop("could not install SPARKS")
}
