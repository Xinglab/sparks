sparks_installed <- require("SPARKS")
if (!sparks_installed) {
  install.packages("./", repos=NULL)
}
sparks_installed <- require("SPARKS")
if (!sparks_installed) {
  stop("could not install SPARKS")
}
