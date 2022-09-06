# this is borrowed from GOstats in order to handle properly
# GOTERM and GOBPPARENTS and get rid of the "invisible binding" NOTES

GOenv <- function(what) {
    getAnnMap(what, "GO", load=TRUE,
              type=c("db", "env"))
}

.onLoad <- function(libname, pkgname) {
  options(Matrix.warnDeprecatedCoerce = 2)
}
