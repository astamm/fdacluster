.onUnload <- function(libpath) {
  library.dynam.unload("fdacluster", libpath)
}
