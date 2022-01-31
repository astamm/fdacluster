## Test environments
* local macOS R installation, R 4.1.2
* continuous integration via GH actions:
  * macOS latest release
  * windows latest release
  * ubuntu 20.04 latest release and devel
* [win-builder](https://win-builder.r-project.org/) (release and devel)
* [R-hub](https://builder.r-hub.io)
  - Windows Server 2022, R-devel, 64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There was no ERROR and no WARNINGs.

There was 1 NOTE:

    * Maintainer: 'Aymeric Stamm <aymeric.stamm@math.cnrs.fr>'

      New submission

      Possibly mis-spelled words in DESCRIPTION:
        affine (33:45)

      Found the following (possibly) invalid URLs:
        URL: https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-8/issue-2/Analysis-of-AneuRisk65-data-k-mean-alignment/10.1214/14-EJS938A.full
          From: man/aneurisk65.Rd
                man/fdacluster.Rd
          Status: 500
          Message: Internal Server Error
        URL: https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-8/issue-2/AneuRisk65-A-dataset-of-three-dimensional-cerebral-vascular-geometries/10.1214/14-EJS938.full
          From: man/aneurisk65.Rd
          Status: 500
          Message: Internal Server Error

This is indeed a new submission. The word 'affine' is well defined and the two URLs are in fact valid.
