## Resubmission
This is a resubmission. In this version, I have:

* Fixed undefined behavior sanitizer issues spotted by UBSAN.
* Added reference to published work related to the package in `DESCRIPTION`.

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

There was 2 NOTEs:

    * Maintainer: 'Aymeric Stamm <aymeric.stamm@math.cnrs.fr>'

      Days since last update: 2
      
      Possibly misspelled words in DESCRIPTION:
        Sangalli (37:5)
        Secchi (37:20)
        Vantini (37:31)
        Vitelli (37:43)

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

Possibly misspelled words are in fact names and the two URLs are in fact valid.

    * installed size is 10.6Mb
      sub-directories of 1Mb or more:
        data   2.7Mb
        libs   7.5Mb
