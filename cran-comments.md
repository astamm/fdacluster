## Quick patch release


I am submitting a patch release days ago from last release because I did not
set up future workers correctly (i.e. with **fdacluster** loaded in them) which 
breaks reverse dependency **squat**. This is now fixed. I plan submitting a new
release of **squat** as soon as this patch release is accepted.

## Test environments
* local macOS R installation, R 4.4.2
* continuous integration via GH actions:
  * macOS latest release
  * windows latest release
  * ubuntu 20.04 latest release and devel
* [win-builder](https://win-builder.r-project.org/) (release, devel)
* [R-hub](https://builder.r-hub.io)
  - Windows Server 2022, R-devel, 64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There was no ERROR and no WARNINGs.

There was 1 NOTE:

    * checking installed package size ... NOTE
        installed size is  9.7Mb
        sub-directories of 1Mb or more:
          data   4.3Mb
          doc    2.6Mb
          help   2.0Mb
