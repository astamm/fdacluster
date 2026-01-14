## Test environments
* local macOS R installation, R 4.5.2
* continuous integration via GH actions:
  * macOS latest release
  * windows latest release
  * ubuntu 22.04 latest release and devel
* [win-builder](https://win-builder.r-project.org/) (release, devel)
* [R-hub](https://builder.r-hub.io): all platforms passed except some which missed C/C++ libraries to install nloptr.

## R CMD check results
There was no ERROR and no WARNINGs.

There was 1 NOTE:

    * checking installed package size ... NOTE
        installed size is  9.7Mb
        sub-directories of 1Mb or more:
          data   4.3Mb
          doc    2.6Mb
          help   2.0Mb
