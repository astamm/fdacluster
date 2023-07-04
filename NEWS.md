# fdacluster 0.3.0

* Added median centroid type;
* Median and mean centroid types are now defined on the union of individual grids;
* Simplified `caps` class to avoid storing objects multiple times under different names;
* Added vignette on initialization strategies for k-means;
* Added article on use case about the Berkeley growth study;
* Added article on supported input formats.

# fdacluster 0.2.2

* Make sure one can use **fdacluster** with namespace notation.
* Make sure not to use **fda** or **funData** before checking it is available.

# fdacluster 0.2.1

* Add DBSCAN clustering;
* Fix C++ compiler issues that errored when accessing empty vectors.

# fdacluster 0.2.0

* Add hierarchical clustering;
* Enforce `n_clusters` in output via linear programming (LP) using the 
**lpSolve** package;
* New [`caps`](https://astamm.github.io/fdacluster/reference/caps.html) class 
for storing results from functional **C**lustering with **A**mplitude and 
**P**hase **S**eparation in a consistent way;
* Add tools for comparing clustering results (`mcaps` objects, `autoplot` and 
`plot` specialized method implementations);
* Add seeding strategies for kmeans (via hierarchical clustering or k-means++ or 
k-means++ with exhaustive search of the first center or exhaustive search of all 
the centers);
* Add within-cluster domain auto-extension via mean imputation;
* Add possibility to cluster according to phase variability instead of amplitude 
variability.
* Renaming of functions: to perform k-means with alignment, now use 
[`fdakmeans()`](https://astamm.github.io/fdacluster/reference/fdakmeans.html), 
to perform HAC with alignment, now use 
[`fdahclust()`](https://astamm.github.io/fdacluster/reference/fdahclust.html).

# fdacluster 0.1.1

* Fixed undefined behavior sanitizer issues spotted by UBSAN.
* Added reference to published work related to the package in `DESCRIPTION`.

# fdacluster 0.1.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
