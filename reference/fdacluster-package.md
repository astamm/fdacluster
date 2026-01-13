# fdacluster: Joint Clustering and Alignment of Functional Data

Implementations of the k-means, hierarchical agglomerative and DBSCAN
clustering methods for functional data which allows for jointly aligning
and clustering curves. It supports functional data defined on
one-dimensional domains but possibly evaluating in multivariate
codomains. It supports functional data defined in arrays but also via
the 'fd' and 'funData' classes for functional data defined in the 'fda'
and 'funData' packages respectively. It currently supports shift,
dilation and affine warping functions for functional data defined on the
real line and uses the SRVF framework to handle boundary-preserving
warping for functional data defined on a specific interval. Main
reference for the k-means algorithm: Sangalli L.M., Secchi P., Vantini
S., Vitelli V. (2010) "k-mean alignment for curve clustering"
[doi:10.1016/j.csda.2009.12.008](https://doi.org/10.1016/j.csda.2009.12.008)
. Main reference for the SRVF framework: Tucker, J. D., Wu, W., &
Srivastava, A. (2013) "Generative models for functional data using phase
and amplitude separation"
[doi:10.1016/j.csda.2012.12.001](https://doi.org/10.1016/j.csda.2012.12.001)
.

## References

Arthur, D., and S. Vassilvitskii. 2007. “K-Means++ the Advantages of
Careful Seeding.” In Proceedings of the Eighteenth Annual ACM-SIAM
Symposium on Discrete Algorithms, 1027–35.

Marron, J. S., J. O. Ramsay, L. M. Sangalli, and A. Srivastava. 2014.
“Statistics of Time Warpings and Phase Variations.”

Marron, J. S., J. O. Ramsay, L. M. Sangalli, and A. Srivastava. 2015.
“Functional Data Analysis of Amplitude and Phase Variation.” Statistical
Science, 468–84.

Ramsay, J., and B. W. Silverman. 2005. Functional Data Analysis.
Springer Series in Statistics. Springer.

Sangalli, L. M., P. Secchi, S. Vantini, and V. Vitelli. 2010. “K-Mean
Alignment for Curve Clustering.” Computational Statistics & Data
Analysis 54 (5): 1219–33.

Tucker, J. D., W. Wu, and A. Srivastava. 2013. “Generative Models for
Functional Data Using Phase and Amplitude Separation.” Computational
Statistics & Data Analysis 61: 50–66.

Vantini, S. 2012. “On the Definition of Phase and Amplitude Variability
in Functional Data Analysis.” Test 21 (4): 676–96.
https://doi.org/10.1007/s11749-011-0268-9.

## See also

Useful links:

- <https://astamm.github.io/fdacluster/>

- <https://github.com/astamm/fdacluster>

## Author

**Maintainer**: Aymeric Stamm <aymeric.stamm@cnrs.fr>
([ORCID](https://orcid.org/0000-0002-8725-3654))

Other contributors:

- Laura Sangalli <laura.sangalli@polimi.it> \[contributor\]

- Piercesare Secchi <piercesare.secchi@polimi.it> \[contributor\]

- Simone Vantini <simone.vantini@polimi.it> \[contributor\]

- Valeria Vitelli <valeria.vitelli@medisin.uio.no> \[contributor\]

- Alessandro Zito <zito.ales@gmail.com> \[contributor\]
