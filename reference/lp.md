# Linear and integer programming

Interface to `lp_solve` linear/integer programming system.

## Usage

``` r
lp(
  direction = "min",
  objective.in,
  const.mat,
  const.dir,
  const.rhs,
  transpose.constraints = TRUE,
  int.vec,
  presolve = 0,
  compute.sens = 0,
  binary.vec,
  all.int = FALSE,
  all.bin = FALSE,
  scale = 196,
  dense.const,
  num.bin.solns = 1,
  use.rw = FALSE,
  timeout = 0L
)
```

## Arguments

- direction:

  Character string giving direction of optimization: "min" (default) or
  "max."

- objective.in:

  Numeric vector of coefficients of objective function

- const.mat:

  Matrix of numeric constraint coefficients, one row per constraint, one
  column per variable (unless transpose.constraints = FALSE; see below).

- const.dir:

  Vector of character strings giving the direction of the constraint:
  each value should be one of "\<," "\<=," "=," "==," "\>," or "\>=".
  (In each pair the two values are identical.)

- const.rhs:

  Vector of numeric values for the right-hand sides of the constraints.

- transpose.constraints:

  By default each constraint occupies a row of const.mat, and that
  matrix needs to be transposed before being passed to the optimizing
  code. For very large constraint matrices it may be wiser to construct
  the constraints in a matrix column-by-column. In that case set
  transpose.constraints to FALSE.

- int.vec:

  Numeric vector giving the indices of variables that are required to be
  integer. The length of this vector will therefore be the number of
  integer variables.

- presolve:

  Numeric: presolve? Default 0 (no); any non-zero value means "yes."
  Currently ignored.

- compute.sens:

  Numeric: compute sensitivity? Default 0 (no); any non-zero value means
  "yes."

- binary.vec:

  Numeric vector like int.vec giving the indices of variables that are
  required to be binary.

- all.int:

  Logical: should all variables be integer? Default: FALSE.

- all.bin:

  Logical: should all variables be binary? Default: FALSE.

- scale:

  Integer: value for lpSolve scaling. Details can be found in the
  lpSolve documentation. Set to 0 for no scaling. Default: 196

- dense.const:

  Three column dense constraint array. This is ignored if const.mat is
  supplied. Otherwise the columns are constraint number, column number,
  and value; there should be one row for each non-zero entry in the
  constraint matrix.

- num.bin.solns:

  Integer: if all.bin=TRUE, the user can request up to num.bin.solns
  optimal solutions to be returned.

- use.rw:

  Logical: if TRUE and num.bin.solns \> 1, write the lp out to a file
  and read it back in for each solution after the first. This is just to
  defeat a bug somewhere. Although the default is FALSE, we recommend
  you set this to TRUE if you need num.bin.solns \> 1, until the bug is
  found.

- timeout:

  Integer: timeout variable in seconds, defaults to 0L which means no
  limit is set.

## Value

An [lpSolve::lp.object](https://rdrr.io/pkg/lpSolve/man/lp.object.html)
object.

## Details

This function calls the `lp_solve` 5.5 solver. That system has many
options not supported here. The current version is maintained at
https://lpsolve.sourceforge.net/5.5/.

Note that every variable is assumed to be \>= 0!
