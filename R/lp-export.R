#' Linear and integer programming
#'
#' Interface to `lp_solve` linear/integer programming system.
#'
#' This function calls the `lp_solve` 5.5 solver. That system has many options
#' not supported here. The current version is maintained at
#' https://lpsolve.sourceforge.net/5.5/.
#'
#' Note that every variable is assumed to be >= 0!
#'
#' @inheritParams lpSolve::lp
#'
#' @return An [lpSolve::lp.object] object.
#'
#' @keywords internal
#' @export
lp <- function(direction = "min",
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
               timeout = 0L) {
  lpSolve::lp(
    direction = direction,
    objective.in = objective.in,
    const.mat = const.mat,
    const.dir = const.dir,
    const.rhs = const.rhs,
    transpose.constraints = transpose.constraints,
    int.vec = int.vec,
    presolve = presolve,
    compute.sens = compute.sens,
    binary.vec = binary.vec,
    all.int = all.int,
    all.bin = all.bin,
    scale = scale,
    dense.const = dense.const,
    num.bin.solns = num.bin.solns,
    use.rw = use.rw,
    timeout = timeout
  )
}
