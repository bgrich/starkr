#' Construct Zero Field Energy Matrix.
#'
#' \code{zero_field_energy_mat} constructs a matrix with the zero field
#' energies on the diagonal.
#'
#' This function takes vectors for n, l, and j and returns a matrix with the
#' energies for those n, l, and j states on the diagonal. To create the
#' vectors for n, l, and j,\code{\link{state_list}} should be used. This will
#' create a matrix with all the combinations of n, l, and j needed to create
#' the Stark matrix. The energy is determined by:
#' \deqn{E = -1 / (2* (n - \delta(n, l, j))^2)}
#' where \eqn{\delta(n, l, j)} is determined by \code{\link{quantum_defect}}
#'
#' @param n A numeric. A vector that contains a series of principle quantum
#' numbers.
#' @param l A numeric. A vector that contains a series of orbital angular
#' momentum quantum numbers.
#' @param j A numeric. A vector that contains a series of total angular
#' momentum quantum numbers.
zero_field_energy_mat <- function(n, l, j){
  # Initializes a matrix to hold the zero field energy
  Energy <- matrix(0, nrow = length(n), ncol = length(n))

  # For loop puts zero field energies on the diagonal of the matrix
  for(i in 1:length(n)){
    Energy[i, i] <- -1 / (n[i] - quantum_defect(n[i], l[i], j[i])) ^ 2 / 2
  }

  Energy

}

