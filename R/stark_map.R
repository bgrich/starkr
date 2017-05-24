#' Construct Zero Field Energy Matrix.
#'
#' \code{zero_field_energy_mat} constructs a matrix with the zero field
#' energies on the diagonal.
#'
#' This function takes vectors for n, l, and j and returns a matrix with the
#' energies for those n, l, and j states on the diagonal. To create the
#' vectors for n, l, and j, \code{\link{state_list}} should be used. This will
#' create a matrix with all the combinations of n, l, and j needed to create
#' the Stark matrix. The energy is determined by:
#' \deqn{E = -1 / (2* (n - \delta(n, l, j))^2)}
#' where \eqn{\delta(n, l, j)} is determined by \code{\link{quantum_defect}}.
#'
#' @param n A numeric. A vector that contains a series of principle quantum
#' numbers.
#' @param l A numeric. A vector that contains a series of orbital angular
#' momentum quantum numbers.
#' @param j A numeric. A vector that contains a series of total angular
#' momentum quantum numbers.
#'
#' @export
zero_field_energy_mat <- function(n, l, j){
  # Initializes a matrix to hold the zero field energy
  Energy <- matrix(0, nrow = length(n), ncol = length(n))

  # For loop puts zero field energies on the diagonal of the matrix
  for(i in 1:length(n)){
    Energy[i, i] <- -1 / (n[i] - quantum_defect(n[i], l[i], j[i])) ^ 2 / 2
  }

  Energy

}

#' Construct Zero Field Energy Data Frame.
#'
#' \code{zero_field_energy_df} constructs a data frame with the zero field
#' energies and information about the states.
#'
#' This function takes vectors for n, l, j, and a single number mj and returns a
#' matrix with the energies for those (n, l, j, mj) states on the diagonal. To
#' create the vectors for n, l, and j, \code{\link{state_list}} should be used.
#' This will create a matrix with all the combinations of n, l, and j needed to
#' create the Stark matrix. The energy is determined by: \deqn{E = -1 / (2* (n -
#' \delta(n, l, j))^2)} where \eqn{\delta(n, l, j)} is determined by
#' \code{\link{quantum_defect}}.
#'
#' @param n A numeric. A vector that contains a series of principle quantum
#'   numbers.
#' @param l A numeric. A vector that contains a series of orbital angular
#'   momentum quantum numbers.
#' @param j A numeric. A vector that contains a series of total angular momentum
#'   quantum numbers.
#' @param mj A numeric. A single number that represents the magnetic momentum
#'   quantum number.
#'
#'   @export
zero_field_energy_df <- function(n, l, j, mj){

  # Initializes a data frame to hold the zero field energy and information
  # about the states
  Energy_df <- data.frame(E0 = numeric(), n = numeric(), l = numeric(), j = numeric(), mj = numeric(), state = character())

  # For loop puts all of the zero energies and states into the data frame
  for(i in 1:length(n)){
    new_row <- data.frame(E0 = -1 / (n[i] - quantum_defect(n[i], l[i], j[i])) ^ 2 / 2, n = n[i], l = l[i], j = j[i], mj = mj, state = paste(n[i],l[i], j[i], mj, sep = ','))

    Energy_df <- rbind(Energy_df, new_row)
  }
  #Turns the data frame in to a dplyr table.
  Energy_df <- dplyr::tbl_df(Energy_df)

  Energy_df
}

#' Builds a Stark Matrix.
#'
#' \code{stark_matrix} builds a Stark matrix.
#'
#' This function builds a Stark matrix to be used in the creation of a Stark
#' map. The function accepts vectors n, l, and j that describe all of the
#' quantum numbers to be used in building the Stark matrix. These vectors can
#' be obtained using \code{\link{state_list}} and are the columns of its
#' result. Only one value of mj should be provided as the Stark matricies
#' are mj independent. The Stark matrix is treated as symmetric and is only
#' computed for the uper right triangle of the matrix. Copies of those
#' elements are made in to their symmetric counterparts on the lower left
#' trinagle of the matrix.
#'
#' @param n A numeric. The vector of principle quantum numbers.
#' @param l A numeric. The vector of orbital angular momentum quantum numbers.
#' @param j A numeric. The vector of total angular momentum quantum numbers.
#' @param mj A numeric. The magnetic momentum quantum number.
#'
#' @export
stark_matrix <- function(n, l, j, mj){

  #Initializes and computes the Stark Matrix

  StarkMatrix <- matrix(, nrow = length(n), ncol = length(n))

  #Fills the Stark matrix. Treats the Stark matrix as symmetric and computes
  #only the elements for the upper right triangle of the matrix. Copies those
  #into the symmetric terms on the lower left triangle of the matrix.

  for(i in 1:length(n)){
    print(paste("Current row being processed: ", i, sep = ""))
    for(k in i:length(n)){
      StarkMatrix[i, k] <- stark_matrix_elem(n[i], n[k], l[i], l[k], j[i], j[k], mj, mj)
      StarkMatrix[k, i] <- StarkMatrix[i,k]
    }
  }

  StarkMatrix

}


#' Eigenvalues of the Stark Matrix.
#'
#' \code{stark_eigen} calculates the eigenvalues of the Stark matrix for a
#' sequence of electric field sizes.
#'
#' This function calculates the eigenvalues of the sum of the the zero electric
#' field matrix and the Stark matrix times the electric field. The eigenvalues
#' are output as a matrix with the different states as the columns and each
#' field step as the rows. The columns are in descending order and the rows
#' are in the initial state order of the zero electric field matrix and Stark
#' matrix.
#'
#' @param stark_matrix A matrix created using the function
#'   \code{\link{stark_matrix}}.
#' @param zero_field_mat A matrix created using the function
#'   \code{\link{zero_field_energy_mat}}.
#' @param field_min A numeric that sets the mininum of the electric field.
#' @param field_max A numeric that sets the maximum of the electric field.
#' @param step_size A numeric that sets the step size of the electric field.
#'
#' @export
stark_eigen <- function(stark_matrix, zero_field_mat, field_min, field_max, step_size){

  # Determines the number of electric field points and the step size.
  field <- seq(field_min, field_max, by = step_size)
  field.au <- field / 5.142e9
  Energy <- numeric()

  # Diagonalizes a matrix of the zero field energy + the Stark matrix times the
  # electric field. Also builds a matrix of the eigenvalues
  for(k in 1:length(field)){
    new_row <- eigen(zero_field_mat + stark_matrix * field.au[k])$values
    Energy <- rbind(Energy, new_row)
  }

  Energy

}


#' Tidy Stark Energy Data Frame.
#'
#' \code{tidy_stark_energy} creates a tidy data frame of the Stark energy
#'
#' This function creates a tidy data frame of the Stark energy. The zero energy
#' data frame is provided by the function \code{\link{zero_field_energy_df}}.
#' The Stark energy matrix is provided by the function
#' \code{\link{stark_eigen}}. The field_min, field_max, and step_size should
#' be set to same values as when \code{stark_eigen} was run.
#'
#' @param zero_frame A data frame that contains the zero electric field
#' energies and state information.
#' @param stark_energy A matrix that contains the eigen energies for the
#' states in the columns and the field levels in the rows.
#' @param field A vector containing the field positions. Must match that used
#' for \code{stark_eigen}.
#'
#' @export
tidy_stark_energy <- function(zero_frame, stark_energy, field){
  # field <- seq(field_min, field_max, by = step_size)

  #Turns the data frame in to a dplyr table.
  ZeroFrame <- dplyr::tbl_df(zero_frame)

  #Determines which direction the Energy eigenvalues are going and arranges the zero energy data frame to match. If the minimum zero field energy is the last column, then the Zero Field data frame is put in descending order. If the min zero field energy is the first column, then the Zero Field data frame is put in ascending order.
  if(which.min(stark_energy[1,])>1){
    ZeroFrame <- ZeroFrame %>%
      dplyr::arrange(desc(E0), desc(l))
  } else{
    ZeroFrame <- ZeroFrame %>%
      dplyr::arrange(E0, l)
  }

  EnergyDataFrame <- data.frame(Field = numeric(), E = numeric(), E0 = numeric(), n = numeric(), l = numeric(), j = numeric(), mj = numeric(), state = character())

  #Creates a tidy data frame of the Energy eigen states at all fields
  for(k in 1:length(stark_energy[1,])){
    new_chunk <- data.frame(Field = field, E = stark_energy[,k],E0 = ZeroFrame$E0[k], n = ZeroFrame$n[k], l = ZeroFrame$l[k], j = ZeroFrame$j[k], mj = ZeroFrame$mj[k], state = ZeroFrame$state[k])
    EnergyDataFrame <- rbind(EnergyDataFrame, new_chunk)
  }

  EnergyDataFrame <- dplyr::tbl_df(EnergyDataFrame)

  EnergyDataFrame

}
