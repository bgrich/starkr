#' Calculate Radial Matrix Elements Involving Stark States.
#'
#' \code{stark_radial_mat_elem} calculates the radial matrix element between
#' a given initial state and a Stark state.
#'
#' This function calculates the radial matrix elements between a given initial
#' state and a Stark state at some arbitrary field F. The Stark state that is
#' given as input is the state that a particular line would correspond with at
#' zero field. The Stark matrix required for the function is produced by
#' \code{\link{stark_matrix}} and must correspond to the provided n_min and
#' n_max.
#'
#' @param stark_mat A matrix that is produced using \code{\link{stark_matrix}}
#' that must correspond to n_min and m_max.
#' @param initial_state A character representing the initial state in the form
#' "n,l,j,m_j".
#' @param stark_state A character representing the state at zero field that
#' corresponds to the Stark state and has the form "n,l,j,m_j".
#' @param field A numeric representing the field strength at which the Stark
#' state is evaluated.
#' @param n_min A numeric giving the minimum n used in determining the size
#' of the Stark matrix.
#' @param n_max A numeric giving the maximum n used in determining the size
#' of the Stark matrix.
#' @param ... Additional parameters to be sent to other functions.
#'
#' @export
stark_radial_mat_elem <- function(stark_mat, initial_state, stark_state, field, n_min, n_max, n_add_min, n_add_max){

  # Converts the field to atomic units
  field_au <- field/5.142e9

  # Splits the initial and Stark state strings into numerics
  initial_state_num <- as.numeric(unlist(strsplit(initial_state, split = ",")))
  final_state_num <- as.numeric(unlist(strsplit(stark_state, split = ",")))

  # Creates the number matrix and sets the size
  number_matrix <- state_list(n_min, n_max, final_state_num[4], n_add_min, n_add_max)
  size <- nrow(number_matrix)

  # Initializes the zero-field energy matrix
  zero_field_energy <- matrix(0,nrow = size, ncol = size)

  # Initializes the zero-field data frame
  zero_field_frame <- data.frame(E0 = numeric(),
                                 n = numeric(),
                                 l = numeric(),
                                 j = numeric(),
                                 mj = numeric(),
                                 state = character())

  # Builds the zero-field matrix and data frame
  for(i in 1:size){

    zero_field_energy[i, i] <- -1 / (number_matrix[i, 1]- quantum_defect(number_matrix[i, 1],  number_matrix[i, 2], number_matrix[i, 3])) ^ 2 / 2

    new_Row <- data.frame(E0 = -1 / (number_matrix[i, 1] - quantum_defect(number_matrix[i, 1],  number_matrix[i, 2], number_matrix[i, 3])) ^ 2 / 2,
                          n = number_matrix[i, 1],
                          l = number_matrix[i, 2],
                          j = number_matrix[i, 3],
                          mj = final_state_num[4],
                          state = paste(number_matrix[i, 1],
                                        number_matrix[i, 2],
                                        number_matrix[i, 3],
                                        final_state_num[4], sep = ','))

    zero_field_frame <- rbind(zero_field_frame, new_Row)
  }

  #Turns the data frame in to a dplyr table.
  unordered_data_frame <- tbl_df(zero_field_frame)
  ordered_data_frame <- tbl_df(zero_field_frame)

  # Computes the eigen-energies for the zero field states
  energy_order <- eigen(zero_field_energy)$values

  # Computes the eigenvectors for an arbitrary field
  field_vectors <- eigen(zero_field_energy + stark_mat*field_au)$vectors

  # Orders the zero field data frame to match the energy order from eigen
  if(which.min(energy_order)>1){
    ordered_data_frame <- ordered_data_frame%>%
      arrange(desc(E0), desc(l))
  } else{
    ordered_data_frame <- ordered_data_frame%>%
      arrange(E0, l)
  }

  # Mutates the ordered and unordered data frames to have an index
  ordered_data_frame <- ordered_data_frame %>%
    mutate(index = c(1:nrow(ordered_data_frame)))

  unordered_data_frame <- unordered_data_frame %>%
    mutate(index = c(1:nrow(unordered_data_frame)))

  # Determines the column of the Stark state eigen vector
  column <- as.numeric(ordered_data_frame %>%
                         filter(state %in% stark_state) %>%
                         select(index))

  # Sets the Stark state eigenvector
  stark_vector <- field_vectors[, column]

  # Splits the state character into numerics for later use
  unordered_data_frame <- unordered_data_frame %>%
    group_by(index)%>%
    mutate(charstates = strsplit(as.character(state), split = ","))%>%
    mutate(n1 = as.numeric(unlist(charstates)[1]),
           l1 = as.numeric(unlist(charstates)[2]),
           j1 = as.numeric(unlist(charstates)[3]),
           mj1 = as.numeric(unlist(charstates)[4]))%>%
    ungroup()

  # Matches states with the l +/- 1 states of the initial state
  state_index <- unordered_data_frame%>%
    filter(l1 %in% (initial_state_num[2]-1) | l1 %in% (initial_state_num[2] + 1))%>%
    select(index)

  # Determines the matrix element for each possible state and sums them
  matrix_elem <- numeric()
  for(i in 1:length(state_index$index)){
    rad_mat_elem <- radial_matrix_element(initial_state_num[1],
                                          unordered_data_frame$n1[state_index$index[i]],
                                          initial_state_num[2],
                                          unordered_data_frame$l1[state_index$index[i]],
                                          initial_state_num[3],
                                          final_state_num[3],
                                          1)

    matrix_elem <- c(matrix_elem, rad_mat_elem * stark_vector[state_index$index[i]])
  }
  sum(matrix_elem)

}
