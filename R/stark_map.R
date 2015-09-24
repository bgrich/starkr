zero_field_energy_mat <- function(n, l, j){
  # Initializes a matrix to hold the zero field energy
  Energy <- matrix(0, nrow = length(n), ncol = length(n))

  # For loop puts zero field energies on the diagonal of the matrix
  for(i in 1:length(n)){
    Energy[i, i] <- -1 / (n[i] - QuantumDefect(n[i], l[i], j[i])) ^ 2 / 2
  }

  Energy

}
