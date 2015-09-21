# Requires as input the principal quantum number n, orbital angular momentum
# number l, and total angular momentum number j.
# For example, if l = 1 and j = 1/2, this is the np_1/2 state.

#' Quantum Defect Calculation.
#'
#' \code{quantum_defect} calculates the quantum defect for an arbitrary state
#' (n, l, j) for Rubidium-85
#'
#'  This function calculates the quantum defect for an arbitrary state
#'  (n, l, j) in Rubidium-85. Depending on the choices of l and j, the proper
#'  parameters will be selected and the quantum defect will be calculated
#'  according to equation (16.19) on pg. 351 of Rydberg Atoms by Gallagher.
#'
#'  @param n A numeric. The principle quantum number n.
#'  @param l A numeric. The orbital angular momentum quantum number l.
#'  @param j A numeric. The total angular momentum quantum number j.
#'
#'  @examples
#'  quantum_defect(32, 1, 1.5)
#'  quantum_defect(5, 0, 0.5)

quantum_defect <- function(n,l,j){

  #Chooses the quantum defect parameters based on the l and j quantum numbers.
  #The quantum defect parameters come from updated values from a variety of
  #papers.

  if(l == 0){
    delta_0 <- 3.1311804
    delta_2 <- 0.1784
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 1 & j == 1/2){
    delta_0 <- 2.6548849
    delta_2 <- 0.2900
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if (l == 1 & j == 3/2){
    delta_0 <- 2.6416737
    delta_2 <- 0.2950
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 2 & j == 3/2){
    delta_0 <- 1.348091
    delta_2 <- -0.60286
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 2 & j == 5/2){
    delta_0 <- 1.34646572
    delta_2 <- -0.59600
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 3 & j == 5/2){
    delta_0 <- 0.0165192
    delta_2 <- -0.085
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 3 & j == 7/2){
    delta_0 <- 0.0165437
    delta_2 <- -0.086
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  } else if(l == 4){
    delta_0 <- 0.00400
    delta_2 <- 0
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  }else{
    delta_0 <- 0
    delta_2 <- 0
    delta_4 <- 0
    delta_6 <- 0
    delta_8 <- 0
  }

  # Calculates the quantum defect delta_nlj based on Equation 16.19 from
  # Rydberg Atoms by Gallagher (Pg. 351)

  delta <- delta_0 + delta_2 / (n - delta_0) ^ 2 + delta_4 / (n - delta_0) ^ 4 + delta_6 / (n - delta_0) ^ 6 + delta_8 / (n - delta_0) ^ 8

  delta
}
