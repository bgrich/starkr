#' Stark Matrix Element.
#'
#' \code{stark_matrix_element} calculates the elements of a Stark matrix
#'
#' This function calculates contribution to the energy due to the Stark effect
#' for an arbitrary set of states (n1, l1, j1, mj1) and (n2, l2, j2, mj2). The
#' element is determined based on equation 10 from Zimmerman et al., PRA, 20,
#' 2251 (1979). For a non-zero matrix element, the conditions mj1 = mj2 and l2 =
#' l1 +/- 1 must be met.
#'
#' @param n1 A numeric. The principle quantum number of state 1.
#' @param n2 A numeric. The principle quantum number of state 2.
#' @param l1 A numeric. The orbital angular momentum number of state 1.
#' @param l2 A numeric. The orbital angular momentum number of state 2.
#' @param j1 A numeric. The total angular momentum number of state 1.
#' @param j2 A numeric. The total angular momentum number of state 2.
#' @param mj1 A numeric. The magnetic angular momentum number of state 1.
#' @param mj2 A numeric. The magnetic angular momentum number of state 2.
stark_matrix_elem <- function(n1, n2, l1, l2, j1, j2, mj1, mj2){

  # Determines if the two mj terms are the same, if they are not then the matrix
  # element is set to zero
  if(mj1 != mj2){
    StarkElem <- 0
  } else{
    # Determines if l2 = l1 +/-1. If not, the matrix element is set to zero.
    if((l2 == (l1 + 1)) | (l2 == (l1 - 1))){

      # Calculates the spherical harmonic matrix element for mj1+1/2. If it is
      # equal to zero, sets that term in the summation equal to zero.
      if(sphere_mat_element(l1, l2, mj1 + 1 / 2) == 0){

        SumPlus <- 0

      } else{

        SumPlus <- clebsch_gordan(l1, 1/2, mj1 + 1 / 2, -1 / 2, j1, mj1) * clebsch_gordan(l2, 1 / 2, mj1 + 1 / 2, -1 / 2, j2, mj1) * sphere_mat_element(l1, l2, mj1 + 1 / 2)

      }
      # Calculates the spherical harmonic matrix element for mj1-1/2. If it is
      # equal to zero, sets that term in the summation equal to zero.
      if(sphere_mat_element(l1, l2, mj1 - 1 / 2) == 0){

        SumMinus <- 0

      } else{

        SumMinus <- clebsch_gordan(l1, 1 / 2, mj1 - 1 / 2, 1 / 2, j1, mj1) * clebsch_gordan(l2, 1 / 2, mj1 - 1 / 2, 1 / 2, j2, mj1) * sphere_mat_element(l1, l2, mj1 - 1 / 2)

      }

      # Calculates the stark matrix element.
      StarkElem <- radial_matrix_element(1, n1, n2, l1, l2, j1, j2) * (SumPlus + SumMinus)

    } else{
      StarkElem <- 0
    }
  }

  StarkElem

}
