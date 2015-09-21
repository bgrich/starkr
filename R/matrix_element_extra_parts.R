#' Clebsch-Gordan Coeffcient.
#'
#' \code{clebsch_gordan} calculates the Clebsch-Gordan coefficient
#'
#' This function calculates the Clebsch-Gordan coeffcient for an arbitrary
#' selection of quantum numbers j_1, j_2, m_1, m_2, j, and m_j. The condition
#' that must be met is (m_1 + m_2) = m_j. The analytic form used to calculate
#' the Clebsch-Gordan coefficient can be found in Edmonds or Lindgren and
#' Morrison.
#'
#' @param j1 A numeric. The total angular momentum for state 1.
#' @param j2 A numeric. The total angular momentum for state 2.
#' @param m1 A numeric. The magnetic quantum number for state 1.
#' @param m2 A numeric. The magnetic quantum number for state 2.
#' @param j A numeric. The combined total angular momentum.
#' @param mj A numeric. The combined magnetic quantum number.
#'
#' @examples
#' clebsch_gordan(1, 1 / 2, 0, 1 / 2, 3 / 2, 1 / 2)
clebsch_gordan <- function(j1, j2, m1, m2, j, mj){

  #If case for the conditions that m1+m2 = mj
  if((m1 + m2)!= mj){
    CG <- 0
  } else {
    #Calculates the terms that are not in the summation
    NotSumNumerator <- sqrt((2 * j + 1) * factorial(j1 + j2 - j) * factorial(j1 - m1) * factorial(j2 - m2) * factorial(j + mj) * factorial(j - mj))

    NotSumDenominator <- sqrt(factorial(j1 + j2 + j + 1) * factorial(j1 - j2 + j) * factorial(-j1 + j2 + j) * factorial(j1 + m1) * factorial(j2 + m2))

    #Settings for the summation. Summation only counts up to s = 4999. If j's and m's are such that one of the terms would be beyond N, then the code will not handle all possible cases.
    N <- 5000
    s <- 0
    sum <- 0
    sum_s <- 0
    for(i in 1:N){
      D1 <- s
      D2 <- j1 - m1 - s
      D3 <- j - mj - s
      D4 <- j2 - j + m1 + s

      #If any of the denominator terms (D1 to D4) are negative then the term returns zero
      if((D1 < 0) | (D2 < 0) | (D3 < 0) | (D4 < 0)){
        sum_s <- 0
      } else{

        #Calculates the terms in the summation at this step of s
        sum_s <- (-1) ^ (s + j1 - m1) * factorial(j1 + m1 + s) * factorial(j2 + j - m1 - s) / (factorial(D1) * factorial(D2) * factorial(D3) * factorial(D4))

      }
      #Adds current term to the sum and increments s
      sum <- sum + sum_s
      s <- s + 1
    }
    CG <- (NotSumNumerator / NotSumDenominator) * sum
  }

  CG
}

#' Wigner 3j Symbol.
#'
#' \code{wigner_3j} calculates the Wigner 3j symbol using the Clebsch-Gordan
#' coefficient calculated in \code{\link{clebsch_gordan}}.
#'
#' This function calculates the Wigner 3j symbol using the Clebsch-Gordan
#' coefficient calculated in \code{\link{clebsch_gordan}}. The Wigner 3j symbol
#' can be written as follows:
#' \tabular{ccc}{
#' j1 \tab j2 \tab j3\cr
#' m1 \tab m2 \tab m3
#' }
#'
#' @param j1 A numeric.
#' @param j2 A numeric.
#' @param j3 A numeric.
#' @param m1 A numeric.
#' @param m2 A numeric.
#' @param m3 A numeric.
wigner_3j <- function(j1,j2,j3,m1,m2,m3){
  output <- (-1)^(j1-j2-m3)*Clebsch_Gordan(j1,j2,m1,m2,j3,-m3)/sqrt(2*j3+1)
  output
}

#Calculates the Clebsch-Gordan coefficient using the analytic form found in Cornwell's Group Theory of Physics for arbitrary j1, j2, m1, m2, j, and mj.
#' Clebsch-Gordan Coeffcient Alternate Version.
#'
#' \code{clebsch_gordan} calculates the Clebsch-Gordan coefficient
#'
#' This function calculates the Clebsch-Gordan coeffcient for an arbitrary
#' selection of quantum numbers j_1, j_2, m_1, m_2, j, and m_j. The condition
#' that must be met is (m_1 + m_2) = m_j. The analytic form used to calculate
#' the Clebsch-Gordan coefficient can be found in Cornwell's Group Theory of
#' Physics.
#'
#' @param j1 A numeric. The total angular momentum for state 1.
#' @param j2 A numeric. The total angular momentum for state 2.
#' @param m1 A numeric. The magnetic quantum number for state 1.
#' @param m2 A numeric. The magnetic quantum number for state 2.
#' @param j A numeric. The combined total angular momentum.
#' @param mj A numeric. The combined magnetic quantum number.
#'
#' @examples
#' clebsch_gordan2(1, 1 / 2, 0, 1 / 2, 3 / 2, 1 / 2)
clebsch_gordan2 <- function(j1, j2, m1, m2, j, mj){

  #If case for the conditions that m1+m2 = mj
  if((m1 + m2) != mj){
    CG <- 0
  } else {
    #Calculates the terms that are not in the summation
    NotSumTerm1 <- sqrt((2 * j + 1) * factorial(j1 + j2 - j) * factorial(j1 - j2 + j) * factorial(-j1 + j2 + j) / factorial(j1 + j2 + j + 1))

    NotSumTerm2 <- sqrt(factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) * factorial(j2 - m2) * factorial(j + mj) * factorial(j - mj))

    #Settings for the summation. Summation only counts up to s = 4999. If j's and m's are such that one of the terms would be beyond N, then the code will not handle all possible cases.
    N <- 5000
    s <- 0
    sum <- 0
    sum_s <- 0
    for(i in 1:N){
      D1 <- s
      D2 <- j1 + j2 - j - s
      D3 <- j1 - m1 - s
      D4 <- j2 + m2 - s
      D5 <- j - j2 + m1 + s
      D6 <- j - j1 - m2 + s

      #If any of the denominator terms (D1 to D6) are negative then the term returns zero
      if((D1 < 0) | (D2 < 0) | (D3 < 0) | (D4 < 0) | (D5 < 0) | (D6 < 0)){
        sum_s <- 0
      } else{
        #Calculates the terms in the summation at this step of s
        sum_s <- (-1) ^ s * (factorial(D1) * factorial(D2) * factorial(D3) * factorial(D4) * factorial(D5) * factorial(D6)) ^ (-1)

      }
      #Adds current term to the sum and increments s
    sum <- sum + sum_s
    s <- s + 1
    }
    CG <- (NotSumTerm1 * NotSumTerm2) * sum
  }

  CG
}

#Calculates the Spherical Harmonic contribution to the matrix element: <l,m|cos(theta)|l',m>. This is based on Zimmerman et al, PRA 20, 2251 (1979).

#' Spherical Matrix Element.
#'
#' \code{sphere_mat_element} calculates the spherical matrix element
#'
#' This function calculates the spherical matrix element for a given l, l',
#' and m_l. It finds the matrix element <l, m|cos(theta)|l', m> based on a
#' paper by Zimmerman et al, PRA 20, 2251 (1979).
#'
#' @param l A numeric. The orbital angular momentum of the starting state.
#' @param lprime A numeric. The orbital angular momentum of the final state.
#' @param ml A numeric. The magnetic momentum of both the starting and final
#' state. The m_l of each state must be equal.
sphere_mat_element <- function(l, lprime, ml){
  if((lprime==(l+1))|(lprime==(l-1))){
    #If lprime is equal to l-1, then this if statement is processed
    if(lprime == (l-1)){

      MatElem <- sqrt((l^2 - ml^2)/((2*l+1)*(2*l-1)))

    }
    #If lprime is equal to l+1, then this if statement is processed
    if(lprime == (l+1)){

      MatElem <- sqrt(((l+1)^2 - ml^2)/((2*l+3)*(2*l+1)))

    }
  } else{

    MatElem <- "Error. Must meet condition l\' = l +/- 1."

  }

  MatElem

}
