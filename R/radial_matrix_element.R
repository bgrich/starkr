#' Radial Matrix Element.
#'
#' \code{radial_matrix_element} calculates the radial matrix element.
#'
#' This function calculates the radial matrix element for two arbitrary states
#' (n1, l1, j1) and (n2, l2, j2). A Numerov algorithm is used to compute the
#' radial matrix elements as done in Appendix A of Zimmerman et al, PRA, 20, 2251
#' (1979). The scaling used in this function is \eqn{\xi = \sqrt r},
#' \eqn{\Psi = r^(3/4) R(r)} as done by Bhatti, Cromer, and Cooke, PRA, 24, 161
#' (1981).
#'
#' @param k A numeric. The power of r to be calculated over. To get a dipole
#' matrix element, k must be equal to 1. Default k = 1.
#' @param n1 A numeric. The principle quantum number of state 1.
#' @param n2 A numeric. The principle quantum number of state 2.
#' @param l1 A numeric. The orbital angular momentum number of state 1.
#' @param l2 A numeric. The orbital angular momentum number of state 2.
#' @param j1 A numeric. The total angular momentum number of state 1.
#' @param j2 A numeric. The total angular momentum number of state 2.
radial_matrix_element <- function(n1, n2, l1, l2, j1, j2, k = 1){
  # Number of Electrons
  Z <- 1

  # Quantum numbers for states 1 and 2. Primary quantum number n, orbital
  # angular momentum l, total angular momentum j
  delta1 <- quantum_defect(n1, l1, j1)
  delta2 <- quantum_defect(n2, l2, j2)
  E1 <- -1 / (2 * (n1 - delta1) ^ 2)
  E2 <- -1 / (2 * (n2 - delta2) ^ 2)

  # Inner and outer turning points, core radius for both states
  r1_O <- 2 * n1 * (n1 + 15)
  r1_I <- (n1 ^ 2 - n1 * sqrt(n1 ^ 2 - l1 * (l1 + 1)))
  r2_O <- 2 * n2 * (n2 + 15)
  r2_I <- (n2 ^ 2 - n2 * sqrt(n2 ^ 2 - l2 * (l2 + 1)))

  # Core radius as a function of core polarizability r = (a_c)^(1/3)
  # Core polarizability is a_c = 9.0760 for Rubidium
  core.radius <- (9.0760) ^ (1/3)

  # Determine which outer turning point is the larger to set as starting point
  r_0 <- max(r1_O, r2_O)

  # Defining scaled x-axis ksi = sqrt(r), step size h, and starting
  # point ksi_0 = sqrt(r_0)
  ksi_0 <- sqrt(r_0)
  h <- 0.01
  ksi_1 <- ksi_0 - h
  ksi_2 <- ksi_1 - h

  # Initial wavefunction guesses
  Psi1_0 <- 10 ^ -15
  Psi1_1 <- 10 ^ -14

  Psi2_0 <- 10 ^ -15
  Psi2_1 <- 10 ^ -14

  # Defining terms to be used in Numerov Algorithm
  ksi_iminus1 <- ksi_0
  ksi_i <- ksi_1
  ksi_iplus1 <- ksi_2

  Psi1_iminus1 <- Psi1_0
  Psi1_i <- Psi1_1

  Psi2_iminus1 <- Psi2_0
  Psi2_i <- Psi2_1

  # Establishing Numerov integration data frame
  ksi <- numeric()
  Psi1 <- numeric()
  Psi2 <- numeric()
  N1_i <- numeric()
  N2_i <- numeric()
  Psi12 <- numeric()

  if (r1_O < r2_O) {
    ksi <- c(ksi, ksi_0, ksi_1)
    Psi1 <- c(Psi1, 0, 0)
    Psi2 <- c(Psi2, Psi2_0, Psi2_1)
    N1_i <- c(N1_i, 0, 0)
    N2_i <- c(N2_i, 2 * ksi_0 ^ 2 * Psi2_0 ^ 2, 2 * ksi_1 ^ 2 * Psi2_1 ^ 2 )
    Psi12 <- c(Psi12, 0, 0)
  } else {
    ksi <- c(ksi, ksi_0, ksi_1)
    Psi1 <- c(Psi1, Psi1_0, Psi1_1)
    Psi2 <- c(Psi2, 0, 0)
    N1_i <- c(N1_i, 2 * ksi_0 ^ 2 * Psi1_0 ^ 2, 2 * ksi_1 ^ 2 * Psi1_1 ^ 2)
    N2_i <- c(N2_i, 0, 0)
    Psi12 <- c(Psi12, 0, 0)
  }

  # Numerov Algorithm
  # Iterates algorithm until the condition of ksi_(i+1) < sqrt(r_I)
  # or ksi_(i+1) < sqrt(core.radius) is met
  repeat {
    # When ksi_i is larger than the smallest of the starting points, the
    # normalization for the larger outer turning point accumulates while the
    # other does not.
    if (ksi_i > sqrt(min(r1_O, r2_O))) {
      #First statement is case when r1_O < r2_O, second statement is r2_O < r1_O
      if (r1_O < r2_O){
        g2_iplus1  <- -8 * (ksi_iplus1 ^ 2 * E2 + Z - (l2 + 1/4) * (l2 + 3 / 4) / (2 * ksi_iplus1 ^ 2))
        g2_i       <- -8 * (ksi_i ^ 2 * E2 + Z - (l2 + 1 / 4) * (l2 + 3/4) / (2 * ksi_i ^ 2))
        g2_iminus1 <- -8 * (ksi_iminus1 ^ 2 * E2 + Z - (l2 + 1 / 4) * (l2 + 3 / 4) / (2 * ksi_iminus1 ^ 2))

        Psi2_iplus1 <- (Psi2_iminus1 * (g2_iminus1 - 12 / h^2) + Psi2_i * (10 * g2_i + 24 / h ^ 2)) / (12 / h ^ 2 - g2_iplus1)

        N1_iplus1 <- 0
        N2_iplus1 <- 2 * ksi_iplus1 ^ 2 * Psi2_iplus1 ^ 2 * h

        if (ksi_iplus1 < sqrt(max(r1_I, r2_I)) | ksi_iplus1 < sqrt(core.radius)){
          break
        } else {
          ksi <- c(ksi, ksi_iplus1)
          Psi1 <- c(Psi1, 0)
          Psi2 <- c(Psi2, Psi2_iplus1)
          N1_i <- c(N1_i, N1_iplus1)
          N2_i <- c(N2_i, N2_iplus1)
          Psi12 <- c(Psi12, 0)
        }

      } else {
        g1_iplus1  <- -8 * (ksi_iplus1 ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_iplus1 ^ 2))
        g1_i       <- -8 * (ksi_i ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_i ^ 2))
        g1_iminus1 <- -8 * (ksi_iminus1 ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_iminus1 ^ 2))

        Psi1_iplus1 <- (Psi1_iminus1 * (g1_iminus1 - 12 / h^2) + Psi1_i * (10 * g1_i + 24 / h ^ 2)) / (12 / h ^ 2 - g1_iplus1)

        N1_iplus1 <- 2 * ksi_iplus1 ^ 2 * Psi1_iplus1 ^ 2 * h
        N2_iplus1 <- 0

        if (ksi_iplus1 < sqrt(max(r1_I, r2_I)) | ksi_iplus1 < sqrt(core.radius)) {
          break
        } else {
          ksi <- c(ksi, ksi_iplus1)
          Psi1 <- c(Psi1, Psi1_iplus1)
          Psi2 <- c(Psi2, 0)
          N1_i <- c(N1_i, N1_iplus1)
          N2_i <- c(N2_i, N2_iplus1)
          Psi12 <- c(Psi12, 0)

        }
      }


      if (r1_O < r2_O) {
        Psi2_iminus1 <- Psi2_i
        Psi2_i <- Psi2_iplus1
      } else {
        Psi1_iminus1 <- Psi1_i
        Psi1_i <- Psi1_iplus1
      }

    } else {

      g1_iplus1  <- -8 * (ksi_iplus1 ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_iplus1 ^ 2))
      g1_i       <- -8 * (ksi_i ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_i ^ 2))
      g1_iminus1 <- -8 * (ksi_iminus1 ^ 2 * E1 + Z - (l1 + 1 / 4) * (l1 + 3 / 4) / (2 * ksi_iminus1 ^ 2))
      g2_iplus1  <- -8 * (ksi_iplus1 ^ 2 * E2 + Z - (l2 + 1 / 4) * (l2 + 3 / 4) / (2 * ksi_iplus1 ^ 2))
      g2_i       <- -8 * (ksi_i ^ 2 * E2 + Z - (l2 + 1 / 4) * (l2 + 3 / 4) / (2 * ksi_i ^ 2))
      g2_iminus1 <- -8 * (ksi_iminus1 ^ 2 * E2 + Z - (l2 + 1 / 4) * (l2 + 3 / 4) / (2 * ksi_iminus1 ^ 2))

      Psi1_iplus1 <- (Psi1_iminus1 * (g1_iminus1 - 12 / h ^ 2) + Psi1_i * (10 * g1_i + 24 / h ^ 2)) / (12 / h ^ 2 - g1_iplus1)
      Psi2_iplus1 <- (Psi2_iminus1 * (g2_iminus1 - 12 / h ^ 2) + Psi2_i * (10 * g2_i + 24 / h ^ 2)) / (12 / h ^ 2 - g2_iplus1)

      Psi12_iplus1 <- 2 * Psi1_iplus1 * Psi2_iplus1 * ksi_iplus1 ^ (2 + 2 * k) * h
      N1_iplus1 <- 2 * ksi_iplus1 ^ 2 * Psi1_iplus1 ^ 2 * h
      N2_iplus1 <- 2 * ksi_iplus1 ^ 2 * Psi2_iplus1 ^ 2 * h

      if (ksi_iplus1 < sqrt(max(r1_I, r2_I)) | ksi_iplus1 < sqrt(core.radius)) {
        break
      } else {
        ksi <- c(ksi, ksi_iplus1)
        Psi1 <- c(Psi1, Psi1_iplus1)
        Psi2 <- c(Psi2, Psi2_iplus1)
        N1_i <- c(N1_i, N1_iplus1)
        N2_i <- c(N2_i, N2_iplus1)
        Psi12 <- c(Psi12, Psi12_iplus1)
      }
      Psi1_iminus1 <- Psi1_i
      Psi1_i <- Psi1_iplus1
      Psi2_iminus1 <- Psi2_i
      Psi2_i <- Psi2_iplus1
    }

    ksi_iminus1 <- ksi_i
    ksi_i <- ksi_iplus1
    ksi_iplus1 <- ksi_iplus1 - h
  }

  RadialMatrixElement <- sum(Psi12) / (sqrt(sum(N1_i)) * sqrt(sum(N2_i)))
  RadialMatrixElement
}
