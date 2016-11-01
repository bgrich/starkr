#' List of States for Stark Matrix.
#'
#' \code{state_list} creates a matrix containing the states to be calculated in
#' the Stark matrix calculation.
#'
#' This function outputs a matrix with the n, l, and j states for a given m_j
#' and range of n's. It gives teh full manifold between nmin and max. It also
#' returns the low angular momentum states (l < 6) for a range of nmin to
#' n_add_min and nmax to n_add_max. If n_add_min or n_add_max are <= 0, then
#' those sections will be ignored and no extra states will be added.
#'
#' The output takes the form of a matrix with columns n, l, and j.
#'
#' @param nmin A numeric. The minimum principle quantum number for the whole
#' manifold
#' @param nmax A numeric. THe maximum principle quantum number for the whole
#' manifold
#' @param mj A numeric. The magnetic momentum quantum number.
#' @param n_add_min A numeric. The minimum principle quantum number for low
#' angular momentum states.
#' @param n_add_max A numeric. The maximum principle quantum number for low
#' angular momentum states.
#'
#' @export
state_list <- function(nmin, nmax, mj, n_add_min = 7, n_add_max = 7){

  #Populates a vector of n's in the range nmin to nmax.
  n <- c(nmin:nmax)

  #Initializes and fills a matrix with all of the n, l, and j states for the Stark matrix
  StateList <- numeric()

  #The lowest l level is determined by l = mj - 1/2. For that l level, there is only a single j state j = mj. For all higher l's, there are two j states, j = l +/- 1/2.
  for(i in 1:length(n)){
    l0 <- mj - 1/2
    l <- l0
    for(j in (l0):(n[i]-1)){
      if(l == l0){
        StateList <- rbind(StateList, c(n[i],l,mj))
      } else{
        StateList <- rbind(StateList, c(n[i],l,l-1/2))
        StateList <- rbind(StateList, c(n[i],l,l+1/2))
      }
      l <- l + 1
    }
  }

  #Adds in the lower l states for some additional nearby n's.
  if((n_add_min <= 0) & (n_add_max <= 0)){

  } else{
    if(n_add_min <= 0){
      nadd <- c((nmax+1):n_add_max)
    } else if(n_add_max <= 0){
      nadd <- c(n_add_min:(nmin-1))
    } else {
      nadd <- c((nmin - n_add_min):(nmin-1),(nmax+1):(nmax + n_add_max))
    }

    for(i in 1:length(nadd)){
      l0 <- mj-1/2
      l <- l0
      for(j in (l0):(5)){
        if(l == l0){
          StateList <- rbind(StateList, c(nadd[i],l,mj))
        } else{
          StateList <- rbind(StateList, c(nadd[i],l,l-1/2))
          StateList <- rbind(StateList, c(nadd[i],l,l+1/2))
        }
        l <- l + 1
      }
    }
  }

  #Returns the number matrix
  StateList
}
