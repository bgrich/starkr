#Function that outputs a matrix with the n, l, and j states for a given mj, and range of n's. Gives the full manifold for nmin to nmax. Also takes range of nmin-7 to nmin-1 and nmax+1 to nmax+7 and gives all of the low angular momentum states (l < 6) and adds those in as well.
#Output takes the form n, l, j

#' List of States for Stark Matrix.
#'
#' \code{state_list} creates a matrix containing the states to be calculated in
#' the Stark matrix calculation.
#'
#' This function outputs a matrix with the n, l, and j states for a given m_j
#' and range of n's. It gives teh full manifold between nmin and max. It also
#' returns the low angular momentum states (l < 6) for a range of nmin to
#' n_add_min and nmax to n_add_max.
#'
#'
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
      nadd <- c(n_add_min:(nmin-1),(nmax+1):n_add_max)
    }

    for(i in 1:length(nadd)){
      l0 <- mj-1/2
      l <- l0
      for(j in (l0):(5)){
        if(l == l0){
          NumberMatrix <- rbind(NumberMatrix, c(nadd[i],l,mj))
        } else{
          NumberMatrix <- rbind(NumberMatrix, c(nadd[i],l,l-1/2))
          NumberMatrix <- rbind(NumberMatrix, c(nadd[i],l,l+1/2))
        }
        l <- l + 1
      }
    }
  }

  #Returns the number matrix
  NumberMatrix
}
