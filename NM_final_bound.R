##########################################################################
#Main function: Nelder-Mead 
#Initial simplex: Pfeffer's method
#Bound handling: Nearst method
#Argument: Alpha = 1, GAmma = 2, rho = 0.5, sigma = 0.5
#Stopping criterion: sd 
##########################################################################

Nelder_Mead <- function(objective, bound = NULL, dimension, tolerance = 10^(-6), verbose = F){
  if(!is.null(bound)){
    UpperBound <- bound[2] ; LowerBound <- bound[1] 
  }
  #Pfeffers_method
  Pfeffers_method <- function(dimension, bound = NULL, deltausual = 0.05, deltazero = 0.0075){
    simplex_num <- dimension+1
    vertice <- matrix(NA, nrow = simplex_num, ncol = dimension) 
    #initial guess
    if(is.null(bound)){
      x0 <- rnorm(dimension)
    }else{
      x0 <- runif(dimension, min = bound[1], max = bound[2])
    }
    vertice[, 1:dimension] <- matrix(rep(x0, simplex_num), nrow = simplex_num, byrow = T)
    for(j in 2:simplex_num){
      if(x0[j-1]!=0){
        vertice[j, j-1] = vertice[j, j-1]+ deltausual *  x0[j-1] 
      }else{
        vertice[j, j-1] = deltazero
      }
    }
    return(vertice)
  }
  Initial_Simplex <- Pfeffers_method(dimension, bound)
  Vertices <- Initial_Simplex
  ###Iteration
  #Hyperparameter
  Alpha <- 1
  GAmma <- 2
  rho <- 1/2
  sigma <- 1/2
  trigger <- T
  iteration <- 1 
  while(trigger){
    #Evaluation
    evaluation <- apply(Vertices, 1, objective)
    ###Reflection
    #index
    vertices_high_index <- which.max(evaluation)[1]
    vertices_sehigh_index <- which(evaluation == max(evaluation[-which.max(evaluation)]))[1]
    vertices_low_index <-  which.min(evaluation)[1]
    #vertice
    vertices_high <- Vertices[vertices_high_index,]
    vertices_cent <- apply(Vertices[-vertices_high_index, ], 2, mean)
    vertices_refl <- (1+Alpha) * vertices_cent - Alpha * vertices_high
    #value
    value_high <- objective(vertices_high)
    value_sehigh <-objective(Vertices[vertices_sehigh_index,])
    value_low <- objective(Vertices[vertices_low_index,])
    value_refl <- objective(vertices_refl)
    
    if(value_low <= value_refl & value_refl < value_sehigh){
      iteration <- iteration + 1
      if(verbose){cat("The", iteration-1, "th iteration is finished. \n")}
      Vertices[vertices_high_index, ] <- vertices_refl
    }else if(value_refl < value_low){
      iteration <- iteration + 1
      if(verbose){cat("The", iteration-1, "th iteration is finished. \n")}
      ###Expansion
      vertices_exp <- GAmma * vertices_refl + (1 - GAmma) * vertices_cent
      value_exp <- objective(vertices_exp)
      if(value_exp < value_refl){ 
        Vertices[vertices_high_index, ] <- vertices_exp
      }else{ Vertices[vertices_high_index, ] <- vertices_refl}
    }else{
      iteration <- iteration + 1
      if(verbose){cat("The", iteration-1, "th iteration is finished. \n")}
      ###Contraction
      if(value_high >= value_refl & value_refl > value_sehigh){
        replace <- vertices_refl
      }else{ replace <- vertices_high}
      Vertices[vertices_high_index, ] <- replace
      vertices_con <- rho * replace + (1-rho) * vertices_cent
      value_con <- objective(vertices_con)
      if(value_con <= objective(replace)){
        Vertices[vertices_high_index, ] <- vertices_con
      }else{
        ###Shrink
        vertices_best <- Vertices[vertices_low_index,] #vector
        #updating exclude the best
        Vertices[-vertices_low_index, ] <- sigma * Vertices[-vertices_low_index, ] + (1-sigma) * vertices_best
      }
      #Bound check
      Vertices[Vertices>UpperBound] = UpperBound
      Vertices[Vertices<LowerBound] = LowerBound
    }
    #Stopping criterion
    Sc <- sd(apply(Vertices, 1, objective))
    if(Sc <= tolerance){trigger <- F}
  }
  evaluation_final <- apply(Vertices, 1, objective)
  optimal_point <- Vertices[which.min(evaluation_final), ] 
  optimal_value <- objective(optimal_point)
  return(list(optimal_point = optimal_point, optimal_value = optimal_value, iteration = iteration))
}

#####
f3 <- function(x){
  temp <- sum(x*x)
  temp2 <- sum(0.5*x)
  out <- c(temp, temp2^2, temp2^4)
  return(sum(out))
}
tolerance = 10^(-6); verbose = T
Nelder_Mead(objective=f3,bound = c(-100,100),dimension=5, tolerance,verbose )
