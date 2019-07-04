##########################################################################
#Main function: SPSO
#Initial point: Random LHQ
#Bound handling: Nearst method
#Alpha:2.05, Beta = 2.05, Theta = 0.72984
##########################################################################

#Main function
SPSO <- function(objective, bound = c(-100,100), particle_num, dimension, iteration, verbose = FALSE){
  #Bound
  LowerBound <- bound[1] ; UpperBound <- bound[2] ; interval <- abs(UpperBound - LowerBound)
  #Tuning parameter
  Alpha = 2.05
  Beta = 2.05
  #initial points generation, LHQ
  LHD <- function(level, dimension, bound = c(-100, 100)){
    if(is.null(bound)){
      print("Stop!")
    }
    initial_LHD <- function(level, dimension, bound){
      range <- abs(bound[2]-bound[1])
      bin <- range/level
      LB <-seq(bound[1], bound[2], by =bin)
      LB <- LB[-length(LB)]
      UB <- (LB + bin)
      interval <- rbind(LB, UB)
      LHD <- matrix(NA, nrow = level, ncol = dimension)
      for(i in 1:ncol(interval)){
        LHD[i,] <- runif(dimension, min = interval[1,i], max = interval[2, i])
      }
      #shuffe
      for(i in 1:ncol(LHD)){
        LHD[, i] <-  sample(LHD[,i])
      }
      return(LHD)
    }
    Phi_p <- function(D, p = 50){
      dist_set <- as.vector(dist(t(D)))
      index_set <- rep(NA,length(dist_set))
      for(i in 1:length(dist_set)){
        index_set[i] <-sum(dist_set ==dist_set[i])
      }
      phi_p <- sum((index_set * dist_set^(-p)))^(1/p)
      return(phi_p)
    }
    ini_LHD <- initial_LHD(level, dimension, bound)
    phi_p <- Phi_p(D = ini_LHD)
    return(list(LHD = ini_LHD, Phi_p_value = phi_p))
  }
  population <- LHD(particle_num, dimension, bound)$LHD
  #[[]]:particle index []:Current or So-far best point index
  #Create two rows  to store (1)Current point (2) Current best 
  particle_path = lapply(seq_len(nrow(population)), function(i) matrix(rep(population[i,],2), nrow = 2,byrow =T))
  #Velocity
  velocity =  lapply(seq_len(nrow(population)), function(i) rep(0, dimension))
  #Objective function evaluation
  evaluationo_ini = apply(population, 1, objective)
  initial_global_index = which.min(evaluationo_ini)
  global_point = population[initial_global_index, ] #first global point
  #Iteration 
  for(i in seq_len(iteration)){
    #Update velocity,current point, current best point in each particle
    for( j in seq_len(particle_num)){
      #Random number
      e1 = runif(dimension)
      e2 = runif(dimension)
      #inertia function Theta(t)
      theta = 0.72984
      current = particle_path[[j]][1, ]
      current_best = particle_path[[j]][2, ]
      velocity[[j]] = theta * velocity[[j]] + Alpha*e1*(global_point - current) + Beta *e2*(current_best - current)
      #Next point
      next_point = current + velocity[[j]]
      
      #Check boundary
      #Nearset method
      next_point[next_point>UpperBound] <- UpperBound 
      next_point[next_point<LowerBound] <- LowerBound 
      #Evaluation
      eval_current_best = objective(current_best)
      eval_next = objective(next_point)
      #Decide current best
      if(eval_next <= eval_current_best){
        particle_path[[j]][2, ] = next_point
      }
      particle_path[[j]][1, ] = next_point
    }
    #Find global
    current_best_set = lapply(seq_len(length(particle_path)), function(i) particle_path[[i]][2,])
    evaluation = sapply(current_best_set, objective)
    global_index = which.min(evaluation)
    if(length(global_index)!=1){ global_index = sample(global_index,1) }
    global_point = as.vector(current_best_set[[global_index]])
    if(verbose == T){cat("The", i, "th iteration is finished. \n")}
  }
  global = global_point
  global_optimal = objective(global)
  result = list(global_optimal_point = global, optimal_value = global_optimal)
  return(result)
}
#Example
f3 <- function(x){
  temp <- sum(x*x)
  temp2 <- sum(0.5*x)
  out <- c(temp, temp2^2, temp2^4)
  return(sum(out))
}
particle_num = 100;dimension = 10; iteration = 1000; bound = c(-100,100)
SPSO(Compo1, bound , particle_num, dimension, iteration, verbose = T)
