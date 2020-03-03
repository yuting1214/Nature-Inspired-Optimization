##########################################################################
#Main function: QPSO
#Initial point: Random LHQ
#Bound handling:Nearst method
#Beta = 1/0.96
##########################################################################

QPSO <- function(objective, bound = c(-100,100), particle_num, dimension, iteration,  verbose = FALSE){
  #Bound
  LowerBound <- bound[1] ; UpperBound <- bound[2] ; interval <- abs(UpperBound - LowerBound)
  #Tuning parameter
  Beta = 1/0.96
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
  
  #Objective function evaluation
  evaluationo_ini = apply(population, 1, objective)
  initial_global_index = which(evaluationo_ini == min(evaluationo_ini))
  global_point = population[initial_global_index,] #first global point
  #Iteration 
  for(i in seq_len(iteration)){
    #Update each particle
    for( j in seq_len(particle_num)){
      #Update current best
      current_point <- particle_path[[j]][1, ]
      current_best_point <- particle_path[[j]][2, ]
      check_pbest <- objective(current_point) < objective(current_best_point)
      if(check_pbest){
        particle_path[[j]][2, ] <- current_point
        current_best_point <-current_point
      }
      check_gbest <- objective(current_best_point) < objective(global_point)
      if(check_gbest){global_point <- current_best_point}
      #Find mbest
      current_best_list <- sapply(1:length(particle_path), function(i) particle_path[[i]][2,])
      mbest <- apply(current_best_list, MARGIN = 1, mean)
      #Update each particle's dimension
      for( k  in seq_len(dimension)){
        phi_1 <- runif(1);phi_2 <- runif(1)
        #local attractor
        p_la <- (phi_1 * current_best_point[k] + phi_2* global_point[k]) / (phi_1+phi_2)
        #Direction
        u <- runif(1)
        if(u >=0.5){
          temp_point = p_la - Beta * abs(mbest[k] - current_point[k]) * log(1/u)
        }else{
          temp_point = p_la + Beta * abs(mbest[k] - current_point[k]) * log(1/u)
        }
        #Check boundary
        if(temp_point>UpperBound){temp_point = UpperBound }
        if(temp_point<LowerBound){temp_point = LowerBound }
        particle_path[[j]][1,k] = temp_point
      }
    }
    if(verbose == T){cat("The", i, "th iteration is finished. \n")}
  }
  #Find global
  global = global_point
  global_optimal = objective(global)
  result = list(global_optimal_point = global, optimal_value = global_optimal)
  return(result)
}
testing_21 <- function(x){
  f_ID = c(4, 11, 5);s = c(10, 20, 30);l = c(1, 1e-6, 1);b = c(0, 100, 200)
  C_i(x, f_ID,s,l,b)
}
particle_num = 100;dimension = 2; iteration = 100; bound = c(-100,100)
QPSO(testing_21, bound , particle_num, dimension, iteration,verbose = T)

