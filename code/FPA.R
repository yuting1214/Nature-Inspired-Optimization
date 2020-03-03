
FPA <- function(objective, bound = c(-100,100), flower_num, dimension, iteration, verbose = FALSE){
  #Bound
  LowerBound <- bound[1] ; UpperBound <- bound[2] ; interval <- abs(UpperBound - LowerBound)
  #Tuning parameter
  p <- 0.5
  lambda <- 1.5
  sigma_square <- ( (gamma(lambda+1)/(lambda * gamma( (1+lambda)/2 ))) * 
               ( sin(pi*lambda/2)/2^((lambda-1)/2) ))^(1/lambda)
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
  current_point <- LHD(flower_num, dimension, bound)$LHD
  #Objective function evaluation
  evaluationo_ini <- apply(current_point, 1, objective)
  initial_global_index <- which.min(evaluationo_ini)
  global_point <- current_point[initial_global_index, ] #first global point
  #Iteration 
  for(i in seq_len(iteration)){
    for ( j in seq_len(flower_num)){
      prob <- runif(1)
      #Global pollinate
      if(prob < p){
        #Levy distribution (Mantegna algorithm)
        U <- rnorm(1,0,sigma_square); V <- rnorm(1,0,1)
        s <- U/(abs(V))^(1/lambda)
        L_direction <- runif(dimension, LowerBound, UpperBound)
        next_point <- current_point[j,] + s * L_direction * (global_point - current_point[j,])
      }else{ #Local pollinate
        e <- runif(1)
        #flower constancy index
        fci <- sample(1:flower_num,2)
        next_point <- current_point[j,] + e * (current_point[fci[1],] - current_point[fci[2],])
      }
      #Check boundary
      #Nearset method
      next_point[next_point>UpperBound] <- UpperBound 
      next_point[next_point<LowerBound] <- LowerBound 
      #Update
      eval_current <- objective(current_point[j,])
      eval_next <- objective(next_point)
      if(eval_next < eval_current){
        current_point[j,] <- next_point
      }
    }
    #Find global
    global_eval <- apply(current_point, 1 ,objective)
    global_index <- which.min(global_eval)
    global_point <- current_point[global_index,]
    if(verbose == T){cat("The", i, "th iteration is finished. \n")}
  }
  global <- global_point
  global_optimal <- objective(global)
  result <- list(global_optimal_point = global, optimal_value = global_optimal)
  return(result)
}

#Testing function

f3 <- function(x){
  temp <- sum(x*x)
  temp2 <- sum(0.5*x)
  out <- c(temp, temp2^2, temp2^4)
  return(sum(out))
}
flower_num = 100;dimension = 10; iteration = 500; bound = c(-100,100)
FPA(Compo3, bound , flower_num, dimension, iteration, verbose = T) 
