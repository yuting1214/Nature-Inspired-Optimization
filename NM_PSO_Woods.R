##########################################################################
#Main function: Nelder-Mead + PSO
# Nelder-Mead :
#Initial simplex: Pfeffer's method
#Bound handling: Nearst method
#Argument: Alpha = 1, GAmma = 2, rho = 0.5, sigma = 0.5

# PSO :
#Initial point: Random LHQ
#Bound handling: Nearst method
#Argument: Alpha = 2.05, Beta = 2.05, Theta = 0.72984

#Stopping criterion: Woods
##########################################################################

NM_PSO <-function(objective, bound = c(-100,100), dimension, iteration_UB, 
                  tolerance = 10^(-4), verbose = T){
  #Initial NM_population, 
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
    vertice[vertice> bound[2]] =  bound[2]
    vertice[vertice < bound[1]] =  bound[1]
    return(vertice)
  }
  Population_num_NM <-  dimension + 1
  Population_NM <- Pfeffers_method(dimension = dimension , bound = bound)
  #Initial PSO_population,
  #LHD
  LHD <- function(level, dimension, bound = bound){
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
  Population_num_PSO <-  (2 * dimension )
  Population_PSO <- LHD(Population_num_PSO, dimension, bound)$LHD
  #Initial PSO_velocity,
  velocity <- matrix(0, nrow = Population_num_PSO, ncol = dimension )
  #Total_population
  Population_total <- rbind(Population_NM, Population_PSO)
  #Sorted_population 
  evaluation_total <- apply(Population_total, 1, objective)
  Sorted_population <- Population_total[order(evaluation_total),]
  #Hyperparameter, NM
  Alpha <- 1
  GAmma <- 2
  rho <- 1/2
  sigma <- 1/2
  #Hyperparameter, PSO
  #Bound
  LowerBound <- bound[1] ; UpperBound <- bound[2]
  #Tuning parameter
  Alpha_PSO = 2.05
  Beta_PSO = 2.05
  theta = 0.72984
  #Iteraion
  trigger <- T
  iteration <- 1 
  while(trigger){
    #Nelder-mead {
    #Find N+1 best for  NM
    Vertices <- Sorted_population[1:Population_num_NM,]
    #Evaluation
    evaluation_NM <- apply(Vertices, 1, objective)
    ###Reflection
    #index
    vertices_high_index <- which.max(evaluation_NM)[1]
    vertices_sehigh_index <- which(evaluation_NM == max(evaluation_NM[-which.max(evaluation_NM)]))[1]
    vertices_low_index <-  which.min(evaluation_NM)[1]
    #vertice
    vertices_high <- Vertices[vertices_high_index,]
    vertices_cent <- apply(Vertices[-vertices_high_index, ], 2, mean)
    vertices_refl <- (1+Alpha) * vertices_cent - Alpha * vertices_high
    #value
    value_high <- objective(vertices_high)
    value_sehigh <- objective(Vertices[vertices_sehigh_index,])
    value_low <- objective(Vertices[vertices_low_index,])
    value_refl <- objective(vertices_refl)
    
    if(value_low <= value_refl & value_refl < value_sehigh){
      Vertices[vertices_high_index, ] <- vertices_refl
    }else if(value_refl < value_low){
      ###Expansion
      vertices_exp <- GAmma * vertices_refl + (1 - GAmma) * vertices_cent
      value_exp <- objective(vertices_exp)
      if(value_exp < value_refl){ 
        Vertices[vertices_high_index, ] <- vertices_exp
      }else{ Vertices[vertices_high_index, ] <- vertices_refl}
    }else{
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
    }
    #Check bound
    Vertices[Vertices>UpperBound] = UpperBound
    Vertices[Vertices<LowerBound] = LowerBound
    #Update
    Sorted_population[1:Population_num_NM,] <- Vertices
    # }
    ### PSO{
    #Sorted population before runing PSO
    evaluation_total <- apply(Sorted_population, 1, objective)
    Sorted_population <- Sorted_population[order(evaluation_total),]
    #Find 2N worst in whole population for  PSO
    population <- Sorted_population[-(1:Population_num_NM),]
    global_point_pso <- Sorted_population[1, ]
    #Find neighbor best
    neighbor_best_index <- c()
    evaluation_PSO <- apply(population, 1 , objective)
    evaluation_PSO_matrix <- rbind(1:length(evaluation_PSO), evaluation_PSO)
    for( j in seq_len(nrow(population)/2) ){
      index <- which.min(evaluation_PSO_matrix[2,(2*j-1):(2*j)])
      neighbor_best_index <- c(neighbor_best_index, evaluation_PSO_matrix[1,(2*j-1):(2*j)][index])
    } 
    #
    neighbor_best <- population[neighbor_best_index, ]
    neighbor_best_list <- lapply(seq_len(nrow(neighbor_best)), function(i) neighbor_best[i,])
    neighbor_best_double <- lapply(seq_len(length(neighbor_best_list)), function(i) rep(neighbor_best_list[[i]],2))
    neighbor_best_final <- matrix(unlist(neighbor_best_double), ncol = dimension, byrow = T)
    #Velocity
    e1 <- matrix(runif(dimension * nrow(population)),  nrow(population))
    e2 <- matrix(runif(dimension * nrow(population)),  nrow(population))
    velocity <- theta * velocity + Alpha_PSO*e1*(t(replicate(nrow(population),global_point_pso)) - population) +
      Beta_PSO *e2*(neighbor_best_final - population)
    
    #Check bound
    population <- population + velocity
    population[population > UpperBound ] <- UpperBound
    population[population < LowerBound ] <- LowerBound
    #Update
    Sorted_population[-(1:Population_num_NM),] <- population
    
    #Sort whole population after update
    evaluation_total <- apply(Sorted_population, 1, objective)
    Sorted_population <- Sorted_population[order(evaluation_total),]
    
    #Stopping criterion
    NM_population <- Sorted_population[1:Population_num_NM,]
    criterion_quant <-sqrt(apply((t(replicate(nrow(NM_population), NM_population[1,])) - NM_population)^2,1, sum))
    criterion <- max(criterion_quant)
    delta <- max(1, sqrt(sum(Sorted_population[1,]^2)) )
    Sc <- ( (1/delta) * criterion ) 
    if(Sc <= tolerance){trigger <- F}
    if(iteration == iteration_UB){trigger <- F}
    iteration <- iteration + 1
    if(verbose){cat("The", iteration-1, "th iteration is finished. \n ")}
  }
  global_point <- Sorted_population[1,]
  evaluation_final <- objective(global_point)
  return(list(optimal_point = global_point, optimal_value = evaluation_final, iteration = iteration-1))
}

f3 <- function(x){
  temp <- sum(x*x)
  temp2 <- sum(0.5*x)
  out <- c(temp, temp2^2, temp2^4)
  return(sum(out))
}
dimension = 10; iteration_UB = 10000
NM_PSO(Compo1, bound = c(-10,10), dimension, iteration_UB,tolerance = 10^(-4), verbose = T)


dimension = 2; iteration_UB = 10000
NM_PSO(Compo1, bound = c(-100,100), dimension, iteration_UB,tolerance = 10^(-4), verbose = T)


