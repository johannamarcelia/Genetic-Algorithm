data_fitting <- function(x, crm){
#INPUTS
# x is time series data
# crm is a vector of changepoint info eg:(n, l_1, l_2,..., l_n)
#OUTPUTS
# logval is the logliklihood value of specified arima model using crm changepoints 

  n <- length(x)
  t <- 1:n
  y <- x/100
  cos.t1 <- cos(2*pi*t/12)
  sin.t1 <- sin(2*pi*t/12)
  cos.t2 <- cos(4*pi*t/12)
  sin.t2 <- sin(4*pi*t/12)
  no.seg <- crm[1] + 1
  mod.crm <- c(1, crm[-1], n+1)
  R <- matrix(0, nrow=n, ncol=no.seg) 
  for (j in 1:(length(mod.crm)-1)) {
    for (i in 1:n){
      if (i==mod.crm[length(mod.crm)] & j==(length(mod.crm)-1)){
        R[i,j] <- 1
      } else if (i>=mod.crm[j]&i<mod.crm[j+1]) {
        R[i,j] <- 1
      } else
        R[i,j] <- 0
    }
  }   
  S <- t * R/120  #divide by 12 or 120 to scale for arima function
  #model <- tryCatch(
  #  {arima(y, order=c(1,0,1), seasonal = c(1,0,0), include.mean=FALSE, 
  #	xreg = cbind(R, S, sin.t1, cos.t1, sin.t2, cos.t2), method='CSS') 
  #    }, error=function(msg){
  #    arima(y, order=c(1,0,1), seasonal = c(1,0,0), include.mean=FALSE, 
  #	xreg = cbind(R, S, sin.t1, cos.t1, sin.t2, cos.t2), method='CSS-ML') 
  #   }, error = function(msg){
  #      arima(y, order=c(1,0,1), seasonal = c(1,0,0), include.mean=FALSE, 
  #	xreg = cbind(R, S, sin.t1, cos.t1, sin.t2, cos.t2), method='ML') 
  #   })
        

  model <- arima(y, order=c(1,0,1), seasonal = c(1,0,0), include.mean=FALSE, 
	xreg = cbind(R, S, sin.t1, cos.t1, sin.t2, cos.t2), method='CSS') # the method was CSS-ML, but wouldn't run
  logval <- model$log
  fitvals <- y-model$resid
  dfv <- c(logval, fitvals)
  return(dfv)
}


penalty <- function(x, crm){
#INPUTS
# x is time series data
# crm is a vector of changepoint info eg:(n, l_1, l_2,..., l_n)
#OUTPUTS
# pnt is the penalty tern for the specified crm changepoints 

  n.obs <- length(x)
  m <- crm[1]
  if (m == 0) {
    pnt <- 0
  } else {
    tau <- c(crm[-1],n.obs+1)
    pnt <- log(m+1)+sum(log(diff(tau)))+sum(log(tau[-1]))
  }
  return(pnt)
}


fitness <- function(x, crm){
#Computes combination of datafitting and penalty function
  return(2 * (-data_fitting(x, crm)[1]) + penalty(x, crm))
}



##### GA CODE ######
changepoints <- function(x){
#INPUTS
#x = time series data
#OUTPUTS
#min_crm[p] = the number and changepoint location that minimizes the fitness function

  #Generate Initial CRM Population#
  no.obs <- length(x)
  no.crm <- 125 # This number can be changed
  sl <- 6   #segment length--min number of datapoints in a segment 
  crm_list <- c()
  no.cng <- (0)
  crm_list <- list()
  i <- 0
  mu <- mean(x) 
  sig <- sd(x)

  while(length(crm_list) < (no.crm-1)){
    i <- i+1
    if (i>1000){
      break
    } else{
      m <- sample(1:4, 1)                   # m is number of changepoints in ith crm
      pts <- sample(12:(no.obs-sl), m, replace = F) # randomly choose m points in data w/o replacement
      print(pts)
      if (length(pts) > 1){
      #print(abs(diff(sort(pts))))
      spacing = any(abs(diff(sort(pts)))<sl) #if pts are too close together they get thrown out and new pts are generated    
      if (spacing == TRUE){
        print("next")
        i = i-1
        next
      }
    }
    }
    check_dup <- length(c(crm_list, list(c(m, sort(pts))))) == length(unique(c(crm_list, list(c(m, sort(pts))))))
    print(check_dup)
    if (check_dup == TRUE){  
      crm_list[i] <- list(c(m, sort(pts)))
    } else {
      i <- i-1
      next
    }
  }
  crm_list <- (append(crm_list, no.cng))
  print(crm_list)
  #Initial CRM Population has been generated 

  
  #Set up iterations, matrices, and lists to store subsequent answers
  iterations <- 100 #this is arbitrary, can change
  solutions <- as.list(numeric(no.crm * iterations)) # nrow = iterations, ncol = no.crm)
  dim(solutions) <- c(iterations,no.crm)
  min_fit <- list()
  min_crm <- list()

  ### MAIN LOOP ###
  for (p in 1:iterations){
  print("iteration number")
  print(p)
    #Calulate Fitness Scores & Rank
    fit_list <- c()
    for (i in 1:length(crm_list)){
      #init <- c(rep(mu, length(unlist(crm_list[i]))), sig)
      #fit.optm_norm <- optim(par=init, fitness, x=x, crm=unlist(crm_list[i]), 
		#method="BFGS", hessian=TRUE)
 	print(unlist(crm_list[i]))
      fit_list[i]<- fitness(x, unlist(crm_list[i])) 
    }
    min_fit[p] <- min(fit_list)
    min_crm[p] <- crm_list[match(min(fit_list), fit_list)]

    rank_fit <- rank(-fit_list)
    rank <- c()                               #rank is the list of the probabilities for each crm
    for (i in 1:length(fit_list)){
      rank[i] <- rank_fit[i]/sum(rank_fit)
    }
                                              #store fitness scores and ranks in matrix
    matrix_data <- rbind(crm_list, fit_list)
    entry_list <- lapply(seq_len(ncol(matrix_data)), function(i) matrix_data[,i])
    for (i in 1:no.crm){
      solutions[p,i] <- entry_list[i]
    }

    #### INTERIOR LOOP ####
    # Create new generation of crms to run through main loop
    new_gen <- list()
    k <- 0 
    while (length(new_gen) < no.crm){
      if (k==0) {
        k = k + 1
        new_gen[k] <- list(unlist(crm_list[match(min(fit_list), fit_list)]))
      } else {
      ### Chosing the parents ###
      parents <- sample(crm_list, 2, replace=FALSE, prob=rank)
      p1 <- parents[1]
      p2 <- parents[2]
      p12_cng <- sort(c(unlist(p1)[-1], unlist(p2)[-1]))

      ### Crossover ###
      c_cng <- c()
      for (i in 1:length(p12_cng)){
        uni <- runif(1)
        if (uni < 0.5){
          c_cng <- sort(c(p12_cng[i], c_cng)) 
        }
      } 
      if (length(c_cng) >= 1){                  #check that changepoints are far enough apart,if not skip iteration
        spacing <- any(diff(c(1, c_cng, no.obs))<sl)
        if (spacing == TRUE) {
          next
        }
      }

      if (length(c_cng) == 0) {                #avoids changepoint of NULL size
        child <- c(0)
        } else {
        child <- c(length(c_cng), c_cng)
      }

      ### Fine Tune Shifting ###
      if (child[1] != 0){
      child_times <- child[-1]
  
      for (i in 1:length(child_times)){
        uni <- runif(1)
        if (uni < 0.3){                           #moves up with 30% prob
          child_times[i] <- child_times[i] + 1
        } else if (uni <0.6){
          child_times[i] <- child_times[i] - 1    #moves down with 30% prob
        }
      }
      child <- c(length(child_times), child_times)
      if (length(child_times) >= 1){                  #check that changepoints are far enough apart,if not skip iteration 
        spacing2 <- any(diff(c(1, c_cng, no.obs))<sl)
          if (spacing2 == TRUE) {
          print("NEXT IT IS WORKING WOHOOOO")
           next
          } 
      }


      ### Mutation ### 
      obs <- seq(sl,(no.obs-sl))
      uni <- runif(1)
      if (uni < 0.1){
        child_range <- c()
	  for (i in 1:child[1]){
	    child_range <- c(child_range,(-sl:sl) + child[i+1])
	  }	
        avail_times <- obs[(obs) %in% child_range == FALSE]
        mut_pt <- sample(avail_times, 1)
        times <- sort(c(unlist(child[-1]), mut_pt))
        child <- c(length(times), times)
      }
      }
       
      check <- length(c(new_gen, list(child))) == length(unique(c(new_gen, list(child))))
      #print(check)
	if (check == TRUE){
        k <- k + 1
        new_gen[k] <- list(child)
      } else {
	  next
      }      
    }
    }
    
    crm_list <- new_gen
  }
  print("min crm")
  print(min_crm[p])  
  ### SOLUTION ### 
  n <- length(x)
  t <- 1:n
  y <- x/100
  cos.t1 <- cos(2*pi*t/12)
  sin.t1 <- sin(2*pi*t/12)
  cos.t2 <- cos(4*pi*t/12)
  sin.t2 <- sin(4*pi*t/12)
  no.seg <- unlist(min_crm[p])[1] + 1
  mod.crm <- c(1, unlist(min_crm[p])[-1], n+1)
  R <- matrix(0, nrow=n, ncol=no.seg) 
  for (j in 1:(length(mod.crm)-1)) {
    for (i in 1:n){
      if (i==mod.crm[length(mod.crm)] & j==(length(mod.crm)-1)){
        R[i,j] <- 1
      } else if (i>=mod.crm[j]&i<mod.crm[j+1]) {
        R[i,j] <- 1
      } else
        R[i,j] <- 0
    }
  }   
  S <- t * R/120
  model <- arima(y, order=c(1,0,1), seasonal = c(1,0,0), include.mean=FALSE, 
	xreg = cbind(R, S, sin.t1, cos.t1, sin.t2, cos.t2), method='CSS')
  model$coef
  B_0 <- model$coef[4:(4+no.seg-1)]
  B_1 <- model$coef[(4+no.seg):(3+2*no.seg)]
  params <- c(list(B_0), list(B_1), (min_crm[p]))
  return(params)
}