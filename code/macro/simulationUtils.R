cumIncF <- function(Data, timeGrid, method = "Nelson-Aalen"){
  time = Data$time_
  delta = Data$status_
  nsamp <- length(time)
  ordinds <- order(time, -delta)
  
  time <- time[ordinds]
  delta <- delta[ordinds]
  
  ## Calculate the indices at which there is a jump (event)
  jumpinds <- c(1:nsamp)[delta==1]
  na <- rep(0,nsamp)
  oldsa <- rep(0,nsamp)
  
  # Rx = 1:
  risk <- nsamp
  jna <- ifelse(delta[1]==1,1.0/risk,0.0)
  na[1] <- jna
  
  for (i in 2:nsamp) {
    risk<-nsamp-i+1
    jna <- ifelse(delta[i]==1,1.0/risk,0.0)
    na[i] <- na[i-1]+jna
  }
  
  cumInc <- 1-exp(-na)
  
  cumIncStepF<- stepfun(time, c(0,cumInc))
  return(data.frame(timeGrid = timeGrid, cumInc = cumIncStepF(timeGrid)))
}


simulate_trial_exposure_mix <- function(n, 
                                        scenario = c("homogeneous","moderate","extreme"),
                                        arm = c("placebo","vrc01"),
                                        model = c("A","B","C","D","mixCD1", "mixCD2"),
                                        seed = NULL){                  
  
  if (!is.null(seed)) set.seed(seed)
  
  scenario <- match.arg(scenario)
  arm <- match.arg(arm)
  model <- match.arg(model)
 
  
  t_max <- 104
  t80 <- 80
  
  # ----------------------------
  # Lambda
  # ----------------------------
  lambda_i <- rep(NA,n)
  if (scenario=="homogeneous"){
    lambda_i <- rep(0.03,n)/365.5*7
  } else if (scenario=="moderate"){
    idx <- rbinom(n,1,0.10)
    lambda_i[idx==1] <- 0.12/365.5*7
    lambda_i[idx==0] <- 0.022/365.5*7
    #0.0252*exp(0.458*c)
    # x_norm <- rnorm(n, mean = 0.2, sd = 1)
    # lambda_i <- 0.0252*exp(0.458*x_norm)/365.5*7
  } else {
    idx <- rbinom(n,1,0.047)
    lambda_i[idx==1] <- 1/365.5*7
    lambda_i[idx==0] <- 0
  }
  
  if(model=="mixCD1"){
    mix_D_prob <- 0.25
  } else if(model=="mixCD2"){
    mix_D_prob <- 0.5
  }
  
  Y <- T <- rep(NA,n)
  C <- rep(0,n)
  D_tilde <- rep(NA,n)
  
  for (i in 1:n){
    
    lam <- lambda_i[i]
    
    if (lam==0){
      T[i] <- Inf
      Y[i] <- t_max
      C[i] <- 0
      next
    }
    
    # ----------------------------
    # exposures
    # ----------------------------
    t <- 0
    times <- c()
    strains <- c()
    
    while(TRUE){
      z <- rexp(1,lam)
      t <- t + z
      if (t > t_max) break
      times <- c(times,t)
      strains <- c(strains,rbinom(1,1,0.7))
    }
    
    if (length(times)==0){
      T[i] <- Inf
      Y[i] <- t_max
      C[i] <- 0
      next
    }
    
    # ----------------------------
    # PLACEBO
    # ----------------------------
    if (arm=="placebo"){
      T[i] <- times[1]
      Y[i] <- min(T[i],t_max)
      C[i] <- as.integer(T[i] < t_max)
      D_tilde[i] <- ifelse(C[i]==1,strains[1],NA)
      next
    }
    
    # ----------------------------
    # VRC01
    # ----------------------------
    delayed_flag <- FALSE
    T_found <- FALSE
    
    for (k in seq_along(times)){
      
      t_k <- times[k]
      d_k <- strains[k]
      
      # ==========================
      # AFTER WEEK 80
      # ==========================
      if (t_k >= t80){
        T[i] <- t_k
        D_tilde[i] <- d_k
        T_found <- TRUE
        break
      }
      
      # ==========================
      # RESISTANT BEFORE 80
      # ==========================
      if (d_k==1){
        
        if (model=="A"){
          T[i] <- t_k
        } else {
          T[i] <- t_k + rexp(1,1/8)
        }
        
        D_tilde[i] <- 1
        T_found <- TRUE
        break
      }
      
      # ==========================
      # SENSITIVE BEFORE 80
      # ==========================
      u <- runif(1)
      
      if (model %in% c("mixCD1", "mixCD2")){
        
        if (u < 0.25){
          # immediate detection
          T[i] <- t_k
          D_tilde[i] <- 0
          T_found <- TRUE
          break
        }
        
        # mixture of C (block) and D (delay)
        if (runif(1) < mix_D_prob){ # mix_D_prob = p_d/(p_d + p_b)
          # DELAY
          delayed_flag <- TRUE # all viruses that were not detected right away were delayed with mix_D_prob
        } else {
          # BLOCK
          next
        }
        
      } else {
        
        pb <- ifelse(model %in% c("A","C"),0.75,0)
        pd <- ifelse(model %in% c("B","D"),0.75,0) # delayed
        
        if (u < 0.25){
          T[i] <- t_k
          D_tilde[i] <- 0
          T_found <- TRUE
          break
        } else if (pb == 0.75){ # all viruses that were not detected right away were blocked
          next
        } else {
          delayed_flag <- TRUE # all viruses that were not detected right away were delayed
        }
      }
      
      # ==========================
      # Resolve delayed crossing 560
      # ==========================
      if (delayed_flag && k < length(times)){ #only if there is another exposure before week 104
        
        t_next <- times[k+1]
        
        if (t_next >= t80){
          
          detect_time <- t80 + runif(1,0,t_max - t80)
          
          if (detect_time < t_next){
            T[i] <- detect_time
            D_tilde[i] <- ifelse(model %in% c("B"), 0, 1) # delayed sensitive infection is detected as sensitive in Model B; delayed sensitive infection is detected as resistant in Model D or mixCD
          } else {
            T[i] <- t_next
            D_tilde[i] <- strains[k+1]
          }
          
          T_found <- TRUE
          break
        }else{
          # if the next exposure is before week 80, then the delayed crossing is not detected and the next exposure determines infection
          next
        }
      }
    }
    
    # ----------------------------
    # Final unresolved
    # ----------------------------
    if (!T_found){ 
      if (delayed_flag){# the first 1:k sensitive viruses are before Week 80 and the k+1 th infection is greater than week 104; the previous loop did not consider
        #print(c(k, times, strains, model))
        T[i] <- t80 + runif(1,0,t_max - t80)
        D_tilde[i] <- ifelse(model %in% c("B"), 0, 1)
      } else {
        T[i] <- Inf
      }
    }
    
    # ----------------------------
    # Censoring
    # ----------------------------
    Y[i] <- min(T[i],t_max)
    C[i] <- as.integer(T[i] < t_max)
    
    if (C[i]==0){
      D_tilde[i] <- NA
    }
  }
  
  return(data.frame(
    id=1:n,
    lambda=lambda_i,
    Y=Y,
    T=T,
    C=C,
    D_tilde=D_tilde
  ))
}