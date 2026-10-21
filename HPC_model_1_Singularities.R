
args <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(args[1])
#set.seed(seed)
set.seed(1003)

library(tictoc)
tic()
library(MASS) 
library(lmerTest) 
library(truncnorm)


# Function to create the base dataset per batch
create_dataset_per_batch <- function(Therapist_Part,Therapist_Full,n_R_Part,n_R_Full,N_Batch) {
  
  # --- Remove NA therapists and NA rep counts for Part-time ---
  valid_part <- !(is.na(Therapist_Part) | is.na(n_R_Part))
  Therapist_Part     <- Therapist_Part[valid_part]
  n_R_Part           <- n_R_Part[valid_part]
  
  # --- Remove NA therapists and NA rep counts for Full-time ---
  valid_full <- !(is.na(Therapist_Full) | is.na(n_R_Full))
  Therapist_Full     <- Therapist_Full[valid_full]
  n_R_Full           <- n_R_Full[valid_full]
  
  
  # Create Therapist categories for Part-time
  Therapist_Part_seq <- unlist(mapply(function(t, r) rep(t, each = r), 
                                      Therapist_Part, n_I * n_R_Part))
  
  # Create Therapist categories for Full-time
  Therapist_Full_seq <- unlist(mapply(function(t, r) rep(t, each = r), 
                                      Therapist_Full, n_I * n_R_Full))
  # Combine both therapist categories sequentially
  Therapist <- c( Therapist_Part_seq, Therapist_Full_seq)
  
  # Now decide N_Batch = length(Therapist)
  N_Batch <- length(Therapist)
  
  # Patient IDs AFTER therapist vector is fully built
  Patient <- 1:N_Batch
  
  # Build dataset
  data <- data.frame(Patient = Patient,
                     Therapist = factor(Therapist))
  
  return(data)
}

########################################################################################################################################
############################################################################################################
#####################################################################

## Intervention generation for design with non random therapist allocations


# Function to generate Intervention
intervention_generation <- function(data, n_I) {
  Int_name = "Int"
  # copy input
  newdata <- data
  
  # total rows
  n <- nrow(newdata)
  
  # check divisibility
  if (n %% n_I == 0) {
    
    # assign treatments evenly
    newdata$Int_name <- rep(seq_len(n_I), each = n / n_I)
    # Convert Int to factor
    newdata$Int_name <- factor(newdata$Int_name)
    #Randomly Shuffle the interventions
    newdata$Int_name<- sample(newdata$Int_name)
  }
  
  # check divisibility
  if (n %% n_I == 1) {
    
    # assign treatments evenly
    newdata$Int_name <- rep(seq_len(n_I), each = (n+1) / n_I)
    # Convert Int to factor
    newdata$Int_name <- factor(newdata$Int_name[1:(length(newdata$Int_name)-1)])
    #Randomly Shuffle the interventions
    
    newdata$Int_name<- sample(newdata$Int_name)
    
  }
  
  newdata
}

#####################################
# Intervention_generation_randomised

# Generates all combinations of Int = 1,2,...n_I x Therapist = 1,2,...n_T 
# Then randomise IT combinations to patients based on their replications
#####################################
intervention_generation_randomised <- function(data,n_I,Therapist_Part_1,Therapist_Full_1,n_R_Part,n_R_Full,
                                               Int_name       = "Int",
                                               Therapist_name = "Therapist") {
  # Generates all combinations of Int = 1,2,...n_I x Therapist = 1,2,...n_T such that all combinations appear an equal number of times within the block.
  
  # Count how many rows we have
  n <- nrow(data)
  
  
  n_Int <- as.character((1:n_I))
  
  # --- Clean Part-time inputs ---
  valid_part <- !(is.na(Therapist_Part_1) | is.na(n_R_Part))
  Therapist_Part_1_clean <- Therapist_Part_1[valid_part]
  n_R_Part_clean         <- n_R_Part[valid_part]
  
  # --- Clean Full-time inputs ---
  valid_full <- !(is.na(Therapist_Full_1) | is.na(n_R_Full))
  Therapist_Full_1_clean <- Therapist_Full_1[valid_full]
  n_R_Full_clean         <- n_R_Full[valid_full]
  
  
  # --- Generate Part-time IT combinations ---
  IT_Part <- as.vector(outer(n_Int, Therapist_Part_1_clean, paste, sep = "."))
  IT_Part_rep <- unlist(
    mapply(function(t, r) rep(t, each = r),
           IT_Part, n_R_Part_clean)
  )
  
  # --- Generate Full-time IT combinations ---
  IT_Full <- as.vector(outer(n_Int, Therapist_Full_1_clean, paste, sep = "."))
  IT_Full_rep <- unlist(
    mapply(function(t, r) rep(t, each = r),
           IT_Full, n_R_Full_clean)
  )
  
  # --- Combine ---
  IT_All <- c(IT_Part_rep, IT_Full_rep)
  
  
  
  if (length(IT_All) != n) {
    stop(
      sprintf(
        "Error: IT combinations replicates are not equal to the total patients= %d,"
      )
    )
  }
  
  # Randomise IT combinations
  IT_randomised <- sample(IT_All)
  
  # Extract Interventions and Therapist part from randomised IT combinations)
  Int <- sub("\\..*", "",IT_randomised)
  Therapist <- sub(".*\\.", "", IT_randomised)
  
  # Add them to data
  data$Int <- factor(Int)  
  data$Therapist <- factor(Therapist) 
  
  
  data$IT <- factor(IT_randomised)
  
  return(data)
}

########################################################################################################################################
############################################################################################################
#####################################################################

Create_Data_Therapist_Selected_From_Available <- function(
    Therapist_Avail_Part, Therapist_Avail_Full,
    n_I, n_RP, n_RF, N_B, frac,local_indices_part,local_indices_full) {
  
  # Number of available therapists
  n_Avail_Part <- length(Therapist_Avail_Part)
  n_Avail_Full <- length(Therapist_Avail_Full)
  
  # Initialize selected vectors
  Therapist_Select_Part <- c()
  n_RP_Selected <- c()
  Therapist_Select_Full <- c()
  n_RF_Selected <- c()
  
  # Batch size with all available therapists
  N_B_1 <- n_I * (sum(n_RP) + sum(n_RF))
  
  # Case 1: batch smaller than target
  if (N_B_1 <= N_B) {
    Therapist_Select_Part <- Therapist_Avail_Part
    n_RP_Selected <- n_RP
    Therapist_Select_Full <- Therapist_Avail_Full
    n_RF_Selected <- n_RF
    N_Batch_Selected <- N_B_1
    
    print(paste("All the therapists are selected for this batch. The batch size is", 
                N_Batch_Selected, "which is less than the target", N_B))
    
    df <- create_dataset_per_batch(Therapist_Select_Part, Therapist_Select_Full,
                                   n_RP_Selected, n_RF_Selected, N_Batch_Selected)
  }
  
  # Case 2: batch slightly larger (within fraction)
  if ((N_B < N_B_1) & (N_B_1 <= ((1 + frac) * N_B))) {
    Therapist_Select_Part <- Therapist_Avail_Part
    n_RP_Selected <- n_RP
    Therapist_Select_Full <- Therapist_Avail_Full
    n_RF_Selected <- n_RF
    N_Batch_Selected <- N_B_1
    
    print(paste("All the therapists are selected for this batch. The batch size is", 
                N_Batch_Selected, "which is slightly more than the target", N_B))
    
    df <- create_dataset_per_batch(Therapist_Select_Part, Therapist_Select_Full,
                                   n_RP_Selected, n_RF_Selected, N_Batch_Selected)
  }
  
  # Case 3: batch too large, need to select subset
  if (N_B_1 > ((1 + frac) * N_B)) {
    
    # Combine replication counts
    x <- c(n_RP, n_RF)
    n <- length(x)
    
    # Target sum per intervention
    target <- ((1 + frac) * N_B) / n_I
    
    # Dynamic programming subset-sum solver
    subset_sum_dp <- function(x, target) {
      total <- sum(x)
      dp <- rep(FALSE, total + 1)
      dp[0 + 1] <- TRUE
      prev <- vector("list", total + 1)
      for (i in 0:total) prev[[i + 1]] <- integer(0)
      for (i in seq_along(x)) {
        for (s in total:0) {
          if (dp[s + 1]) {
            new_s <- s + x[i]
            if (!dp[new_s + 1]) {
              dp[new_s + 1] <- TRUE
              prev[[new_s + 1]] <- c(prev[[s + 1]], i)
            }
          }
        }
      }
      achievable <- which(dp)
      diffs <- abs(seq(0, total) - target)
      closest_sum <- achievable[which.min(diffs[achievable])]
      selected_indices <- prev[[closest_sum]]
      return(list(indices = selected_indices, sum = closest_sum))
    }
    
    res <- subset_sum_dp(x, target)
    random_sol <- rep(0, n)
    random_sol[res$indices] <- 1
    N_Batch_Selected <- n_I * sum(x[as.logical(random_sol)])
    
    # Assign selected therapists using logical indexing
    Therapist_Select_Part <- Therapist_Avail_Part[as.logical(random_sol[1:length(n_RP)])]
    n_RP_Selected <- n_RP[as.logical(random_sol[1:length(n_RP)])]
    
    Therapist_Select_Full <- Therapist_Avail_Full[as.logical(random_sol[(length(n_RP)+1):n])]
    n_RF_Selected <- n_RF[as.logical(random_sol[(length(n_RP)+1):n])]
    
    print(paste("Only a subset of the available therapists are selected for this batch. The batch size is", 
                N_Batch_Selected, "which is a lot more than the target", N_B))
    
    
    df <- create_dataset_per_batch(Therapist_Select_Part, Therapist_Select_Full,
                                   n_RP_Selected, n_RF_Selected, N_Batch_Selected)
  }
  
  # Return results
  return(list(
    df = df,
    Therapist_Select_Part = Therapist_Select_Part,
    Therapist_Select_Full = Therapist_Select_Full,
    n_RP_Selected = n_RP_Selected,
    n_RF_Selected = n_RF_Selected,
    N_Batch_Selected = N_Batch_Selected
  ))
}

########################################################################################################################################
############################################################################################################
#####################################################################



Randomise_Therapists_stage_1_Shapiro <- function(Therapist_Avail_Part, Therapist_Avail_Full, 
                                                 n_RP, n_RF, N_B1) {
  
  # Combine all available therapists
  Therapist_Avail_all <- as.numeric(c(Therapist_Avail_Part, Therapist_Avail_Full))
  
  # Combine capacities per therapist (can be unequal)
  Capacity_all_per_batch <- 2*c(n_RP, n_RF)
  names(Capacity_all_per_batch) <- Therapist_Avail_all
  
  # Initialize empty data frame
  df_T <- data.frame(Patient = numeric(0), Therapist = numeric(0))
  
  # Patients to assign
  Patient <- 1:N_B1
  
  for (i in Patient) {
    assigned <- FALSE
    
    # Only consider therapists with remaining capacity
    remaining_T <- Therapist_Avail_all[sapply(Therapist_Avail_all, function(x) {
      sum(df_T$Therapist == x) < Capacity_all_per_batch[as.character(x)]
    })]
    
    while (!assigned && length(remaining_T) > 0) {
      T <- sample(remaining_T, 1)
      cap_T <- Capacity_all_per_batch[as.character(T)]
      
      if (!is.na(cap_T) && sum(df_T$Therapist == T) < cap_T) {
        df_T <- rbind(df_T, data.frame(Patient = i, Therapist = T))
        assigned <- TRUE
      } else {
        remaining_T <- setdiff(remaining_T, T)
      }
    }
    
    if (!assigned) {
      message("No therapist available for patient ", i)
      break
    }
  }
  
  
  return(df_T)
}


########################################################################################################################################
############################################################################################################
#####################################################################

# Incomplete Block Design Selection and randomise the therapists in stage 1
Randomise_Therapists_stage_1_IBD <- function(Therapist_Avail_Part,Therapist_Avail_Full,
                                             n_I, n_RP, n_RF, n_T_Batch){
  
  
  
  Therapist_Avail_All<- c(Therapist_Avail_Part,Therapist_Avail_Full)
  
  if(n_T_Batch>length(Therapist_Avail_All)){
    
    print(paste("n_T_Batch:", n_T_Batch))
    print(paste("Available therapists:", length(Therapist_Avail_All)))
    
    stop("number of therapists available is less than the number of therapists to be chosen for Incomplete Block Design")
  }
  
  
  Therapist_Select_IBD <-Therapist_Avail_All[sample.int(length(Therapist_Avail_All), n_T_Batch)] #randomly select n_T_Batch therapists from all available ones
  
  # Initialize empty vectors
  Therapist_Select_Part <- c()
  Therapist_Select_Full <- c()
  n_RP_Selected <- c()
  n_RF_Selected <- c()
  
  # Assign therapists to part/full and get replication counts
  
  for (therapist in Therapist_Select_IBD) {
    
    if (length(Therapist_Avail_Part) > 0 &&
        therapist %in% Therapist_Avail_Part) {
      
      idx <- match(therapist, Therapist_Avail_Part)
      Therapist_Select_Part <- c(Therapist_Select_Part, therapist)
      n_RP_Selected <- c(n_RP_Selected, n_RP[idx])
      
    } else if (length(Therapist_Avail_Full) > 0 &&
               therapist %in% Therapist_Avail_Full) {
      
      idx <- match(therapist, Therapist_Avail_Full)
      Therapist_Select_Full <- c(Therapist_Select_Full, therapist)
      n_RF_Selected <- c(n_RF_Selected, n_RF[idx])
      
    } 
  }
  
  
  # Batch size with selected therapists
  sum_RP <- if (length(n_RP_Selected) > 0) sum(n_RP_Selected) else 0
  sum_RF <- if (length(n_RF_Selected) > 0) sum(n_RF_Selected) else 0
  
  N_Batch_IBD <- n_I * (sum_RP + sum_RF)
  
  
  
  print(paste(n_T_Batch, "Therapists are selected for this batch using Incomplete Block Design and",
              N_Batch_IBD, "patients are recruited." ))
  
  
  df_IBD<-create_dataset_per_batch(Therapist_Select_Part,Therapist_Select_Full,
                                   n_RP_Selected,n_RF_Selected,N_Batch_IBD)
  df_IBD<- df_IBD[sample(nrow(df_IBD)), ]
  rownames(df_IBD) <- NULL
  df_IBD$Patient<- 1:length(df_IBD$Patient)
  
  
  # Return results as a list
  return(list(
    df_IBD = df_IBD,
    Therapist_Select_Part = Therapist_Select_Part,
    Therapist_Select_Full = Therapist_Select_Full,
    n_RP_Selected = n_RP_Selected,
    n_RF_Selected = n_RF_Selected,
    N_Batch_IBD = N_Batch_IBD
  ))
  
}



# Function to generate treatment assignment for N patients using stage 2 intervention randomisation
Randomise_Interventions_Stage_2 <- function(ds,n_I) {
  
  # Add Intervention column
  ds$Int <- NA
  
  # Get unique therapists
  therapists <- unique(ds$Therapist)
  
  for (T in therapists) {
    # Patients for this therapist
    patients_idx <- which(ds$Therapist == T)
    n_patients <- length(patients_idx)
    
    # Generate balanced interventions (1,2,1,2,...)
    interventions <- rep(1:n_I, length.out = n_patients)
    
    # Randomize the order within this therapist
    interventions <- sample(interventions)
    
    # Assign to the data frame
    ds$Int[patients_idx] <- interventions
  }
  
  return(ds)
}

########################################################################################################################################
############################################################################################################
#####################################################################


samples_from_Design <- function(ds, del_0, del_1, del_2, del_3, del_4,
                                sig2_T, sig2_B, sig2_C,
                                sig2_IT, sig2_IB, sig2_IC, sig2_TB, sig2_CB, 
                                sig2_ITB, sig2_ICB, sig2_T_learn, sig2_C_learn, sig2_IT_learn, sig2_IC_learn,
                                sig2_E, time, matching) {
  
  ds <- tryCatch({
    
    # ----- Data preparation -----
    ds$Int       <- factor(ds$Int)
    ds$Therapist <- factor(ds$Therapist)
    ds$Batch     <- factor(ds$Batch)
    ds$Centre    <- factor(ds$Centre)
    
    ds$contrast <- ifelse(ds$Int == levels(ds$Int)[1], 1, -1)
    
    n_patients <- nrow(ds)
    centres    <- levels(ds$Centre)
    n_C        <- length(centres)
    
    # ----- time == 0 -----
    if(time == 0) {
      
      if (n_C == 1) {
        
        J_T   <- model.matrix(~ 0 + Therapist, ds) 
        J_B   <- model.matrix(~ 0 + Batch, ds)
        J_IT  <- model.matrix(~ 0 + Int:Therapist, ds)
        J_IB  <- model.matrix(~ 0 + Int:Batch, ds)
        J_TB  <- model.matrix(~ 0 + Therapist:Batch, ds)
        J_ITB <- model.matrix(~ 0 + Int:Therapist:Batch, ds)
        
        Cov_mat <- sig2_T   * (J_T  %*% t(J_T)) +
          sig2_B   * (J_B  %*% t(J_B)) +
          sig2_IT  * (J_IT %*% t(J_IT)) +
          sig2_IB  * (J_IB %*% t(J_IB)) +
          sig2_TB  * (J_TB %*% t(J_TB)) +
          sig2_ITB * (J_ITB %*% t(J_ITB)) +
          sig2_E   * diag(n_patients)
      }
      
      if (n_C > 1) {
        if (length(sig2_T) != n_C) stop("sig2_T must have exactly one variance per centre.")
        if (is.null(names(sig2_T)) || any(names(sig2_T) == "")) names(sig2_T) <- centres
        sig2_T <- sig2_T[centres]
        
        J_B   <- model.matrix(~ 0 + Batch, ds)
        J_C   <- model.matrix(~ 0 + Centre, ds)
        J_IT  <- model.matrix(~ 0 + Int:Therapist, ds)
        J_IB  <- model.matrix(~ 0 + Int:Batch, ds)
        J_IC  <- model.matrix(~ 0 + Int:Centre, ds)
        J_TB  <- model.matrix(~ 0 + Therapist:Batch, ds)
        J_CB  <- model.matrix(~ 0 + Centre:Batch, ds)
        J_ITB <- model.matrix(~ 0 + Int:Therapist:Batch, ds)
        J_ICB <- model.matrix(~ 0 + Int:Centre:Batch, ds)
        
        J_T_list <- list()
        for (c in centres) {
          idx_c <- ds$Centre == c
          Jc <- model.matrix(~ 0 + Therapist, droplevels(ds[idx_c, , drop = FALSE]))
          Jtmp <- matrix(0, n_patients, ncol(Jc))
          Jtmp[idx_c, ] <- Jc
          J_T_list[[c]] <- Jtmp
        }
        
        Cov_T <- matrix(0, n_patients, n_patients)
        for (c in centres) {
          Cov_T <- Cov_T + sig2_T[c] * (J_T_list[[c]] %*% t(J_T_list[[c]]))
        }
        
        Cov_mat <- Cov_T +
          sig2_B    * (J_B    %*% t(J_B)) +
          sig2_C    * (J_C    %*% t(J_C)) +
          sig2_IT   * (J_IT   %*% t(J_IT)) +
          sig2_IB   * (J_IB   %*% t(J_IB)) +
          sig2_IC   * (J_IC   %*% t(J_IC)) +
          sig2_TB   * (J_TB   %*% t(J_TB)) +
          sig2_CB   * (J_CB   %*% t(J_CB)) +
          sig2_ITB  * (J_ITB  %*% t(J_ITB)) +
          sig2_ICB  * (J_ICB  %*% t(J_ICB)) +
          sig2_E    * diag(n_patients)
      }
      
      J_fix <- cbind(intercept = 1, contrast = ds$contrast)
      mu <- J_fix %*% c(del_0, del_1)
      
      if(matching == 1){
        therapist_effects <- del_4 + rnorm(length(levels(ds$Therapist)))
        patSeverity <- therapist_effects[ds$Therapist] + rnorm(nrow(ds))
        J_fix <- cbind(J_fix, patSeverity)
        mu <- J_fix %*% c(del_0, del_1, 1)
      }
      
      ds$y <- as.vector(MASS::mvrnorm(n = 1, mu = mu, Sigma = Cov_mat))
    }
    
    # ----- time == 1 -----
    if(time == 1){
      
      ds <- ds[order(ds$Therapist, ds$Patient), ]
      ds$Patients_Before <- ave(ds$Patient, ds$Therapist, FUN = seq_along)
      ds$t_s <- as.numeric(scale(ds$Patients_Before, center = TRUE, scale = FALSE))
      
      if (n_C == 1) {
        J_T   <- model.matrix(~ 0 + Therapist, ds)
        J_IT  <- model.matrix(~ 0 + Int:Therapist, ds)
        J_T_learn  <- J_T * ds$t_s
        J_IT_learn <- J_IT * (ds$t_s * ds$contrast)
        
        Cov_mat <- sig2_T        * (J_T %*% t(J_T)) +
          sig2_IT       * (J_IT %*% t(J_IT)) +
          sig2_T_learn  * (J_T_learn %*% t(J_T_learn)) +
          sig2_IT_learn * (J_IT_learn %*% t(J_IT_learn)) +
          sig2_E        * diag(n_patients)
      }
      
      if (n_C > 1) {
        if (length(sig2_T) != n_C) stop("sig2_T must have exactly one variance per centre.")
        if (is.null(names(sig2_T)) || any(names(sig2_T) == "")) names(sig2_T) <- centres
        sig2_T <- sig2_T[centres]
        
        J_C   <- model.matrix(~ 0 + Centre, ds)
        J_IT  <- model.matrix(~ 0 + Int:Therapist, ds)
        J_IC  <- model.matrix(~ 0 + Int:Centre, ds)
        
        J_T_learn   <- model.matrix(~ 0 + Therapist, ds) * ds$t_s
        J_C_learn   <- J_C * ds$t_s
        J_IT_learn  <- J_IT * (ds$t_s * ds$contrast)
        J_IC_learn  <- J_IC * (ds$t_s * ds$contrast)
        
        J_T_list <- list()
        for (c in centres) {
          idx_c <- ds$Centre == c
          Jc <- model.matrix(~ 0 + Therapist, droplevels(ds[idx_c, , drop = FALSE]))
          Jtmp <- matrix(0, n_patients, ncol(Jc))
          Jtmp[idx_c, ] <- Jc
          J_T_list[[c]] <- Jtmp
        }
        
        Cov_T <- matrix(0, n_patients, n_patients)
        for (c in centres) {
          Cov_T <- Cov_T + sig2_T[c] * (J_T_list[[c]] %*% t(J_T_list[[c]]))
        }
        
        Cov_mat <- Cov_T +
          sig2_C        * (J_C        %*% t(J_C)) +
          sig2_IT       * (J_IT       %*% t(J_IT)) +
          sig2_IC       * (J_IC       %*% t(J_IC)) +
          sig2_T_learn  * (J_T_learn  %*% t(J_T_learn)) +
          sig2_C_learn  * (J_C_learn  %*% t(J_C_learn)) +
          sig2_IT_learn * (J_IT_learn %*% t(J_IT_learn)) +
          sig2_IC_learn * (J_IC_learn %*% t(J_IC_learn)) +
          sig2_E        * diag(n_patients)
      }
      
      J_fix <- cbind(
        intercept = 1,
        contrast = ds$contrast,
        t_s = ds$t_s,
        t_s_contrast = ds$t_s * ds$contrast
      )
      
      mu <- J_fix %*% c(del_0, del_1, del_2, del_3)
      
      if(matching == 1){
        therapist_effects <- del_4 + rnorm(length(levels(ds$Therapist)))
        patSeverity <- therapist_effects[ds$Therapist] + rnorm(nrow(ds))
        J_fix <- cbind(J_fix, patSeverity)
        mu <- J_fix %*% c(del_0, del_1, del_2, del_3, 1)
      }
      
      ds$y <- as.vector(MASS::mvrnorm(n = 1, mu = mu, Sigma = Cov_mat))
    }
    
    ds  # return dataset if no error
    
  }, error = function(e){
    message("Simulation failed in samples_from_Design: ", e$message)
    ds$y <- NA
    ds
  })
  
  return(ds)
}

#####################################################################################################################################
#####################################################################################################################################

generate_designs <- function(N,n_B,n_C,n_I,n_T,n_P,Therapist_Part,Therapist_Full,n_T_IBD,
                             n_RP,n_RF,frac,a,Add_Batch,Ex_C_P,Ex_C_F){
  
  N_B <- rep(0,n_C) # Initialize target batch sizes for flexible randmisation
  
  total_exp_cap <- sum(Ex_C_P*n_P+(Ex_C_F*(n_T - n_P)))  ### total expected capacity combining all the centres
  
  
  for (i in 1:n_C) {
    N_B[i]<- round(N*((Ex_C_P*n_P[i])+(Ex_C_F*(n_T[i] - n_P[i])))/(total_exp_cap*n_B))  # Target sample size in every centre 
  }
  
  
  ds_Non_rand<- data.frame()
  ds_orthogonal<- data.frame()      #initialise the final data frame for orthogonal design
  ds <- data.frame()      #initialise the final data frame
  ds_Shapiro <- data.frame()      #initialise the data frame for shapiro Therapist randomisation
  ds_Flex_Shapiro <- data.frame()      #initialise the data frame for shapiro Therapist randomisation
  ds_IBD <- data.frame()      #initialise the final data frame for Incomplete-Block Design
  
  N_Batch_values_Non<- vector("list", n_C)
  N_Batch_values_orthogonal<- vector("list", n_C)  #initialise the final batch sizes for Orthogonal design
  N_Batch_values<-vector("list", n_C)   #initialise the final batch sizes for Flexible design
  N_Batch_values_Shapiro<-vector("list", n_C) #initialise the final batch sizes for Shapiro Design
  N_Batch_values_Flex_Shapiro<-vector("list", n_C) #initialise the final batch sizes for Shapiro Design
  N_Batch_values_IBD<-vector("list", n_C)  #initialise the final batch sizes for IBD
  n_F<-rep(0,n_C)
  
  
  if(Add_Batch == 0){
    
    for(i in  1:n_C){   
      
      n_F[i]<- n_T[i]-n_P[i]  # number of Full-time therapists
      
      
      for(k in 1:n_B) {
        
        frac_P<-runif(1,min = a, max = 1) # fraction of part-time therapists available in a batch
        frac_F<-runif(1,min = a, max = 1) # fraction of full-time therapists available in a batch
        
        # Sample the available part and full-time therapists for current batch
        Therapist_Avail_Part<- sort(as.numeric(sample(Therapist_Part[[i]], size = ceiling(n_P[i]*frac_P), replace = FALSE)))
        Therapist_Avail_Full<- sort(as.numeric(sample(Therapist_Full[[i]], size = ceiling(n_F[i]*frac_F), replace = FALSE)))
        
        
        
        offset <- if (i > 1) sum(n_T[1:(i-1)]) else 0  # To adjust for the local indices of therapists 
        local_indices_part <- Therapist_Avail_Part - offset  # Local indices of part-time therapists
        local_indices_full <- Therapist_Avail_Full - (offset + n_P[i])  # Local indices of full-time therapists
        
        
        ## Generate data frame for non random therapist allocation
        
        df_Non_rand<- create_dataset_per_batch(Therapist_Avail_Part,Therapist_Avail_Full,
                                               n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full], sum(n_RP[[i]][local_indices_part])+ sum(n_RF[[i]][local_indices_full]) )
        
        # Function to generate Intervention
        ds1_Non_rand<- intervention_generation(df_Non_rand, n_I)
        
        # Add batch number to the data 
        ds1_Non_rand$Batch<- as.factor(rep(k,nrow(ds1_Non_rand)))
        # Add Centre number to the designs
        ds1_Non_rand$Centre<- as.factor(rep(i,nrow(ds1_Non_rand)))
        
        # A count of the batch size
        N_Batch_values_Non[[i]][k]<- nrow(ds1_Non_rand) # Batch sizes
        
        ds_Non_rand <- rbind(ds_Non_rand, ds1_Non_rand)    # Full data after selection of every batch
        
        #print("non random")
        #print(df_Non_rand)
        #print( ds1_Non_rand)
        
        # Generating treatment Randomisation and data frame 
        #treating IT combinations as treatments in the Orthogonal design
        
        # Baseline data structure for current batch
        
        df_orthogonal<- create_dataset_per_batch(Therapist_Part[[i]],Therapist_Full[[i]],
                                                 rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Part[[i]])),
                                                 rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Full[[i]])),min(unlist(c(n_RP, n_RF)), na.rm = TRUE)*n_T[i]*n_I)
        
        
        # randomise patients to IT combinations for this current batch
        ds1_orthogonal<- intervention_generation_randomised(df_orthogonal,n_I,Therapist_Part[[i]],Therapist_Full[[i]],
                                                            rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Part[[i]])),
                                                            rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Full[[i]])))
        # Add batch number to the data 
        ds1_orthogonal$Batch<- as.factor(rep(k,nrow(ds1_orthogonal)))
        # Add Centre number to the designs
        ds1_orthogonal$Centre<- as.factor(rep(i,nrow(ds1_orthogonal)))
        
        
        # A count of the batch size
        N_Batch_values_orthogonal[[i]][k]<- nrow(ds1_orthogonal) # Batch sizes
        
        ds_orthogonal <- rbind(ds_orthogonal, ds1_orthogonal)    # Full data after selection of every batch
        
        
        
        
        
        # Generating treatment Randomisation and data frame 
        #treating IT combinations as treatments 
        
        # Baseline data structure for current batch
        
        Result_Selected<- Create_Data_Therapist_Selected_From_Available(Therapist_Avail_Part,Therapist_Avail_Full,
                                                                        n_I,n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full],N_B[i],frac,local_indices_part,local_indices_full)
        df <- Result_Selected$df                        #Baseline Data frame
        Therapist_Select_Part<- Result_Selected$Therapist_Select_Part # Selected part-time therapists
        Therapist_Select_Full<- Result_Selected$Therapist_Select_Full # Selected full-time therapists
        n_RP_Selected<- Result_Selected$n_RP_Selected    # Selected part-time therapist's replication
        n_RF_Selected<- Result_Selected$n_RF_Selected    # Selected full-time therapist's replication
        N_Batch_Selected<- Result_Selected$N_Batch_Selected  #Selected batch size
        
        
        
        # randomise patients to IT combinations for this current batch
        ds1<- intervention_generation_randomised(df,n_I,Therapist_Select_Part,Therapist_Select_Full,
                                                 n_RP_Selected,n_RF_Selected)
        
        # Add batch number to the data 
        ds1$Batch<- as.factor(rep(k,nrow(ds1)))
        # Add Centre number to the designs
        ds1$Centre<- as.factor(rep(i,nrow(ds1)))
        
        N_Batch_values[[i]][k]<-  N_Batch_Selected # Batch sizes after selection
        
        
        
        ds <- rbind(ds, ds1)    # Full data after selection of every batch
        
        
        
        
        
        
        # Generating Therapist and Intervention randomisation separetly using Shapiro 1990 methods
        
        ds_Shapiro_T<- Randomise_Therapists_stage_1_Shapiro(Therapist_Avail_Part, Therapist_Avail_Full,
                                                            rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Avail_Part)),
                                                            rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Avail_Full)),min(unlist(c(n_RP, n_RF)), na.rm = TRUE)*n_T[i]*n_I)
        
        ds_Shapiro_I<-Randomise_Interventions_Stage_2(ds_Shapiro_T,n_I)
        
        #make factors for Therapists and Interventions
        ds_Shapiro_I$Therapist<-as.factor(ds_Shapiro_I$Therapist)
        ds_Shapiro_I$Int <-as.factor(ds_Shapiro_I$Int)
        
        
        # Add batch number to the data 
        ds_Shapiro_I$Batch<- as.factor(rep(k,nrow(ds_Shapiro_I)))
        # Add Centre number to the designs
        ds_Shapiro_I$Centre<- as.factor(rep(i,nrow(ds_Shapiro_I)))
        
        # Batch size for shapiro
        N_Batch_values_Shapiro[[i]][k]<- nrow(ds_Shapiro_I) # Batch sizes after selection
        
        
        
        ds_Shapiro <- rbind(ds_Shapiro, ds_Shapiro_I)  # Full data structure after 2-stage randomisation
        
        
        
        
        
        
        
        ### Generating Therapist and Intervention randomisation separetly using flexible Shapiro method
        
        ds_Flex_Shapiro_T<- Randomise_Therapists_stage_1_Shapiro(Therapist_Avail_Part, Therapist_Avail_Full,
                                                                 n_RP[[i]],n_RF[[i]], sum(2*c(n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full])) ) 
        ds_Flex_Shapiro_I<-Randomise_Interventions_Stage_2(ds_Flex_Shapiro_T,n_I)
        
        #make factors for Therapists and Interventions
        ds_Flex_Shapiro_I$Therapist<-as.factor(ds_Flex_Shapiro_I$Therapist)
        ds_Flex_Shapiro_I$Int <-as.factor(ds_Flex_Shapiro_I$Int)
        
        
        # Add batch number to the data 
        ds_Flex_Shapiro_I$Batch<- as.factor(rep(k,nrow(ds_Flex_Shapiro_I)))
        # Add Centre number to the designs
        ds_Flex_Shapiro_I$Centre<- as.factor(rep(i,nrow(ds_Flex_Shapiro_I)))
        
        # Batch size for shapiro
        N_Batch_values_Flex_Shapiro[[i]][k]<- nrow(ds_Flex_Shapiro_I) # Batch sizes after selection
        
        
        
        ds_Flex_Shapiro <- rbind(ds_Flex_Shapiro, ds_Flex_Shapiro_I)  # Full data structure after 2-stage randomisation
        
        
        
        
        # Generating treatment Randomisation using Incomplete-Block Design
        
        # Baseline data structure for current batch
        Result_Selected_IBD<- Randomise_Therapists_stage_1_IBD(Therapist_Avail_Part,Therapist_Avail_Full,
                                                               n_I, n_RP[[i]][local_indices_part], n_RF[[i]][local_indices_full], n_T_IBD[i])
        df_IBD <- Result_Selected_IBD$df_IBD               #Baseline Data frame using IBD
        Therapist_Select_Part<- Result_Selected_IBD$Therapist_Select_Part # Selected part-time therapists
        Therapist_Select_Full<- Result_Selected_IBD$Therapist_Select_Full # Selected full-time therapists
        n_RP_Selected<- Result_Selected_IBD$n_RP_Selected    # Selected part-time therapist's replication
        n_RF_Selected<- Result_Selected_IBD$n_RF_Selected    # Selected full-time therapist's replication
        N_Batch_IBD<- Result_Selected_IBD$N_Batch_IBD  #Selected batch size
        
        # randomise patients to IT combinations for this current batch on IBD
        #ds1_IBD<- intervention_generation_randomised(df_IBD,n_I,Therapist_Select_Part,Therapist_Select_Full,
        #                                             n_RP_Selected,n_RF_Selected)
        ds1_IBD<-  Randomise_Interventions_Stage_2(df_IBD,n_I)
        
        #make factors for Therapists and Interventions
        ds1_IBD$Therapist<-as.factor(ds1_IBD$Therapist)
        ds1_IBD$Int <-as.factor(ds1_IBD$Int)
        
        # Add batch number to the data 
        ds1_IBD$Batch<- as.factor(rep(k,nrow(ds1_IBD)))
        # Add Centre number to the designs
        ds1_IBD$Centre<- as.factor(rep(i,nrow(ds1_IBD)))
        
        N_Batch_values_IBD[[i]][k]<-  N_Batch_IBD # Batch sizes after selection
        
        
        
        ds_IBD <- rbind(ds_IBD, ds1_IBD)    # Full data after selection of every batch  
        
        
      } # for loop for batch ending
      
      
    }  # for loop for centre ending
    
    ds_Non_rand$Patient<- 1:nrow(ds_Non_rand)   # Renumber the patients
    ds_orthogonal$Patient<- 1:nrow(ds_orthogonal)   # Renumber the patients
    ds$Patient<- 1:nrow(ds)   # Renumber the patients
    ds_Shapiro$Patient<- 1:nrow(ds_Shapiro)   # Renumber the patients
    ds_Flex_Shapiro$Patient<- 1:nrow(ds_Flex_Shapiro)   # Renumber the patients
    ds_IBD$Patient<- 1:nrow(ds_IBD)   # Renumber the patients for IBD
    
  }   # if statement for Add_Batch = 0 ending
  
  
  if(Add_Batch == 1){
    
    k<-0   # Batch indicator
    #Stopping indicator for designs
    stop_Non_rand<-stop_orthognal<-stop_flex<-stop_Shapiro<-stop_Flex_Shapiro<-stop_IBD<-0
    
    # If under-recruitment in the Orthogonal Design, then add more batches
    while ( stop_Non_rand==0 || stop_orthognal==0 || stop_flex==0 ||  stop_Shapiro==0 ||  stop_Flex_Shapiro==0 || stop_IBD==0) {
      
      k<-k+1 # Batch indicator
      
      frac_P<-runif(1,min = a, max = 1) # fraction of part-time therapists available in a batch
      frac_F<-runif(1,min = a, max = 1) # fraction of full-time therapists available in a batch
      
      
      
      for(i in  1:n_C){   
        
        n_F[i]<- n_T[i]-n_P[i]  # number of Full-time therapists
        
        
        # Sample the available part and full-time therapists for current batch
        Therapist_Avail_Part<- sort(as.numeric(sample(Therapist_Part[[i]], size = ceiling(n_P[i]*frac_P), replace = FALSE)))
        Therapist_Avail_Full<- sort(as.numeric(sample(Therapist_Full[[i]], size = ceiling(n_F[i]*frac_F), replace = FALSE)))
        
        offset <- if (i > 1) sum(n_T[1:(i-1)]) else 0  # To adjust for the local indices of therapists 
        local_indices_part <- Therapist_Avail_Part - offset  # Local indices of part-time therapists
        local_indices_full <- Therapist_Avail_Full - (offset + n_P[i])  # Local indices of full-time therapists
        
        if(stop_Non_rand==0){       
          
          ## Generate data frame for non random therapist allocation
          
          df_Non_rand<- create_dataset_per_batch(Therapist_Avail_Part,Therapist_Avail_Full,
                                                 n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full], sum(n_RP[[i]][local_indices_part])+ sum(n_RF[[i]][local_indices_full]) )
          
          # Function to generate Intervention
          ds1_Non_rand<- intervention_generation(df_Non_rand, n_I)
          
          # Add batch number to the data 
          ds1_Non_rand$Batch<- as.factor(rep(k,nrow(ds1_Non_rand)))
          # Add Centre number to the designs
          ds1_Non_rand$Centre<- as.factor(rep(i,nrow(ds1_Non_rand)))
          
          # A count of the batch size
          N_Batch_values_Non[[i]][k]<- nrow(ds1_Non_rand) # Batch sizes
          
          ds_Non_rand <- rbind(ds_Non_rand, ds1_Non_rand)    # Full data after selection of every batch
          
          #print("non random")
          #print(df_Non_rand)
          #print( ds1_Non_rand)
          
          
        }
        
        if(stop_orthognal==0){
          # Generating treatment Randomisation and data frame 
          #treating IT combinations as treatments in the Orthogonal design
          
          # Baseline data structure for current batch
          
          df_orthogonal<- create_dataset_per_batch(Therapist_Part[[i]],Therapist_Full[[i]],
                                                   rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Part[[i]])),
                                                   rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Full[[i]])),min(unlist(c(n_RP, n_RF)), na.rm = TRUE)*n_T[i]*n_I)
          
          
          # randomise patients to IT combinations for this current batch
          ds1_orthogonal<- intervention_generation_randomised(df_orthogonal,n_I,Therapist_Part[[i]],Therapist_Full[[i]],
                                                              rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Part[[i]])),
                                                              rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Full[[i]])))
          # Add batch number to the data 
          ds1_orthogonal$Batch<- as.factor(rep(k,nrow(ds1_orthogonal)))
          # Add Centre number to the designs
          ds1_orthogonal$Centre<- as.factor(rep(i,nrow(ds1_orthogonal)))
          
          
          # A count of the batch size
          N_Batch_values_orthogonal[[i]][k]<- nrow(ds1_orthogonal) # Batch sizes
          
          ds_orthogonal <- rbind(ds_orthogonal, ds1_orthogonal)    # Full data after selection of every batch
          
          print("Walwyn design is done at this stage")
          
        }
        
        
        # If under-recruitment in the Flexible Design, then add more batches
        if(stop_flex == 0) {
          
          
          # Generating treatment Randomisation and data frame 
          #treating IT combinations as treatments 
          
          # Baseline data structure for current batch
          
          
          Result_Selected<- Create_Data_Therapist_Selected_From_Available(Therapist_Avail_Part,
                                                                          Therapist_Avail_Full,
                                                                          n_I,n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full],N_B[i],frac,local_indices_part,local_indices_full)
          df <- Result_Selected$df                        #Baseline Data frame
          Therapist_Select_Part<- Result_Selected$Therapist_Select_Part # Selected part-time therapists
          Therapist_Select_Full<- Result_Selected$Therapist_Select_Full # Selected full-time therapists
          n_RP_Selected<- Result_Selected$n_RP_Selected    # Selected part-time therapist's replication
          n_RF_Selected<- Result_Selected$n_RF_Selected    # Selected full-time therapist's replication
          N_Batch_Selected<- Result_Selected$N_Batch_Selected  #Selected batch size
          
          
          
          # randomise patients to IT combinations for this current batch
          ds1<- intervention_generation_randomised(df,n_I,Therapist_Select_Part,Therapist_Select_Full,
                                                   n_RP_Selected,n_RF_Selected)
          
          # Add batch number to the data 
          ds1$Batch<- as.factor(rep(k,nrow(ds1)))
          # Add Centre number to the designs
          ds1$Centre<- as.factor(rep(i,nrow(ds1)))
          
          N_Batch_values[[i]][k]<-  N_Batch_Selected # Batch sizes after selection
          
          
          
          ds <- rbind(ds, ds1)    # Full data after selection of every batch
          print("Flexible design is done at this stage")
          
          
        }
        
        
        # If under-recruitment in the Shapiro Design, then add more batches
        if(stop_Shapiro==0) {
          
          
          # Generating Therapist and Intervention randomisation separetly using Shapiro 1990 methods
          
          ds_Shapiro_T<- Randomise_Therapists_stage_1_Shapiro(Therapist_Avail_Part, Therapist_Avail_Full,
                                                              rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Avail_Part)),
                                                              rep(min(unlist(c(n_RP, n_RF)), na.rm = TRUE),length(Therapist_Avail_Full)),min(unlist(c(n_RP, n_RF)), na.rm = TRUE)*n_T[i]*n_I)
          ds_Shapiro_I<-Randomise_Interventions_Stage_2(ds_Shapiro_T,n_I)
          
          #make factors for Therapists and Interventions
          ds_Shapiro_I$Therapist<-as.factor(ds_Shapiro_I$Therapist)
          ds_Shapiro_I$Int <-as.factor(ds_Shapiro_I$Int)
          
          
          # Add batch number to the data 
          ds_Shapiro_I$Batch<- as.factor(rep(k,nrow(ds_Shapiro_I)))
          # Add Centre number to the designs
          ds_Shapiro_I$Centre<- as.factor(rep(i,nrow(ds_Shapiro_I)))
          
          # Batch size for shapiro
          N_Batch_values_Shapiro[[i]][k]<- nrow(ds_Shapiro_I) # Batch sizes after selection
          
          
          
          ds_Shapiro <- rbind(ds_Shapiro, ds_Shapiro_I)  # Full data structure after 2-stage randomisation
          
          print("Shapiro design is done at this stage")
          
        }
        
        # If under-recruitment in the Shapiro Design, then add more batches
        if(stop_Flex_Shapiro==0) {   
          ### Generating Therapist and Intervention randomisation separetly using flexible Shapiro method
          
          ds_Flex_Shapiro_T<- Randomise_Therapists_stage_1_Shapiro(Therapist_Avail_Part, Therapist_Avail_Full,
                                                                   n_RP[[i]],n_RF[[i]], sum(2*c(n_RP[[i]][local_indices_part],n_RF[[i]][local_indices_full])) ) 
          ds_Flex_Shapiro_I<-Randomise_Interventions_Stage_2(ds_Flex_Shapiro_T,n_I)
          
          #make factors for Therapists and Interventions
          ds_Flex_Shapiro_I$Therapist<-as.factor(ds_Flex_Shapiro_I$Therapist)
          ds_Flex_Shapiro_I$Int <-as.factor(ds_Flex_Shapiro_I$Int)
          
          
          # Add batch number to the data 
          ds_Flex_Shapiro_I$Batch<- as.factor(rep(k,nrow(ds_Flex_Shapiro_I)))
          # Add Centre number to the designs
          ds_Flex_Shapiro_I$Centre<- as.factor(rep(i,nrow(ds_Flex_Shapiro_I)))
          
          # Batch size for shapiro
          N_Batch_values_Flex_Shapiro[[i]][k]<- nrow(ds_Flex_Shapiro_I) # Batch sizes after selection
          
          
          
          ds_Flex_Shapiro <- rbind(ds_Flex_Shapiro, ds_Flex_Shapiro_I)  # Full data structure after 2-stage randomisation
          
          print("Flex.Shapiro design is done at this stage")
          
        }  
        
        
        # If under-recruitment in the Shapiro Design, then add more batches
        if( stop_IBD==0) {
          
          
          # Generating treatment Randomisation and data frame 
          #treating IT combinations as treatments using Incomplete-Block Design
          
          # Baseline data structure for current batch
          
          Result_Selected_IBD<- Randomise_Therapists_stage_1_IBD(Therapist_Avail_Part,Therapist_Avail_Full,
                                                                 n_I, n_RP[[i]][local_indices_part], n_RF[[i]][local_indices_full], n_T_IBD[i])
          df_IBD <- Result_Selected_IBD$df_IBD               #Baseline Data frame using IBD
          Therapist_Select_Part<- Result_Selected_IBD$Therapist_Select_Part # Selected part-time therapists
          Therapist_Select_Full<- Result_Selected_IBD$Therapist_Select_Full # Selected full-time therapists
          n_RP_Selected<- Result_Selected_IBD$n_RP_Selected    # Selected part-time therapist's replication
          n_RF_Selected<- Result_Selected_IBD$n_RF_Selected    # Selected full-time therapist's replication
          N_Batch_IBD<- Result_Selected_IBD$N_Batch_IBD  #Selected batch size
          
          # randomise patients to IT combinations for this current batch on IBD
          #ds1_IBD<- intervention_generation_randomised(df_IBD,n_I,Therapist_Select_Part,Therapist_Select_Full,
          #                                             n_RP_Selected,n_RF_Selected)
          ds1_IBD<-  Randomise_Interventions_Stage_2(df_IBD,n_I)
          
          #make factors for Therapists and Interventions
          ds1_IBD$Therapist<-as.factor(ds1_IBD$Therapist)
          ds1_IBD$Int <-as.factor(ds1_IBD$Int)
          
          # Add batch number to the data 
          ds1_IBD$Batch<- as.factor(rep(k,nrow(ds1_IBD)))
          # Add Centre number to the designs
          ds1_IBD$Centre<- as.factor(rep(i,nrow(ds1_IBD)))
          
          N_Batch_values_IBD[[i]][k]<-  N_Batch_IBD # Batch sizes after selection
          
          
          
          ds_IBD <- rbind(ds_IBD, ds1_IBD)    # Full data after selection of every batch
          print("IBD design is done at this stage")
          
        }
        
      } # Centre loop end
      
      if(nrow(ds_Non_rand)>= N){stop_Non_rand<-1}
      if(nrow(ds_orthogonal)>= N){stop_orthognal<-1}
      if(nrow(ds)>= N){stop_flex<-1}
      if(nrow(ds_Shapiro)>= N ){stop_Shapiro<-1}
      if(nrow(ds_Flex_Shapiro)>= N ){stop_Flex_Shapiro<-1}
      if(nrow(ds_IBD)>= N){stop_IBD<-1}
      
    } # While loop end
    
    ds_Non_rand$Patient<- 1:nrow(ds_Non_rand)   # Renumber the patients
    ds_orthogonal$Patient<- 1:nrow(ds_orthogonal)   # Renumber the patients
    ds$Patient<- 1:nrow(ds)   # Renumber the patients
    ds_Shapiro$Patient<- 1:nrow(ds_Shapiro)   # Renumber the patients
    ds_Flex_Shapiro$Patient<- 1:nrow(ds_Flex_Shapiro)   # Renumber the patients
    ds_IBD$Patient<- 1:nrow(ds_IBD)   # Renumber the patients for IBD
    
  } # if condition end 
  
  
  res <- list(
    ds_Non_rand = ds_Non_rand,
    ds_orthogonal = ds_orthogonal,
    ds_Shapiro   = ds_Shapiro,
    ds =  ds,
    ds_Flex_Shapiro = ds_Flex_Shapiro,
    ds_IBD = ds_IBD,
    N_Batch_values_Non = N_Batch_values_Non,
    N_Batch_values_orthogonal = N_Batch_values_orthogonal,
    N_Batch_values = N_Batch_values,
    N_Batch_values_Shapiro = N_Batch_values_Shapiro,
    N_Batch_values_Flex_Shapiro = N_Batch_values_Flex_Shapiro,
    N_Batch_values_IBD = N_Batch_values_IBD
  )
  
}

########################################################################################################################################
############################################################################################################
#####################################################################

# Model Fitting
Analysis_singularity <- function(
    ds,
    del_0, del_1, del_2, del_3, del_4,
    sig2_T, sig2_B, sig2_C,
    sig2_IT, sig2_IB, sig2_IC,
    sig2_TB, sig2_CB,
    sig2_ITB, sig2_ICB,
    sig2_T_learn, sig2_C_learn,
    sig2_IT_learn, sig2_IC_learn,
    sig2_E,
    time, matching,
    tol = 1e-4
) {
  
  ## ---------------------------------------------------------
  ## Generate data
  ## ---------------------------------------------------------
  dat <- samples_from_Design(
    ds,
    del_0, del_1, del_2, del_3, del_4,
    sig2_T, sig2_B, sig2_C,
    sig2_IT, sig2_IB, sig2_IC,
    sig2_TB, sig2_CB,
    sig2_ITB, sig2_ICB,
    sig2_T_learn, sig2_C_learn,
    sig2_IT_learn, sig2_IC_learn,
    sig2_E,
    time, matching
  )
  
  if (all(is.na(dat$y))) {
    # Return NA vector of correct length depending on n_C
    n_C <- length(unique(dat$Centre))
    vc_names <- if (n_C > 1) {
      c(
        "Therapist","Batch","contrast:Therapist","Therapist:Batch",
        "contrast:Batch","contrast:Therapist:Batch",
        "Centre","Centre:Batch","contrast:Centre","contrast:Centre:Batch"
      )
    } else {
      c(
        "Therapist","Batch","contrast:Therapist","Therapist:Batch",
        "contrast:Batch","contrast:Therapist:Batch"
      )
    }
    return(setNames(rep(NA_real_, length(vc_names)), vc_names))
  }
  
  n_C <- length(unique(dat$Centre))
  
  ## ---------------------------------------------------------
  ## Model formulas same as before
  ## ---------------------------------------------------------
  if (n_C > 1) {
    vc_names <- c(
      "Therapist","Batch","Centre", "contrast:Therapist","contrast:Batch",
      "contrast:Centre","Therapist:Batch","Centre:Batch",
      "contrast:Therapist:Batch","contrast:Centre:Batch"
    )
    
    form <- y ~ contrast +
      (1 | Therapist) +
      (1 | Batch) +
      (1 | Therapist:Batch) +
      (1 | Centre) +
      (1 | Centre:Batch) +
      (1 | contrast:Therapist) +
      (1 | contrast:Batch) +
      (1 | contrast:Therapist:Batch) +
      (1 | contrast:Centre) +
      (1 | contrast:Centre:Batch)
    
  } else {
    vc_names <- c(
      "Therapist","Batch","contrast:Therapist", "contrast:Batch",
      "Therapist:Batch", "contrast:Therapist:Batch"
    )
    
    form <- y ~ contrast +
      (1 | Therapist) +
      (1 | Batch) +
      (1 | Therapist:Batch) +
      (1 | contrast:Therapist) +
      (1 | contrast:Batch) +
      (1 | contrast:Therapist:Batch)
  }
  
  ## Template NA output
  out <- setNames(rep(NA_real_, length(vc_names)), vc_names)
  
  ## ---------------------------------------------------------
  ## Fit model safely
  ## ---------------------------------------------------------
  fit <- tryCatch(
    lmer(form, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(out)
  
  ## ---------------------------------------------------------
  ## Extract variance components & check singularity
  ## ---------------------------------------------------------
  vc <- VarCorr(fit)
  
  for (i in seq_along(vc_names)) {
    tmp <- vc[[vc_names[i]]]
    if (!is.null(tmp)) {
      # Check singularity: variance too small
      out[i] <- as.numeric(attr(tmp, "stddev")^2 < tol)
    }
  }
  
  print(summary(fit))
  return(out)
}



#####################################################################################################################################
#####################################################################################################################################
#####PARAMETERS FOR RECRUITMENT


N <- 1200  # Required sample size in the trial known from Pilot studies

n_C <- 15       # number of Centres (1 for single centre trial)

n_B <- 5        # number of batches
n_I <- 2       # number of interventions

E_n_T <- 4        ## Expected number of therapists in a centre for a multi-centre trial
n_T_single <- 4  ## number of therapists for a single-centre trial

Ex_C_P <-15         #expected capacity of a part-timer
Ex_C_F <-20         #expected capacity of a full-timer
sd_P <- 6           #SD of capacity of a part-timer 
sd_F <- 2           #SD of capacity of a full-timer

a_P <- 0           ### proportion of part-timers in a centre
a <- 1             ### parameter for proportion of available therapists. Distribution Uniform[a,1].

frac <- 0.25  # proportion of extra recruitment in a batch permitted
Add_Batch <-1  # Indicator of whether to add new batches or not if under-recruitment. 1- yes/ 0- No

n_sim <- 25      # Number of simulation runs
time <- 0       # Simulate data from : 0 - Block effect model, 1 - Learning effect model 



#####################################################################################################################################
#####################################################################################################################################
##### PARAMETERS FOR ANALYSIS


if(n_C == 1){
  ### Fixed effects parameters
  del_0<- 0        # delta_0 =  The overall mean effect
  del_1<- 0.4      # Treatment effect
  del_2<- 0.1      #Time effect
  del_3<- 0.1      #Time by intervention interaction effect  
  del_4<- 0.3      # Covariate for patient matching
  
  ### Random effects parameters
  sig2_T<- 0.1      # Sigma for Therapist effects per centre
  sig2_B<-0.1      # Sigma for Batch effect
  sig2_C<-0      # Sigma for Centre effect
  sig2_IT<-0.15    # Sigma for Therapist*Intervention effect
  sig2_IB<-0.1     # Sigma for Intervention*Batch effect
  sig2_IC<-0      # Sigma for Intervention*Centre effect
  sig2_TB<-0.1     # Sigma for Therapist*Batch effect
  sig2_CB<-0      # Sigma for Batch*Centre effect
  sig2_ITB<-0.1   # Sigma for Intervention*Therapist*Batch effect
  sig2_ICB<-0      # Sigma for Intervention*Batch*Centre effect
  sig2_T_learn<-0.1     # Sigma for Therapist learning effect
  sig2_C_learn<-0     # Sigma for Centre learning effect
  sig2_IT_learn<-0.1   # Sigma for Intervention*Therapist learning effect
  sig2_IC_learn<-0    # Sigma for Intervention*Centre learning effect
  sig2_E<- 0.35    # Sigma for patient specific error
  
}else{
  ### Fixed effects parameters
  del_0<- 0        # delta_0 =  The overall mean effect
  del_1<- 0.4      # Treatment effect
  del_2<- 0.1      #Time effect
  del_3<- 0.1      #Time by intervention interaction effect  
  del_4<- 0.3      # Covariate for patient matching
  
  ### Random effects parameters
  sig2_T<- rep(0.1,n_C)  # c(0.1,0.05,0.02,0.08,0.1,0.05)      # Sigma for Therapist effects per centre
  sig2_B<-0.1      # Sigma for Batch effect
  sig2_C<-0.1      # Sigma for Centre effect
  sig2_IT<-0.1    # Sigma for Therapist*Intervention effect
  sig2_IB<-0.1     # Sigma for Intervention*Batch effect
  sig2_IC<-0.1      # Sigma for Intervention*Centre effect
  sig2_TB<-0.05     # Sigma for Therapist*Batch effect
  sig2_CB<-0.05      # Sigma for Batch*Centre effect
  sig2_ITB<-0.05   # Sigma for Intervention*Therapist*Batch effect
  sig2_ICB<-0.05      # Sigma for Intervention*Batch*Centre effect
  sig2_T_learn<-0.05     # Sigma for Therapist learning effect
  sig2_C_learn<-0.05     # Sigma for Centre learning effect
  sig2_IT_learn<-0.05    # Sigma for Intervention*Therapist learning effect
  sig2_IC_learn<-0.05    # Sigma for Intervention*Centre learning effect
  sig2_E<- 0.2    # Sigma for patient specific error
  
}


####################################################################################################################################
#####################################################################################################################################
# --------------------------
# Storage for Standard Model
# --------------------------


# ------------------------------------------------------
# Storage
# ------------------------------------------------------
design_names <- c("Non.rand", "Wal25", "Flexible", "SPP2", "IBD", "Flex.SPP2")

singularity_store <- vector("list", length(design_names))
names(singularity_store) <- design_names

for (d in design_names) {
  singularity_store[[d]] <- list()
}

# ------------------------------------------------------
# Simulation loop
# ------------------------------------------------------
for (j in seq_len(n_sim)) {
  
  ## ===============================
  ## 1. Generate therapist structure
  ## ===============================
  Therapist_Part <- list()
  Therapist_Full <- list()
  c_P  <- list()
  c_F  <- list()
  n_RP <- list()
  n_RF <- list()
  
  if (n_C > 1) {
    n_T <- rpois(n_C, lambda = E_n_T)
    n_T[n_T == 0] <- 1
  } else {
    n_T <- n_T_single
  }
  
  n_P <- floor(n_T * a_P)
  
  ## Therapist IDs
  for (i in seq_along(n_T)) {
    offset <- if (i > 1) sum(n_T[1:(i - 1)]) else 0
    
    Therapist_Part[[i]] <-
      if (n_P[i] > 0) as.character((offset + 1):(offset + n_P[i])) else character(0)
    
    Therapist_Full[[i]] <-
      if (n_T[i] - n_P[i] > 0)
        as.character((offset + n_P[i] + 1):(offset + n_T[i]))
    else character(0)
  }
  
  ## Capacities
  for (i in seq_len(n_C)) {
    
    c_P[[i]] <-
      if (n_P[i] > 0)
        round(rtruncnorm(n_P[i], 0, Inf, Ex_C_P, sd_P))
    else numeric(0)
    
    c_F[[i]] <-
      if (n_T[i] - n_P[i] > 0)
        round(rtruncnorm(n_T[i] - n_P[i], 0, Inf, Ex_C_F, sd_F))
    else numeric(0)
    
    n_RP[[i]] <- pmax(round(c_P[[i]] / (2 * n_B)), 1)
    n_RF[[i]] <- pmax(round(c_F[[i]] / (2 * n_B)), 1)
  }
  
  n_T_IBD <- pmax(floor(n_T * a), 1)
  
  ## ===============================
  ## 2. Generate designs (CRITICAL)
  ## ===============================
  designs <- generate_designs(
    N, n_B, n_C, n_I, n_T, n_P,
    Therapist_Part, Therapist_Full,
    n_T_IBD, n_RP, n_RF,
    frac, a, Add_Batch, Ex_C_P, Ex_C_F
  )
  
  datasets <- list(
    Non_rand  = designs$ds_Non_rand,
    Wal25     = designs$ds_orthogonal,
    Flexible  = designs$ds,
    SPP2      = designs$ds_Shapiro,
    IBD       = designs$ds_IBD,
    Flex_SPP2 = designs$ds_Flex_Shapiro
  )
  
  ## ===============================
  ## 3. Fit model for each design
  ## ===============================
  for (i in seq_along(datasets)) {
    
    matching <- ifelse(i == 1, 1, 0)
    
    out <- Analysis_singularity(
      ds = datasets[[i]],
      del_0, del_1, del_2, del_3, del_4,
      sig2_T, sig2_B, sig2_C,
      sig2_IT, sig2_IB, sig2_IC,
      sig2_TB, sig2_CB,
      sig2_ITB, sig2_ICB,
      sig2_T_learn, sig2_C_learn,
      sig2_IT_learn, sig2_IC_learn,
      sig2_E,
      time, matching
    )
    
    singularity_store[[i]][[j]] <- out
  }
}

# ------------------------------------------------------
# Aggregate: proportions (drop full NA sims)
# ------------------------------------------------------
singularity_results <- lapply(singularity_store, function(x) {
  
  mat <- do.call(rbind, x)
  
  keep <- !apply(is.na(mat), 1, all)
  mat  <- mat[keep, , drop = FALSE]
  
  colMeans(mat, na.rm = TRUE)
})

# ------------------------------------------------------
# Final table
# ------------------------------------------------------
# ------------------------------------------------------
# Final numeric matrix (NO data.frame, NO characters)
# ------------------------------------------------------
combined_singularity <- do.call(rbind, singularity_results)

# attach design names properly
rownames(combined_singularity) <- names(singularity_results)

# enforce numeric (safety, not strictly needed)
storage.mode(combined_singularity) <- "double"

print(round(combined_singularity, 3))

saveRDS(
  combined_singularity,
  file = paste0("combined_singularity_", seed, ".rds")
)
