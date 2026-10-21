
library(ggplot2)

library(abind)  # for combining arrays


setwd("C:/Users/fvsk236/Desktop/Results_Therapist_Paper/HPC files/")

### Combine estimation results

files <- paste0("combined_standard_", c(1001:1405), ".rds")

# Keep only existing files
existing_files <- files[file.exists(files)]

cat("Found", length(existing_files), "files.\n")
all_results <- lapply(existing_files, readRDS)
str(all_results[[1]])


all_array <- abind::abind(all_results, along = 3)  # 3rd dimension = file index

# Compute the cell-wise average WITHOUT letting NA explode the result
avg_result <- apply(all_array, c(1, 2), function(x) mean(x, na.rm = TRUE))

# Convert to data.frame
avg_df <- as.data.frame(avg_result)

# Identify columns for 3-decimal formatting
cols_3 <- c("m", "s")

# Apply rounding
avg_df[, cols_3] <- lapply(avg_df[, cols_3, drop = FALSE], function(x) round(x, 3))

# All other columns get 2 decimals
cols_2 <- setdiff(colnames(avg_df), cols_3)
avg_df[, cols_2] <- lapply(avg_df[, cols_2, drop = FALSE], function(x) round(x, 2))


   ### ============================
   ### BUILD LATEX TABLE OUTPUT
   ### ============================
   
# Replace NaN or NA with "-"
# ------------------------
avg_df[is.na(avg_df)] <- "-"

### ============================
### BUILD LATEX TABLE OUTPUT
### ============================

# Number of columns in part 1
k <- 8  

# ---------- PART 1 ----------
avg_df_1 <- avg_df[, 1:k, drop = FALSE]

body_lines_1 <- lapply(1:nrow(avg_df_1), function(i) {
  row_name <- rownames(avg_df_1)[i]
  if (row_name == "True Values") return(NULL)  # skip this row
  paste0("&& ", paste(c(row_name, avg_df_1[i, ]), collapse = " & "), " \\\\")
})
body_lines_1 <- unlist(body_lines_1)


# ---------- PART 2 ----------
avg_df_2 <- avg_df[, (k + 1):ncol(avg_df), drop = FALSE]

body_lines_2 <- lapply(1:nrow(avg_df_2), function(i) {
  row_name <- rownames(avg_df_2)[i]
  if (row_name == "True Values") return(NULL)  # skip this row
  paste0("&& ", paste(c(row_name, avg_df_2[i, ]), collapse = " & "), " \\\\")
})
body_lines_2 <- unlist(body_lines_2)

# Final table (for inspection)
avg_df

# ---------- PRINT ----------
cat(body_lines_1, sep = "\n")

cat("\n\n")  # space between the two parts

cat(body_lines_2, sep = "\n")
################################################################################################################################
###############################################   ###############################################
### Combine time-to-target results
   
   files <- paste0("time_to_target_", 1001:1405, ".rds")
   
   # Keep only existing files
   existing_files <- files[file.exists(files)]
   
   cat("Found", length(existing_files), "files.\n")

   # Read all files; each must be a matrix of same dimension
   all_results <- lapply(existing_files, readRDS)
   
   # Convert list → 3D array where dim = [rows, cols, files]
   arr <- simplify2array(all_results)

   # Print
   print( round(rowMeans(arr),2))
   
   
################################################################################################################################
   ###############################################   ###############################################
### Plot and compare RMSEs
   
   ##### Compare the RMSE per simulations and then check their statistical properties
   
   ## Store in a variable
   
   files <- paste0("MSE_marg_std_", c(1001:1405), ".rds")
   
   # Keep only existing files
   existing_files <- files[file.exists(files)]
   
   cat("Found", length(existing_files), "files.\n")
   
   # Read each RDS and combine them into one data frame
   MSE_all <- do.call(rbind, lapply(existing_files, function(f) {
     df <- readRDS(f)
     # Make sure it's a data frame
     if(!is.data.frame(df)) df <- as.data.frame(df)
     df
   }))   
   
   
   compare_pair <- function(a, b, name_a="A", name_b="B") {
     
     # keep only complete pairs
     idx <- complete.cases(a, b)
     if (!any(idx)) {
       stop("All paired values are NA. Nothing to compare.")
     }
     
     a_clean <- a[idx]
     b_clean <- b[idx]
     d <- a_clean - b_clean
     
     res <- list(
       mean_d   = mean(d),
       sd_d     = sd(d),
       median_d = median(d),
       wins_a   = mean(a_clean < b_clean),
       t_test   = t.test(d, mu = 0, conf.level = 1 - 0.05/10),
       wilcox   = wilcox.test(d, mu = 0, exact = FALSE)
     )
     
     # bootstrap CI
     library(boot)
     boot_mean <- function(x, i) mean(x[i])
     b <- boot(d, statistic = boot_mean, R = 2000)
     res$boot_ci <- boot.ci(b, type = c("perc","bca"))$percent[4:5]
     
     return(res)
   }
   
   

   
   
   cp1 <- compare_pair(MSE_all[,2], MSE_all[,3], name_a="Walwyn", name_b="Flexible")
   cp2 <- compare_pair(MSE_all[,2], MSE_all[,4], name_a="Walwyn", name_b="Shapiro")
   cp3 <- compare_pair(MSE_all[,2], MSE_all[,5], name_a="Walwyn", name_b="IBD")
   cp4 <- compare_pair(MSE_all[,2], MSE_all[,6], name_a="Walwyn", name_b="Flex.SPP2iro")
   cp5 <- compare_pair(MSE_all[,3], MSE_all[,4], name_a="Flexible", name_b="Shapiro")
   cp6 <- compare_pair(MSE_all[,3], MSE_all[,5], name_a="Flexible", name_b="IBD")
   cp7 <- compare_pair(MSE_all[,3], MSE_all[,6], name_a="Flexible", name_b="Flex.SPP2iro")
   cp8 <- compare_pair(MSE_all[,4], MSE_all[,5], name_a="Shapiro", name_b="IBD")
   cp9 <- compare_pair(MSE_all[,4], MSE_all[,6], name_a="Shapiro", name_b="Flex.SPP2iro")
   cp10 <- compare_pair(MSE_all[,5], MSE_all[,6], name_a="IBD", name_b="Flex.SPP2iro")
   
   df <- data.frame(
     Comparison = c("Wal25 vs Flex","Wal25 vs SPP2","Wal25 vs IBD","Wal25 vs Flex.SPP2",
                    "Flex vs SPP2", "Flex vs IBD","Flex vs Flex.SPP2",
                    "SPP2 vs IBD","SPP2 vs Flex.SPP2","IBD vs Flex.SPP2"),
     Mean = c(cp1$t_test$estimate, cp2$t_test$estimate, cp3$t_test$estimate,cp4$t_test$estimate,cp5$t_test$estimate,cp6$t_test$estimate,
              cp7$t_test$estimate,cp8$t_test$estimate,cp9$t_test$estimate,cp10$t_test$estimate),
     Lower = c(cp1$t_test$conf.int[1], cp2$t_test$conf.int[1], cp3$t_test$conf.int[1], cp4$t_test$conf.int[1], cp5$t_test$conf.int[1], cp6$t_test$conf.int[1],
               cp7$t_test$conf.int[1],cp8$t_test$conf.int[1],cp9$t_test$conf.int[1],cp10$t_test$conf.int[1]),
     Upper = c(cp1$t_test$conf.int[2], cp2$t_test$conf.int[2], cp3$t_test$conf.int[2], cp4$t_test$conf.int[2], cp5$t_test$conf.int[2], cp6$t_test$conf.int[2],
               cp7$t_test$conf.int[2], cp8$t_test$conf.int[2], cp9$t_test$conf.int[2], cp10$t_test$conf.int[2])
   )
   
   
   df <- data.frame(
     Comparison = c("Flex vs SPP2", "Flex vs IBD","Flex vs Flex.SPP2",
                    "SPP2 vs IBD","SPP2 vs Flex.SPP2","IBD vs Flex.SPP2"),
     Mean = c( cp5$t_test$estimate, cp6$t_test$estimate,cp7$t_test$estimate,
              cp8$t_test$estimate,cp9$t_test$estimate,cp10$t_test$estimate),
     Lower = c( cp5$t_test$conf.int[1], cp6$t_test$conf.int[1], cp7$t_test$conf.int[1],
               cp8$t_test$conf.int[1],cp9$t_test$conf.int[1],cp10$t_test$conf.int[1]),
     Upper = c( cp5$t_test$conf.int[2], cp6$t_test$conf.int[2], cp7$t_test$conf.int[2], 
                cp8$t_test$conf.int[2], cp9$t_test$conf.int[2], cp10$t_test$conf.int[2])
   )  
   
   
  # If all Walwyn and shapiro results are NA
   
   cp6 <- compare_pair(MSE_all[,3], MSE_all[,5], name_a="Flexible", name_b="IBD")
   cp7 <- compare_pair(MSE_all[,3], MSE_all[,6], name_a="Flexible", name_b="Flex.SPP2iro")
   cp10 <- compare_pair(MSE_all[,5], MSE_all[,6], name_a="IBD", name_b="Flex.SPP2iro")
   
   
  
   df <- data.frame(
     Comparison = c( "Flex vs IBD","Flex vs Flex.SPP2","IBD vs Flex.SPP2"),
     Mean = c(  cp6$t_test$estimate,cp7$t_test$estimate,cp10$t_test$estimate),
     Lower = c( cp6$t_test$conf.int[1], cp7$t_test$conf.int[1],cp10$t_test$conf.int[1]),
   Upper = c( cp6$t_test$conf.int[2], cp7$t_test$conf.int[2], cp10$t_test$conf.int[2])
   )
   
   #### Plot the final figure
   
   ggplot(df, aes(x = Comparison, y = Mean)) +
     geom_point(size = 3) +
     geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, lwd = 1) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     labs(
       title = "Model-I",
       y = "RMSE Difference ± 95% CI",
       x = ""
     ) +
     theme_minimal(base_size = 12) +  # smaller base font
     theme(
       axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # smaller x-axis ticks
       axis.text.y = element_text(size = 10),                           # smaller y-axis ticks
       axis.title.y = element_text(size = 12),
       plot.title = element_text(size = 14, face = "bold")              # optional: smaller title
     )
   
   
   ##############################################################
   ##############################################################################
   #####################################################################
   ## Combine singularities##############################
   

   # ------------------------------------------------------
   # 1. Files
   # ------------------------------------------------------
   files <- paste0("combined_singularity_", 1001:1405, ".rds")
   existing_files <- files[file.exists(files)]
   cat("Found", length(existing_files), "files.\n")
   
   # ------------------------------------------------------
   # 2. Read all files
   # ------------------------------------------------------
   all_results <- lapply(existing_files, readRDS)
   
   # ------------------------------------------------------
   # 3. Combine into 3D array
   # ------------------------------------------------------
   all_array <- abind::abind(all_results, along = 3)
   
   # ------------------------------------------------------
   # 4. Compute cell-wise average ignoring NA
   # ------------------------------------------------------
   avg_result <- apply(all_array, c(1,2), function(x) mean(x, na.rm = TRUE))
   
   # ------------------------------------------------------
   # 5. Convert to numeric matrix
   # ------------------------------------------------------
   avg_df <- matrix(avg_result, nrow = nrow(avg_result), ncol = ncol(avg_result))
   rownames(avg_df) <- rownames(all_results[[1]])
   colnames(avg_df) <- colnames(all_results[[1]])
   
   # ------------------------------------------------------
   # 6. Round numeric values (3 decimals)
   # ------------------------------------------------------
   avg_df <- round(avg_df, 2)
   
   # ------------------------------------------------------
   # 7. Convert to LaTeX-ready strings with '&'
   # ------------------------------------------------------
   latex_table <- apply(avg_df, 1, function(x) paste(x, collapse = " & "))
   names(latex_table) <- paste0("&& ", rownames(avg_df))   
   # ------------------------------------------------------
   # 8. Print
   # ------------------------------------------------------
   
   print(avg_df)
   
   for (i in seq_along(latex_table)) {
     cat(names(latex_table)[i], "&", latex_table[i], "\\\\\n")
   }
   
   
  #####################################################################################################################################
   ##############################################################################################################################
############################## Box plots of RMSE vs Number of Batches###################################
   
   library(ggplot2)
   library(reshape2)
   library(dplyr)
   
   setwd("C:/Users/fvsk236/Desktop/Results_Therapist_Paper/HPC files/")

   ####  Combine MSE values from RDS
   mse_files <- paste0("MSE_marg_std_", 1001:1405, ".rds")
   existing_mse <- mse_files[file.exists(mse_files)]
   cat("Found", length(existing_mse), "MSE files.\n")
   
   # Read and stack
   MSE_marg_std <- do.call(rbind, lapply(existing_mse, function(f) {
     df <- readRDS(f)
     if(!is.data.frame(df)) df <- as.data.frame(df)
     df
   }))
   
   ### Combine time-to-target values
   time_files <- paste0("time_to_target_", 1001:1405, ".rds")
   existing_time <- time_files[file.exists(time_files)]
   cat("Found", length(existing_time), "time-to-target files.\n")
   
   all_time <- lapply(existing_time, readRDS)
   arr <- simplify2array(all_time)
   
   # Average number of batches per design
   time_to_target_avg <- round(rowMeans(arr), 2)
   
   # Attach names to match design columns
   names(time_to_target_avg) <- colnames(MSE_marg_std)
   
   #### Reshape MSE to long format
   df <- melt(MSE_marg_std, variable.name = "Design", value.name = "MSE")
   
   # Add row index
   df$Index <- rep(1:nrow(MSE_marg_std), times = ncol(MSE_marg_std))
   
   # Reorder columns
   df <- df[, c("Index", "Design", "MSE")]
   
   # Attach batch info safely
   df$time_to_target_avg <- time_to_target_avg[as.character(df$Design)]
   
   # Remove NAs
   df <- df %>% filter(!is.na(MSE), !is.na(time_to_target_avg))
   
   design_cols <- c(
     "Non.rand"   = "#F8766D",
     "Wal25"      = "#A6761D",
     "Flexible"   = "#B79F00",
     "SPP2"       = "#00BA38",
     "IBD"        = "#00B0F6",
     "Flex.SPP2"  = "#E76BF3"
   )
   
   df$Design <- factor(df$Design, levels = names(design_cols))
   
   ## Optional: remove Wal25 safely
   #df <- df %>% filter(Design != "Wal25")
   #df$Design <- droplevels(df$Design)
   
   
   ## Compute Y-axis limits using box statistics
   box_stats <- df %>%
     group_by(time_to_target_avg, Design) %>%
     summarise(
       Q1 = quantile(MSE, 0.25),
       Q3 = quantile(MSE, 0.75),
       IQR = IQR(MSE),
       .groups = "drop"
     ) %>%
     mutate(upper_whisker = Q3 + 1.5 * IQR)
   
   y_max <- ceiling(max(box_stats$upper_whisker))
   y_min <- floor(min(df$MSE))
   
   ## Final plot
   p_box <- ggplot(
     df,
     aes(x = factor(time_to_target_avg), y = MSE, fill = Design)
   ) +
     geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
     stat_summary(
       fun = mean, geom = "point", shape = 18, size = 3,
       color = "black", show.legend = FALSE
     ) +
     scale_fill_manual(values = design_cols, drop = TRUE) +
     scale_y_continuous(limits = c(y_min, y_max), breaks = pretty(c(y_min, y_max), 5)) +
     labs(x = "Recruitment Time to Target", y = "RMSE", fill = "Design") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   print(p_box)
   
 
   
   
   
   
     
   # save as ggplot file for later editing
   saveRDS(
     p_box,
     file = "C:/Users/fvsk236/Desktop/Results_Therapist_Paper/Figures/N=240/Model_1_nT_6_nB_5_aP_0/Model_1_nT_6_nB_5_aP_0_a_1"
   )   
   
   
   ### load and edit later the ggplot
   
   p_box <- readRDS("C:/Users/fvsk236/Desktop/Results_Therapist_Paper/Figures/N=240/Model_1_nT_12_nB_5_aP_1/Model_1_nT_12_nB_5_aP_1_a_1")
   
   # Remove old fill scales
   p_box$scales$scales <- lapply(p_box$scales$scales, function(s) {
     if ("fill" %in% s$aesthetics) return(NULL) else return(s)
   })
   p_box$scales$scales <- Filter(Negate(is.null), p_box$scales$scales)
   
   # Add new fill scale and update label
   p_box <- p_box +
     scale_fill_manual(values = design_cols, drop = TRUE) +
     labs(x = "Recruitment Time to Target")
   
   print(p_box)
   
   saveRDS(
     p_box,
     file = "C:/Users/fvsk236/Desktop/Results_Therapist_Paper/Figures/N=240/Model_1_nT_12_nB_5_aP_1/Model_1_nT_12_nB_5_aP_1_a_1"
   )   
   
   
   
   
   
   
   ####### #########################################################################################################
   ####################################################################################################################
   ##################################################################################
  #################### Combine the power results
   
   ####  Combine power values from RDS
   # 1. Get the files
   power_files <- paste0("power_", 1001:1405, ".rds")
   existing_power <- power_files[file.exists(power_files)]
   cat("Found", length(existing_power), "power files.\n")
   
   # 2. Read them in
   all_matrices <- lapply(existing_power, readRDS)
   
   # 3. Convert list-of-lists-of-vectors to 3D array
   # Determine dimensions
   n_rows <- length(all_matrices[[1]])            # number of "rows" = length of outer list element
   n_cols <- length(all_matrices[[1]][[1]])      # number of "columns" = length of inner vectors
   n_files <- length(all_matrices)               # 3rd dim = number of files / seeds
   
   # Initialize array
   mat_array <- array(NA, dim = c(n_rows, n_cols, n_files))
   
   # Fill the array
   for (i in seq_len(n_files)) {
     for (r in seq_len(n_rows)) {
       mat_array[r, , i] <- all_matrices[[i]][[r]]
     }
   }
   
   # 4. Compute mean across 3rd dimension (files/seeds)
   T2_matrix <- apply(mat_array, c(1, 2), mean, na.rm = TRUE)
   
   # 5. Convert to "power"
   power_matrix <- 1 - T2_matrix
   

   
   # Rows = a values, columns = designs
   df_long <- melt(power_matrix)
   colnames(df_long) <- c("a_idx", "Design", "Power")
   
   # Add actual 'a' values
   a_seq <- seq(50, 100, length.out = nrow(power_matrix))  # adjust if needed
   df_long$a <- a_seq[df_long$a_idx]
   
   # Convert Design to factor with proper labels
   design_labels <- c("Non.rand", "Wal25", "Flexible", "SPP2", "IBD", "Flex.SPP2")
   df_long$Design <- factor(df_long$Design, levels = 1:length(design_labels), labels = design_labels)
   
   
   # Filter so Wal25 is only plotted for the last 'a' value
   df_long <- df_long[!(df_long$Design == "Wal25" & df_long$a != tail(a_seq, 1)), ]
   
   
   # Add a shape column: diamond for last Wal25, circle for everything else
   df_long$Shape <- ifelse(df_long$Design == "Wal25", 18, 16)  # 18 = diamond, 16 = circle
   
   # Plot
   p_power<- ggplot(df_long, aes(x = a, y = Power, color = Design, shape = factor(Shape))) +
     geom_line(size = 1.2) +
     geom_point(size = 3) +
     scale_color_manual(values = c(
       "Non.rand" = "#F8766D",
       "Wal25" = "#A6761D",
       "Flexible" = "#B79F00",
       "SPP2" = "#00BA38",
       "IBD" = "#00B0F6",
       "Flex.SPP2" = "#E76BF3"
     )) +
     scale_shape_manual(values = c("16" = 16, "18" = 18), guide = "none") +
     scale_y_continuous(breaks = seq(0, 1, by = 0.02)) +
     labs(x = "Percentage of Therapists available",
          y = "Power",
          color = "Design") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   print(p_power)
   
   # save as ggplot file for later editing
   saveRDS(
     p_power,
     file = "C:/Users/fvsk236/Desktop/Results_Therapist_Paper/Figures/N=240/Power_Plots/Model_1_nT_12/Model_1_nT_12_aP_1"
   )

   
   # Load the power plot fresh
   p_power <- readRDS("C:/Users/fvsk236/Desktop/Results_Therapist_Paper/Figures/N=240/Power_Plots/Model_1_nT_12/Model_1_power_nT_12_aP_0")
   
   # Update x-label and y-ticks
   p_power <- p_power +
     labs(x = "Percentage of Therapists available") +
     scale_y_continuous(breaks = seq(0, 1, by = 0.02))
   
   # Display the plot
   print(p_power)
    
   
   