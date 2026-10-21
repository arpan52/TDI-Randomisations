library(ggplot2)

# ==========================================================================
# ==========================================================================
design_map <- c(
  "Wal25"     = "E1",
  "SPP2"      = "E2",
  "Non.rand"  = "E3",
  "IBD"       = "N1",
  "Flexible"  = "N2",
  "Flex.SPP2" = "N3"
)

new_levels <- c("E1", "E2", "E3", "N1", "N2", "N3")   # legend order

relabel_designs <- function(df) {
  df$Design <- factor(
    unname(design_map[as.character(df$Design)]),
    levels = new_levels
  )
  df
}

# Same colours as before, attached to the new names
design_colors <- c(
  "E1" = "#A6761D",  # was Wal25
  "E2" = "#00BA38",  # was SPP2
  "E3" = "#F8766D",  # was Non.rand
  "N1" = "#00B0F6",  # was IBD
  "N2" = "#B79F00",  # was Flexible
  "N3" = "#E76BF3"   # was Flex.SPP2
)

old_levels <- c("Non.rand", "Wal25", "Flexible", "SPP2", "IBD", "Flex.SPP2")


# ==========================================================================
# Generic reader: works for power, rmse and time files
# (kept in the original structure; 'prefix' selects the file type)
# ==========================================================================
combine_results <- function(
    path,
    N_value,
    prefix,                 # "power_", "rmse_" or "time_to_target_"
    value_name,             # column name in the output
    transform = identity,   # e.g. function(x) 1 - x for power
    a_min = 50,
    a_max = 100,
    design_labels = old_levels
) {
  
  # 1. Locate files
  all_files      <- file.path(path, paste0(prefix, 1001:1405, ".rds"))
  existing_files <- all_files[file.exists(all_files)]
  
  if (length(existing_files) == 0)
    stop("No '", prefix, "' files found in: ", path)
  
  if (length(existing_files) < length(all_files))
    warning(sprintf("%s: only %d of %d files found in %s",
                    prefix, length(existing_files), length(all_files), path))
  
  message(sprintf("%s [N=%d]: %d files read", prefix, N_value, length(existing_files)))
  
  # 2. Read files
  all_matrices <- lapply(existing_files, readRDS)
  
  n_rows  <- length(all_matrices[[1]])
  n_cols  <- length(all_matrices[[1]][[1]])
  n_files <- length(all_matrices)
  
  # 3. Stack into array
  mat_array <- array(NA, dim = c(n_rows, n_cols, n_files))
  
  for (i in seq_len(n_files)) {
    for (r in seq_len(n_rows)) {
      mat_array[r, , i] <- all_matrices[[i]][[r]]
    }
  }
  
  # 4. Average over replicates (and transform, e.g. Type-II error -> power)
  mean_matrix <- transform(apply(mat_array, c(1, 2), mean, na.rm = TRUE))
  
  # 5. Long format
  df <- reshape2::melt(mean_matrix)
  colnames(df) <- c("a_idx", "Design", value_name)
  
  a_seq <- seq(a_min, a_max, length.out = nrow(mean_matrix))
  df$a  <- a_seq[df$a_idx]
  
  df$Design <- factor(
    df$Design,
    levels = seq_along(design_labels),
    labels = design_labels
  )
  
  df$N <- N_value
  
  df
}


# Shared post-processing: keep N=300 only for Flexible / Flex.SPP2,
# drop Wal25 at N=240 except at a_max. Uses OLD names on purpose.
assemble_all <- function(df_240, df_300) {
  df_300 <- df_300[df_300$Design %in% c("Flexible", "Flex.SPP2"), ]
  df_all <- rbind(df_240, df_300)
  
  a_max <- max(df_all$a)
  
  df_all <- df_all[!(
    df_all$Design == "Wal25" &
      df_all$N == 240 &
      df_all$a != a_max
  ), ]
  
  stopifnot(
    all(unique(df_240$a) == unique(df_300$a)),
    all(levels(df_all$Design) == old_levels)
  )
  
  print(table(df_all$Design, df_all$N))
  
  df_all
}


# ==========================================================================
# POWER PLOT
# ==========================================================================
df_240 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC_aP_1_240/",
  N_value    = 240,
  prefix     = "power_",
  value_name = "Power",
  transform  = function(x) 1 - x          # Type-II error -> power
)

df_300 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC_aP_1_300/",
  N_value    = 300,
  prefix     = "power_",
  value_name = "Power",
  transform  = function(x) 1 - x
)

df_all <- assemble_all(df_240, df_300)
df_all <- relabel_designs(df_all)          # rename AFTER filtering

y_min <- min(df_all$Power, na.rm = TRUE)
y_max <- max(df_all$Power, na.rm = TRUE)

pad   <- 0.01
y_min <- max(0, y_min - pad)
y_max <- min(1, y_max + pad)

p_power <- ggplot(
  df_all,
  aes(
    x = a,
    y = Power,
    color = Design,
    linetype = factor(N),
    group = interaction(Design, N)
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  
  scale_color_manual(
    values = design_colors,
    breaks = new_levels
  ) +
  
  scale_linetype_manual(
    values = c("240" = "solid", "300" = "dashed"),
    name = "Sample size (N)"
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(
      from = floor(y_min / 0.02) * 0.02,
      to   = ceiling(y_max / 0.02) * 0.02,
      by   = 0.02
    )
  ) +
  
  labs(
    x = "Percentage of Therapists available",
    y = "Power",
    color = "Strategy"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_power)

# save as ggplot file for later editing
saveRDS(
  p_power,
  file = "C:/Users/fvsk236/Desktop/Results_Paper_13/Figures/N=240/Power_and_Recruitment_Plots/aP_1/Power_plot_aP_1_240_300"
)

# ==========================================================================
# TIME-TO-TARGET PLOT
# ==========================================================================
df_time_240 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC_aP_1_240/",
  N_value    = 240,
  prefix     = "time_to_target_",
  value_name = "time"
)

df_time_300 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC_aP_1_300/",
  N_value    = 300,
  prefix     = "time_to_target_",
  value_name = "time"
)

df_time_all <- assemble_all(df_time_240, df_time_300)
df_time_all <- relabel_designs(df_time_all)

y_min <- min(df_time_all$time, na.rm = TRUE)
y_max <- max(df_time_all$time, na.rm = TRUE)

pad   <- 0.1
y_min <- y_min - pad
y_max <- y_max + pad

p_time <- ggplot(
  df_time_all,
  aes(
    x = a,
    y = time,
    color = Design,
    linetype = factor(N),
    group = interaction(Design, N)
  )
) +
  geom_line(linewidth = 1.2) +
  
  # --- points for all EXCEPT E3 (old Non.rand) ---
  geom_point(
    data = subset(df_time_all, Design != "E3"),
    size = 2.4,
    shape = 16
  ) +
  
  # --- points for E3 ON TOP (diamond) ---
  geom_point(
    data = subset(df_time_all, Design == "E3"),
    size = 3.2,
    shape = 18
  ) +
  
  scale_color_manual(
    values = design_colors,
    breaks = new_levels
  ) +
  
  scale_linetype_manual(
    values = c("240" = "solid", "300" = "dashed"),
    name = "Sample size (N)"
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(floor(y_min), ceiling(y_max), by = 1)
  ) +
  
  labs(
    x = "Percentage of Therapists available",
    y = "Randomisation time to target",
    color = "Strategy"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_time)

# save as ggplot file for later editing
saveRDS(
  p_time,
  file = "C:/Users/fvsk236/Desktop/Results_Paper_13/Figures/N=240/Power_and_Recruitment_Plots/aP_1/Recruitment_plot_aP_1_240_300"
)




# ==========================================================================
# RMSE PLOT
# NOTE: reads from HPC2_* folders (different from power/time). Confirm intended.
# ==========================================================================
df_rmse_240 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC2_aP_1_240/",
  N_value    = 240,
  prefix     = "rmse_",
  value_name = "rmse"
)

df_rmse_300 <- combine_results(
  path       = "C:/Users/fvsk236/Desktop/Results_Paper_13/HPC2_aP_1_300/",
  N_value    = 300,
  prefix     = "rmse_",
  value_name = "rmse"
)

df_rmse_all <- assemble_all(df_rmse_240, df_rmse_300)
df_rmse_all <- relabel_designs(df_rmse_all)

y_min <- min(df_rmse_all$rmse, na.rm = TRUE)
y_max <- max(df_rmse_all$rmse, na.rm = TRUE)

pad   <- 0.1
y_min <- y_min - pad
y_max <- y_max + pad

p_rmse <- ggplot(
  df_rmse_all,
  aes(
    x = a,
    y = rmse,
    color = Design,
    linetype = factor(N),
    group = interaction(Design, N)
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  
  scale_color_manual(
    values = design_colors,
    breaks = new_levels
  ) +
  
  scale_linetype_manual(
    values = c("240" = "solid", "300" = "dashed"),
    name = "Sample size (N)"
  ) +
  
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(
      from = floor(y_min / 0.02) * 0.02,
      to   = ceiling(y_max / 0.02) * 0.02,
      by   = 0.02
    )
  ) +
  
  labs(
    x = "Percentage of Therapists available",
    y = "rmse",
    color = "Strategy"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_rmse)

# save as ggplot file for later editing (was missing before; adjust filename if needed)
saveRDS(
  p_rmse,
  file = "C:/Users/fvsk236/Desktop/Results_Paper_13/Figures/N=240/Power_and_Recruitment_Plots/aP_1/RMSE_plot_aP_1_240_300"
)

