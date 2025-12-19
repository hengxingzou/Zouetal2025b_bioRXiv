library(tidyverse)
library(sf)
library(patchwork)
library(lme4)
library(rstan)

library(foreach)
library(doSNOW)

cl = makeCluster(4, outfile = "")
registerDoSNOW(cl)

# Read map

latlong_shp = read_sf("Data/Maps/BBS_LatLong_strata.shp")

# Read interaction potentials

Interaction_Pairs = read_csv("Data/IP_IIP/adj_Pairs_IP_1-AVONET_0.95_TRUE.csv") %>% 
  separate_wider_delim(grid, delim = "_", names = c("lat", "long"), cols_remove = F) %>% 
  mutate(across(c(lat, long), ~as.numeric(.)))

# Read intrinsic interaction potentials

IIP = read_csv("Data/IP_IIP/Intrinsic_IP_1-AVONET_0.95_TRUE.csv")

# Get all unique pairs

all_unique_pairs = Interaction_Pairs %>% 
  select(Predator, Prey) %>% 
  distinct()


########## Species-level Latitudinal Shifts ##########


# Calculation

All_Communities = read_csv("Data/Relative_Abundance.csv") %>% 
  filter(year >= 1970, lat <= 55, long >= -135) %>% 
  filter(latin %in% c(all_unique_pairs$Predator, all_unique_pairs$Prey))

all_spp = unique(All_Communities$latin)

species_shifts = foreach(spp = all_spp, .packages = "tidyverse", .combine = "rbind") %dopar% {
  
  cat("Fitting model for species", spp, "\n")
  
  spp_abun = All_Communities %>% 
    filter(latin == spp) %>% 
    mutate(year_since_1969 = year - 1969) %>% 
    select(latin, year_since_1969, region, lat, long, Relative_Abundance)
  
  spp_lat = spp_abun %>% 
    group_by(year_since_1969) %>% 
    reframe(weighted_lat = sum(lat*(Relative_Abundance/sum(Relative_Abundance))))
  
  m_lat = lm(weighted_lat ~ year_since_1969, data = spp_lat)
  coeffs = broom::tidy(m_lat) %>% 
    filter(term != "(Intercept)") %>% 
    mutate(latin = spp)
  
  coeffs
  
}

species_shifts_adj = species_shifts %>% 
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) %>% 
  mutate(signif = if_else(p.adjusted >= 0.05, F, T))

# Write data for other analyses

# write_csv(species_shifts_adj, "Data/Latitudinal_Shifts.csv")

# Analysis-predators

predator_shifts = species_shifts_adj %>% 
  filter(latin %in% all_unique_pairs$Predator) %>% 
  rename(Predator = latin) %>% 
  mutate(level = "full_range", 
         sign = if_else(estimate > 0, "pos", "neg"))

# Save stats; this is Table S2

write_csv(predator_shifts, "Data/Pred_Lat_Shifts.csv")

p_pred_shift = ggplot() + 
  geom_pointrange(data = predator_shifts, 
                  aes(x = reorder(Predator, -estimate), y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error, 
                      alpha = signif, color = sign), 
                  size = 0.5, linewidth = 1) + 
  geom_hline(yintercept = 0) + 
  xlab("Predator") + 
  ylab("Temporal trend of latitudinal centroid") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

p_pred_shift

# Tally numbers

predator_shifts %>% 
  filter(signif == T) %>% 
  pull(sign) %>% 
  table()

# Save figure; this is an intermediate figure not in the manuscript

ggsave("Figures/Predator_Lat_Shift.png", p_pred_shift, width = 2000, height = 2000, unit = "px")

# Analysis-preys

prey_shifts = species_shifts_adj %>% 
  filter(latin %in% all_unique_pairs$Prey) %>% 
  rename(Prey = latin) %>% 
  mutate(level = "full_range", 
         sign = if_else(estimate > 0, "pos", "neg"))

p_prey_shift = ggplot() + 
  geom_pointrange(data = prey_shifts, 
                  aes(x = reorder(Prey, -estimate), y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error, 
                      alpha = signif, color = sign), 
                  size = 0.5, linewidth = 1) + 
  geom_hline(yintercept = 0) + 
  xlab("Predator") + 
  ylab("Temporal trend of latitudinal centroid") + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

# Tally numbers

prey_shifts %>% 
  filter(signif == T) %>% 
  pull(sign) %>% 
  table()

# Save figure; this is an intermediate figure not in the manuscript

ggsave("Figures/Prey_Lat_Shift.png", p_prey_shift, width = 2000, height = 5000, unit = "px")

rm(All_Communities)


########## Weighted Mean Latitudes by Interactions ##########


# Calculate weighted average latitudes and their temporal trends for all pairs
# for predator and prey in their overlapping spatial ranges

data_all_pairs = foreach(n = 1:nrow(all_unique_pairs), .packages = "tidyverse") %dopar% {
  
  pred = all_unique_pairs[n, ]$Predator
  prey = all_unique_pairs[n, ]$Prey
  
  cat("Calculating for predator", pred, "and prey", prey, "\n")
  
  # Scale data
  
  single_pair_n = Interaction_Pairs %>% 
    filter(Predator == pred, Prey == prey) %>% 
    group_by(year_since_1969) %>% 
    reframe(weighted_pred_lat = sum(lat*(Relative_Abundance_Predator/sum(Relative_Abundance_Predator))), 
            weighted_prey_lat = sum(lat*(Relative_Abundance_Prey/sum(Relative_Abundance_Prey))), 
            weighted_IP_lat = sum(lat*(IP/sum(IP)))) %>% 
  mutate(Predator = pred, Prey = prey)
  
  m_pred = lm(weighted_pred_lat ~ year_since_1969, data = single_pair_n)
  m_prey = lm(weighted_prey_lat ~ year_since_1969, data = single_pair_n)
  m_IP = lm(weighted_IP_lat ~ year_since_1969, data = single_pair_n)
  
  # Extract parameters
  
  coeff_pred = broom::tidy(m_pred) %>% 
    mutate(group = "pred")
  coeff_prey = broom::tidy(m_prey) %>% 
    mutate(group = "prey")
  coeff_IP = broom::tidy(m_IP) %>% 
    mutate(group = "IP")
  
  coeffs = rbind(coeff_pred, coeff_prey, coeff_IP) %>% 
    filter(term != "(Intercept)") %>%
    mutate(Predator = pred, Prey = prey)
  
  list(single_pair_n, coeffs)
  
}

all_weighted_lat = do.call(rbind, lapply(data_all_pairs, function(x) x[[1]])) %>% 
  left_join(IIP[, c("Predator", "Prey", "Documented")], by = c("Predator", "Prey"))
coeff_all_pairs = do.call(rbind, lapply(data_all_pairs, function(x) x[[2]])) %>% 
  left_join(IIP[, c("Predator", "Prey", "Documented")], by = c("Predator", "Prey")) %>% 
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) %>% 
  mutate(signif = if_else(p.adjusted >= 0.05, F, T), 
         sign = if_else(estimate > 0, "pos", "neg"))

# Tally numbers

xx = coeff_all_pairs %>% 
  filter(group == "IP", signif == T)

table(xx$sign)
table(xx$sign, xx$Documented)
chisq.test(table(xx$sign, xx$Documented))


########## Analyze Pairwise Latitudinal Shifts ##########


# Predators, whole-range and ranges of interactions

coeff_predators = coeff_all_pairs %>% 
  filter(group == "pred") %>% 
  select(term, estimate, std.error, statistic, p.adjusted, Predator, signif, sign) %>% 
  mutate(level = "pairwise_overlap")

p_IP_shift_bypred = ggplot() + 
  geom_point(data = coeff_predators, 
             aes(x = Predator, y = estimate, alpha = signif, color = sign), 
             position = position_jitter(0.15)) + 
  geom_pointrange(data = predator_shifts, 
                  aes(x = reorder(Predator, -estimate), y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error, 
                      alpha = signif), 
                  color = "black", size = 0.5, linewidth = 1) + 
  geom_hline(yintercept = 0) + 
  ylab("Temporal trend of latitudinal centroid") + 
  scale_color_manual(values = c("blue", "red"))
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

p_IP_shift_bypred

# Save figure; this is Figure 4

ggsave("Figures/IP_Lat_Shift_byPred.png", p_IP_shift_bypred, width = 2000, height = 2000, unit = "px")


########## Decomposition of Contributions to Latitudinal Shifts ##########


# Decompose IP shifts by predator and prey shifts, for their spatial overlaps

pl_all = list()
all_pred = unique(all_unique_pairs$Predator)

for (i in 1:length(all_pred)) {
  
  # Species-level shifts of the predator
  
  p = predator_shifts %>% 
    filter(Predator == all_pred[i])
  p_shift = p$estimate
  p_se = p$std.error
  p_signif = p$signif
  
  # Pairwise latitudinal shifts for the overlapping ranges
  
  single_pred = coeff_all_pairs %>% 
    filter(Predator == all_pred[i]) %>% 
    filter(signif == T) %>%
    select(estimate, group, Predator, Prey, Documented) %>% 
    pivot_wider(names_from = group, values_from = estimate) %>% 
    filter(!is.na(pred) & !is.na(IP) & !is.na(prey)) %>% 
    rename(pred_shift = pred, prey_shift = prey, IP_shift = IP)
  
  p_single_main = single_pred %>% 
    ggplot(aes(x = pred_shift, y = prey_shift, color = IP_shift, shape = as.factor(Documented))) + 
    geom_point(size = 3) + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    geom_vline(xintercept = p_shift, color = "red") +
    geom_vline(xintercept = p_shift - p_se, color = "red", linetype = 2) +
    geom_vline(xintercept = p_shift + p_se, color = "red", linetype = 2) +
    scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) + 
    xlab("Predator latitudinal shifts") + 
    ylab("Prey latitudinal shifts") +
    ggh4x::facet_wrap2(.~Predator, scales = "free") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 12), 
          strip.text = element_text(size = 12),
          legend.position = "none")
  
  build_p_single_main = ggplot_build(p_single_main)
  xrange_p_single_main = build_p_single_main$layout$panel_scales_x[[1]]$range$range
  yrange_p_single_main = build_p_single_main$layout$panel_scales_y[[1]]$range$range
  
  if (p_signif) {
    p_single_main = p_single_main + 
      annotate(geom = "text", 
               x = xrange_p_single_main[2] - 0.1*(xrange_p_single_main[2] - xrange_p_single_main[1]), 
               y = yrange_p_single_main[2] - 0.1*(yrange_p_single_main[2] - yrange_p_single_main[1]), 
               label = "*", size = 20)
  }
  
  p_pred = single_pred %>% 
    ggplot(aes(x = pred_shift)) +
    geom_density(fill = "gray75") + 
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = p_shift, color = "red") +
    geom_vline(xintercept = p_shift - p_se, color = "red", linetype = 2) +
    geom_vline(xintercept = p_shift + p_se, color = "red", linetype = 2) +
    xlim(xrange_p_single_main) + 
    theme_void()
  
  p_prey = single_pred %>% 
    ggplot(aes(x = prey_shift)) +
    geom_density(fill = "gray75") + 
    geom_vline(xintercept = 0) +
    coord_flip() +
    xlim(yrange_p_single_main) + 
    theme_void()
  
  p_single_full = wrap_plots(p_pred, plot_spacer(), p_single_main, p_prey, 
                             nrow = 2, widths = c(10, 1), heights = c(1, 10), 
                             guides = "collect")
  
  pl_all[[i]] = p_single_full
  
}

names(pl_all) = all_pred

pl_all

# Save multi-panel figure; this is Figure 5

pl_combined = wrap_plots(pl_all, axis_titles = "collect")

ggsave("Figures/IP_LatShift_All_Overlap_SignifOnly.png", pl_combined, width = 5000, height = 5000, unit = "px")