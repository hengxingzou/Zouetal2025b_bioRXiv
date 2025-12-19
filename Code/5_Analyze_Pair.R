library(tidyverse)
library(sf)
library(patchwork)
library(lmerTest)
library(gstat)
library(spdep)
library(spatialreg)

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

# Read bioclimatic variables by grid

Bioclim = read_csv("Data/Bioclim.csv") %>% 
  select(bio1, bio12, year, ST_12)

# Read human-induced land use changes

Anthromes = read_csv("Data/Anthromes_Grid_Years_Prop.csv") %>% 
  pivot_wider(names_from = Anthrome, values_from = Coverage) %>% 
  rename(year = Year) %>% 
  # bin data into larger land use categories
  mutate(Settlements = `11` + `12` + `22` + `23` + `24`, 
         Agriculture = `31` + `32` + `33` + `34` + `41` + `42` + `43`, 
         Cultured = `51` + `52` + `53` + `54`, 
         Wild = `61` + `62` + `63`) %>% 
  select(-(7:25), -LONGD, -LATD, -AREA_SQ_KM, -Grid_ID)


########## All Pairs, Maps ##########


# Tally by temporal trends and recorded status

tally_all_trends = Interaction_Pairs %>% 
  select(-(Relative_Abundance_Predator:year), -year_since_1969) %>% 
  distinct() %>% 
  mutate(Slope_IP_Sign = if_else(beta_logIP > 0, "pos", "neg")) %>% 
  group_by(grid, Recorded) %>% 
  summarize(Num_Pos = sum(Slope_IP_Sign == "pos"), Num_Neg = sum(Slope_IP_Sign == "neg"), .groups = "drop") %>% 
  pivot_longer(3:4, names_to = "Sign", values_to = "Count") %>% 
  full_join(latlong_shp, by = c("grid" = "ST_12"))

# Filter for only significant temporal trends

trends_lm = Interaction_Pairs %>% 
  select(Predator, Prey, Recorded, IP, year_since_1969, grid) %>% 
  group_by(Predator, Prey, grid, Recorded) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(log(IP) ~ year_since_1969, data = .)), 
         tidy_model = map(model, broom::tidy)) %>% 
  select(Predator, Prey, grid, Recorded, tidy_model) %>% 
  unnest(tidy_model) %>% 
  filter(term == "year_since_1969")

# Correct for the multiple comparisons (Benjamin-Hochberg)

trends_signif = trends_lm %>% 
  ungroup() %>%
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adjusted < 0.05) %>% 
  mutate(Slope_IP_Sign = if_else(estimate > 0, "pos", "neg"))

# Write data for future analysis
# Recommended because the previous step takes a long time

write_csv(trends_signif, "Data/Interaction_Pairs_corr_pval.csv")

# Read from saved data

trends_signif = read_csv("Data/Interaction_Pairs_corr_pval.csv")

# Filter for significant trends only

trends_signif_list = trends_signif %>% 
  unite(ID, c(Predator, Prey, grid), sep = "/", remove = F)

tally_all_signif = trends_signif %>% 
  group_by(grid, Recorded) %>% 
  summarize(Num_Pos = sum(Slope_IP_Sign == "pos"), Num_Neg = sum(Slope_IP_Sign == "neg"), .groups = "drop") %>% 
  pivot_longer(3:4, names_to = "Sign", values_to = "Count") %>% 
  full_join(latlong_shp, by = c("grid" = "ST_12"))

tally_recorded = tally_all_signif %>% 
  filter(Recorded == 1)

tally_unrecorded = tally_all_signif %>% 
  filter(Recorded == 0)

# Simplified visualization (main text)
# Tally net number of increasing trends in each grid for all recorded and predicted interactions

tally_all_sig_simp = tally_all_signif %>% 
  mutate(Count = if_else(Sign == "Num_Neg", -Count, Count)) %>% 
  group_by(grid) %>% 
  mutate(Net_Count = sum(Count)) %>% 
  select(-Recorded, -Sign, -Count) %>% 
  distinct()

pl_map_netcount = ggplot() + 
  geom_sf(data = latlong_shp, aes(geometry = geometry),
          color = "gray", size = 0.1,
          inherit.aes = F) +
  geom_sf(data = tally_all_sig_simp, 
          aes(fill = Net_Count, geometry = geometry), 
          color = "gray", size = 0.1, 
          inherit.aes = F) + 
  coord_sf() + 
  scale_fill_gradient2(name = "Net count of increasing interactions per grid", 
                       low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "bottom")

# Save figure; this is Figure 1A
# There is code below to save the entire Figure 1

ggsave("Figures/Maps_Net_Count.png", pl_map_netcount, width = 2000, height = 1500, unit = "px")

# Detailed visualization (supp)

overall_min = 0
overall_max = sqrt(max(tally_unrecorded$Count, tally_recorded$Count, 
                       na.rm = T))

pl_map_1 = ggplot() + 
  geom_sf(data = latlong_shp, aes(geometry = geometry),
          color = "gray", size = 0.1,
          inherit.aes = F) +
  geom_sf(data = tally_recorded, 
          aes(fill = sqrt(Count), geometry = geometry), 
          color = "gray", size = 0.1, 
          inherit.aes = F) + 
  coord_sf() + 
  scale_fill_gradientn(limits = c(overall_min, overall_max), 
                       colors = viridis::viridis(10), 
                       name = "Count", 
                       breaks = c(overall_min, overall_max),
                       labels = c(overall_min^2, overall_max^2)) +
  facet_grid(Sign ~ .) + 
  ggtitle("Temporal trends of recorded interactions") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "bottom")

pl_map_0 = ggplot() + 
  geom_sf(data = latlong_shp, aes(geometry = geometry),
          color = "gray", size = 0.1,
          inherit.aes = F) +
  geom_sf(data = tally_unrecorded, 
          aes(fill = sqrt(Count), geometry = geometry), 
          color = "gray", size = 0.1, 
          inherit.aes = F) + 
  coord_sf() + 
  scale_fill_gradientn(limits = c(overall_min, overall_max), 
                       colors = viridis::viridis(10), 
                       name = "Count", 
                       breaks = c(overall_min, overall_max),
                       labels = c(overall_min^2, overall_max^2)) +
  facet_grid(Sign ~ .) + 
  ggtitle("Temporal trends of predicted interactions") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "bottom")

# Save Figure S1

pl_map_both = pl_map_1 + pl_map_0
pl_map_both

ggsave("Figures/Maps_Both_Signif.png", pl_map_both, width = 4000, height = 3000, unit = "px")


########## All Pairs, Decomposition of Contributions ##########


# Filter for interactions with significant temporal trends only

all_trends = Interaction_Pairs %>% 
  select(-(Relative_Abundance_Predator:year), -year_since_1969) %>% 
  distinct() %>% 
  unite(ID, c(Predator, Prey, grid), sep = "/", remove = F) %>%
  filter(ID %in% trends_signif_list$ID) %>%
  select(-ID) %>%
  group_by(Predator, Prey) %>% 
  summarize(mean_beta_logA1 = mean(beta_logA1), sd_beta_logA1 = sd(beta_logA1), 
            mean_beta_logA2 = mean(beta_logA2), sd_beta_logA2 = sd(beta_logA2), 
            mean_beta_logIP = mean(beta_logIP), sd_beta_logIP = sd(beta_logIP)) %>% 
  left_join(IIP %>% select(Predator, Prey, IIP), by = c("Predator", "Prey")) %>% 
  mutate(Recorded = if_else(IIP == 1, 1, 0)) %>% 
  filter(!is.na(sd_beta_logA1) & !is.na(sd_beta_logA2) & !is.na(sd_beta_logIP))

# Simplfied visualization: both recorded and predicted

xrange = range(c(all_trends$mean_beta_logA1 - all_trends$sd_beta_logA1, 
                 all_trends$mean_beta_logA1 + all_trends$sd_beta_logA1)) * 1.05
yrange = range(c(all_trends$mean_beta_logA2 - all_trends$sd_beta_logA2, 
                 all_trends$mean_beta_logA2 + all_trends$sd_beta_logA2)) * 1.05
cont = expand_grid(x = seq(xrange[1], xrange[2], by = (xrange[2]-xrange[1])/150), 
                   y = seq(yrange[1], yrange[2], by = (yrange[2]-yrange[1])/150)) %>% 
  mutate(z = x + y)
zrange = range(cont$z)

# Assemble plot for all interactions

p_all = ggplot() + 
  geom_tile(data = cont, aes(x = x, y = y, fill = z)) +
  geom_pointrange(data = all_trends, 
                  aes(x = mean_beta_logA1, xmin = mean_beta_logA1 - sd_beta_logA1, xmax = mean_beta_logA1 + sd_beta_logA1, 
                      y = mean_beta_logA2, shape = as.factor(Recorded)), 
                  alpha = 0.15, size = 0.75) + 
  geom_pointrange(data = all_trends, 
                  aes(y = mean_beta_logA2, ymin = mean_beta_logA2 - sd_beta_logA2, ymax = mean_beta_logA2 + sd_beta_logA2, 
                      x = mean_beta_logA1, shape = as.factor(Recorded)), 
                  alpha = 0.15, size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  geom_abline(slope = -1, intercept = 0, color = "gray25") +
  xlab("Trend of log relative abundance of predator") + 
  ylab("Trend of log relative abundance of prey") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Trend of log(IP)") +
  scale_alpha_manual(values = c(1, 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

build_p_all = ggplot_build(p_all)
xrange_p_all = build_p_all$layout$panel_scales_x[[1]]$range$range
yrange_p_all = build_p_all$layout$panel_scales_y[[1]]$range$range

p_A1 = all_trends %>% 
  ggplot(aes(x = mean_beta_logA1)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  xlim(xrange_p_all) + 
  theme_void()

p_A2 = all_trends %>% 
  ggplot(aes(x = mean_beta_logA2)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  coord_flip() + 
  xlim(yrange_p_all) +
  theme_void()

pl_pairs = wrap_plots(p_A1, plot_spacer(), p_all, p_A2, 
                      nrow = 2, widths = c(10, 1), heights = c(1, 10), 
                      guides = "collect")

pl_pairs

# Save Figure 1B
# There is code below to save the entire Figure 1

ggsave("Figures/AllPairwise_All.png", pl_pairs, width = 2000, height = 2000, unit = "px")

# Specify recorded or unrecorded

trends_data_1 = all_trends %>% filter(Recorded == 1)
trends_data_0 = all_trends %>% filter(Recorded == 0)

# Set up the contour plots for recorded and predicted interactions

xrange_1 = range(c(trends_data_1$mean_beta_logA1 - trends_data_1$sd_beta_logA1, 
                   trends_data_1$mean_beta_logA1 + trends_data_1$sd_beta_logA1)) * 1.05
yrange_1 = range(c(trends_data_1$mean_beta_logA2 - trends_data_1$sd_beta_logA2, 
                   trends_data_1$mean_beta_logA2 + trends_data_1$sd_beta_logA2)) * 1.05
cont_1 = expand_grid(x = seq(xrange_1[1], xrange_1[2], by = (xrange_1[2]-xrange_1[1])/150), 
                     y = seq(yrange_1[1], yrange_1[2], by = (yrange_1[2]-yrange_1[1])/150)) %>% 
  mutate(z = x + y)
zrange_1 = range(cont_1$z)

xrange_0 = range(c(trends_data_0$mean_beta_logA1 - trends_data_0$sd_beta_logA1, 
                   trends_data_0$mean_beta_logA1 + trends_data_0$sd_beta_logA1)) * 1.05
yrange_0 = range(c(trends_data_0$mean_beta_logA2 - trends_data_0$sd_beta_logA2, 
                   trends_data_0$mean_beta_logA2 + trends_data_0$sd_beta_logA2)) * 1.05
cont_0 = expand_grid(x = seq(xrange_0[1], xrange_0[2], by = (xrange_0[2]-xrange_0[1])/150), 
                     y = seq(yrange_0[1], yrange_0[2], by = (yrange_0[2]-yrange_0[1])/150)) %>% 
  mutate(z = x + y)
zrange_0 = range(cont_0$z)

# Assemble plot for recorded interactions

p_all_1 = ggplot() + 
  geom_tile(data = cont_1, aes(x = x, y = y, fill = z)) +
  geom_pointrange(data = trends_data_1, 
                  aes(x = mean_beta_logA1, xmin = mean_beta_logA1 - sd_beta_logA1, xmax = mean_beta_logA1 + sd_beta_logA1, 
                      y = mean_beta_logA2), 
                  alpha = 0.15, size = 0.75) + 
  geom_pointrange(data = trends_data_1, 
                  aes(y = mean_beta_logA2, ymin = mean_beta_logA2 - sd_beta_logA2, ymax = mean_beta_logA2 + sd_beta_logA2, 
                      x = mean_beta_logA1), 
                  alpha = 0.15, size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  geom_abline(slope = -1, intercept = 0, color = "gray25") +
  xlab("Trend of log relative abundance of predator") + 
  ylab("Trend of log relative abundance of prey") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Trend of log(IP)") +
  scale_alpha_manual(values = c(1, 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

build_p_all_1 = ggplot_build(p_all_1)
xrange_p_all_1 = build_p_all_1$layout$panel_scales_x[[1]]$range$range
yrange_p_all_1 = build_p_all_1$layout$panel_scales_y[[1]]$range$range

p_A1_1 = trends_data_1 %>% 
  ggplot(aes(x = mean_beta_logA1)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  ggtitle("Recorded interactions") + 
  xlim(xrange_p_all_1) + 
  theme_void()

p_A2_1 = trends_data_1 %>% 
  ggplot(aes(x = mean_beta_logA2)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  coord_flip() + 
  xlim(yrange_p_all_1) +
  theme_void()

pl_pairs_1 = wrap_plots(p_A1_1, plot_spacer(), p_all_1, p_A2_1, 
                        nrow = 2, widths = c(10, 1), heights = c(1, 10), 
                        guides = "collect")

pl_pairs_1

# Assemble plot for predicted interactions

p_all_0 = ggplot() + 
  geom_tile(data = cont_0, aes(x = x, y = y, fill = z)) +
  geom_pointrange(data = trends_data_0, 
                  aes(x = mean_beta_logA1, xmin = mean_beta_logA1 - sd_beta_logA1, xmax = mean_beta_logA1 + sd_beta_logA1, 
                      y = mean_beta_logA2), 
                  alpha = 0.15, size = 0.75) + 
  geom_pointrange(data = trends_data_0, 
                  aes(y = mean_beta_logA2, ymin = mean_beta_logA2 - sd_beta_logA2, ymax = mean_beta_logA2 + sd_beta_logA2, 
                      x = mean_beta_logA1), 
                  alpha = 0.15, size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  geom_abline(slope = -1, intercept = 0, color = "gray25") +
  xlab("Trend of log relative abundance of predator") + 
  ylab("Trend of log relative abundance of prey") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Trend of log(IP)") +
  scale_alpha_manual(values = c(1, 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

build_p_all_0 = ggplot_build(p_all_0)
xrange_p_all_0 = build_p_all_0$layout$panel_scales_x[[1]]$range$range
yrange_p_all_0 = build_p_all_0$layout$panel_scales_y[[1]]$range$range

p_A1_0 = trends_data_0 %>% 
  ggplot(aes(x = mean_beta_logA1)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  ggtitle("Predicted interactions") + 
  xlim(xrange_p_all_0) + 
  theme_void()

p_A2_0 = trends_data_0 %>% 
  ggplot(aes(x = mean_beta_logA2)) +
  geom_density(fill = "gray75") + 
  geom_vline(xintercept = 0) +
  coord_flip() + 
  xlim(yrange_p_all_0) +
  theme_void()

pl_pairs_0 = wrap_plots(p_A1_0, plot_spacer(), p_all_0, p_A2_0, 
                        nrow = 2, widths = c(10, 1), heights = c(1, 10), 
                        guides = "collect")

pl_pairs_0

# Save Figure S3

ggsave("Figures/AllPairwise_Recorded.png", pl_pairs_1, width = 2000, height = 2000, unit = "px")
ggsave("Figures/AllPairwise_Predicted.png", pl_pairs_0, width = 2000, height = 2000, unit = "px")


########## Decomposition by Species ##########


pl_pairs_bypred = list()
all_pred = unique(Interaction_Pairs$Predator)

for (i in 1:length(all_pred)) {
  
  # Filter for interactions with significant temporal trends only
  
  all_trends_s = Interaction_Pairs %>% 
    filter(Predator == all_pred[i]) %>% 
    select(-(Relative_Abundance_Predator:year), -year_since_1969) %>% 
    distinct() %>% 
    unite(ID, c(Predator, Prey, grid), sep = "/", remove = F) %>%
    filter(ID %in% trends_signif_list$ID) %>%
    select(-ID) %>%
    group_by(Predator, Prey) %>% 
    summarize(mean_beta_logA1 = mean(beta_logA1), sd_beta_logA1 = sd(beta_logA1), 
              mean_beta_logA2 = mean(beta_logA2), sd_beta_logA2 = sd(beta_logA2), 
              mean_beta_logIP = mean(beta_logIP), sd_beta_logIP = sd(beta_logIP)) %>% 
    left_join(IIP %>% select(Predator, Prey, IIP), by = c("Predator", "Prey")) %>% 
    mutate(Recorded = if_else(IIP == 1, 1, 0)) %>% 
    filter(!is.na(sd_beta_logA1) & !is.na(sd_beta_logA2) & !is.na(sd_beta_logIP))

  # Simplfied visualization: both recorded and predicted
  
  xrange = range(c(min(all_trends_s$mean_beta_logA1 - all_trends_s$sd_beta_logA1, 0), 
                   all_trends_s$mean_beta_logA1 + all_trends_s$sd_beta_logA1)) * 1.05
  yrange = range(c(all_trends_s$mean_beta_logA2 - all_trends_s$sd_beta_logA2, 
                   all_trends_s$mean_beta_logA2 + all_trends_s$sd_beta_logA2)) * 1.05
  cont = expand_grid(x = seq(xrange[1], xrange[2], by = (xrange[2]-xrange[1])/150), 
                     y = seq(yrange[1], yrange[2], by = (yrange[2]-yrange[1])/150)) %>% 
    mutate(z = x + y)
  zrange = range(cont$z)
  
  # Assemble plot for all interactions
  
  p_all_s = ggplot() + 
    geom_tile(data = cont, aes(x = x, y = y, fill = z)) +
    geom_pointrange(data = all_trends_s, 
                    aes(x = mean_beta_logA1, xmin = mean_beta_logA1 - sd_beta_logA1, xmax = mean_beta_logA1 + sd_beta_logA1, 
                        y = mean_beta_logA2, shape = as.factor(Recorded)), 
                    alpha = 0.15, size = 0.75) + 
    geom_pointrange(data = all_trends_s, 
                    aes(y = mean_beta_logA2, ymin = mean_beta_logA2 - sd_beta_logA2, ymax = mean_beta_logA2 + sd_beta_logA2, 
                        x = mean_beta_logA1, shape = as.factor(Recorded)), 
                    alpha = 0.15, size = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
    geom_abline(slope = -1, intercept = 0, color = "gray25") +
    ggh4x::facet_wrap2(~Predator, scale = "free") + 
    xlab("Trend of log relative abundance of predator") + 
    ylab("Trend of log relative abundance of prey") + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         name = "Trend of log(IP)") +
    scale_alpha_manual(values = c(1, 1)) +
    theme_bw() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 12), 
          strip.text = element_text(size = 12),
          legend.position = "none")
  
  build_p_all_s = ggplot_build(p_all_s)
  xrange_p_all_s = build_p_all_s$layout$panel_scales_x[[1]]$range$range
  yrange_p_all_s = build_p_all_s$layout$panel_scales_y[[1]]$range$range
  
  p_A1_s = all_trends_s %>% 
    ggplot(aes(x = mean_beta_logA1)) +
    geom_density(fill = "gray75") + 
    geom_vline(xintercept = 0) +
    xlim(xrange_p_all_s) + 
    theme_void()
  
  p_A2_s = all_trends_s %>% 
    ggplot(aes(x = mean_beta_logA2)) +
    geom_density(fill = "gray75") + 
    geom_vline(xintercept = 0) +
    coord_flip() + 
    xlim(yrange_p_all_s) +
    theme_void()
  
  pl_pairs_s = wrap_plots(p_A1_s, plot_spacer(), p_all_s, p_A2_s, 
                          nrow = 2, widths = c(10, 1), heights = c(1, 10), 
                          guides = "collect")
  
  pl_pairs_bypred[[i]] = pl_pairs_s
  
}

names(pl_pairs_bypred) = all_pred

pl_pairs_bypred

# Save multi-panel figure; this is Figure 2

pl_pairs_combined = wrap_plots(pl_pairs_bypred, axis_titles = "collect")

ggsave("Figures/Pairwise_byPred_SignifOnly.png", pl_pairs_combined, width = 5000, height = 5000, unit = "px")


########## Analyze Turnover with Land Use Change ##########


# Calculate changes in environmental factors, then scale the data

env_change = Anthromes %>% 
  left_join(Bioclim, by = c("ST_12", "year")) %>% 
  filter(year %in% c(1970, 2021)) %>% 
  group_by(ST_12) %>% 
  summarize(diff_temp = bio1[year == 2021] - bio1[year == 1970], 
            diff_prec = bio12[year == 2021] - bio12[year == 1970], 
            diff_Settlements = Settlements[year == 2021] - Settlements[year == 1970], 
            diff_Agriculture = Agriculture[year == 2021] - Agriculture[year == 1970], 
            diff_Cultured = Cultured[year == 2021] - Cultured[year == 1970], 
            diff_Wild = Wild[year == 2021] - Wild[year == 1970], 
            diff_nonWild = diff_Settlements + diff_Agriculture + diff_Cultured)

# Get interaction networks in the form of adjacency matrices for each grid,
# then calculate dissimilarity between the beginning and the end with Frobenius norm

all_grids = unique(Interaction_Pairs$grid)

all_dissims = foreach(g = all_grids, .combine = rbind, .packages = "tidyverse", 
                      .errorhandling = "remove") %dopar% {
  
  dissim_vec_rec = c()
  dissim_vec_pre = c()
  dissim_vec_all = c()
  
  end_year = max(Interaction_Pairs$year_since_1969)
  
  single_grid_rec_y = Interaction_Pairs %>% 
    filter(grid == g, Documented == 1, year_since_1969 == 1) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()
  
  single_grid_pre_y = Interaction_Pairs %>% 
    filter(grid == g, Documented == 0, year_since_1969 == 1) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()
  
  single_grid_all_y = Interaction_Pairs %>% 
    filter(grid == g, year_since_1969 == 1) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()
  
  single_grid_rec_yp1 = Interaction_Pairs %>% 
    filter(grid == g, Documented == 1, year_since_1969 == end_year) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()
  
  single_grid_pre_yp1 = Interaction_Pairs %>% 
    filter(grid == g, Documented == 0, year_since_1969 == end_year) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()
  
  single_grid_all_yp1 = Interaction_Pairs %>% 
    filter(grid == g, year_since_1969 == end_year) %>% 
    select(Predator, Prey, IP) %>% 
    pivot_wider(names_from = Prey, values_from = IP, values_fill = 0) %>% 
    column_to_rownames("Predator") %>% 
    as.matrix()

  output = tibble(dissim_rec = norm(single_grid_rec_yp1 - single_grid_rec_y, type = "F"), 
                  dissim_pre = norm(single_grid_pre_yp1 - single_grid_pre_y, type = "F"), 
                  dissim_all = norm(single_grid_all_yp1 - single_grid_all_y, type = "F"),
                  grid = g)
  
  output
  
}

# Save intermediate dissimilarity results; this is recommended because the calculation is slow

saveRDS(all_dissims, "Code_Submission/All_Dissims.rds")

# Read saved intermediate results

all_dissims = readRDS("Code_Submission/All_Dissims.rds")

# Attach environmental change data

all_data_dissim = all_dissims %>% 
  left_join(env_change, by = c("grid" = "ST_12"))

# Log-transform dissimilarity values
# There is one grid with zero dissimilarity for predicted interaction only
# But doesn't matter because we only look at recorded or all interactions
# Also scale all predictor variables

all_data_dissim_t = all_data_dissim %>% 
  mutate(across(c(dissim_all, dissim_rec), ~log(.))) %>% 
  mutate(across(starts_with("diff"), ~c(scale(.))))

# Visualization
# Change here if want to look at all interaction pairs (dissim_all) or just recorded pairs (dissim_rec)

p_envchange = all_data_dissim_t %>% 
  pivot_longer(starts_with("diff"), names_to = "land_use_type", values_to = "diff") %>% 
  filter(!land_use_type %in% c("diff_nonWild", "diff_temp", "diff_prec")) %>% 
  ggplot(aes(x = diff, y = dissim_rec)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  ggh4x::facet_wrap2(. ~ land_use_type, scale = "free", 
                     labeller = as_labeller(c("diff_Agriculture" = "Agriculture", 
                                              "diff_Cultured" = "Cultured", 
                                              "diff_Settlements" = "Settlements", 
                                              "diff_Wild" = "Wild", 
                                              "diff_temp" = "Temperature", 
                                              "diff_prec" = "Precipitation"))) + 
  xlab("Scaled difference in environmental factors") + 
  ylab("Turnover of interactions") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 12),
        legend.position = "none")

p_envchange

# Save Figure S5

ggsave("Figures/Land_Use_Change.png", p_envchange, width = 3000, height = 3000, unit = "px")

# Run linear models with all environmental factors
# Change here if want to look at all interaction pairs (dissim_all) or just recorded pairs (dissim_rec)

m_all_env = lm(dissim_all ~ diff_temp + diff_prec + diff_Settlements + diff_Wild + diff_Cultured + 
                 diff_temp:diff_prec, 
               data = all_data_dissim_t)

summary(m_all_env)
car::ncvTest(m_all_env)
confint(m_all_env)

# Check for spatial autocorrelation in residuals

m_env_resid = tibble(Region = all_data_dissim_t$grid, resid = residuals(m_all_env)) %>% 
  separate_wider_delim(Region, names = c("lat", "long"), delim = "_") %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

env_vario = variogram(resid ~ 1, data = m_env_resid)

p_vario = env_vario %>% 
  ggplot(aes(x = dist, y = gamma)) + 
  geom_point() + 
  xlab("Distance") + 
  ylab("Semivariance") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_vario

# Turns out the residuals are strongly spatially correlated
# Save Figure S6

ggsave("Figures/Variogram.png", p_vario, width = 3000, height = 1500, unit = "px")

# Get neighbors

coords = all_data_dissim_sp[, c("long", "lat")]
nb = knearneigh(coords, k = 24)
nb = knn2nb(nb)
listw = nb2listw(nb, style = "W")

# Fit models with spatial lag (values are affected by neighboring values)

m_env_spatial = lagsarlm(dissim_rec ~ diff_temp + diff_prec + diff_Settlements + diff_Cultured + diff_Wild + 
                           diff_temp:diff_prec, 
                         data = all_data_dissim_sp, listw = listw)

summary(m_env_spatial)
confint(m_env_spatial)

# Pseudo R2

cor(m_env_spatial$y, fitted(m_env_spatial))^2

# Visualize model fitting results; not in the manuscript

m_all_env_params = broom::tidy(m_env_spatial)

m_all_env_confint = confint(m_env_spatial) %>% 
  as.data.frame() %>% 
  rownames_to_column("term")

labels = c("diff_temp" = "Change in temp", 
           "diff_prec" = "Change in prec", 
           "diff_temp:diff_prec" = "Change in temp:Change in prec", 
           "diff_Cultured" = "Change in cultured", 
           "diff_Settlements" = "Change in settlements", 
           "diff_Wild" = "Change in wild")

p_env_params = m_all_env_confint %>% 
  filter(term != "(Intercept)") %>% 
  filter(term != "diff_temp:diff_prec") %>%
  filter(term != "rho") %>% 
  mutate(signif = if_else(p.value >= 0.05, "no", "yes"), 
         term = recode(term, !!!labels)) %>% 
  mutate(term = factor(term, labels)) %>%
  ggplot(aes(x = term, y = estimate)) + 
  # geom_pointrange(aes(ymin = estimate - std.error, ymax = estimate + std.error, 
  #                     alpha = signif)) + 
  geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`,
                      alpha = signif)) +
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

# Save environmental association with turnover; this is currently not in the manuscript

ggsave("Figures/Env_Turnover.png", p_env_params, width = 3000, height = 1500, unit = "px")

# Save stats of the linear model

write_csv(m_all_env_params, "Data/Stats/turnover_envchange_lm.csv")
