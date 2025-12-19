library(tidyverse)
library(broom)

library(foreach)
library(doSNOW)

cl = makeCluster(8, outfile = "")
registerDoSNOW(cl)

# Read species abundance data
# Note--this is using relative abundance adjusted for detectability

All_Communities = read_csv("Data/Relative_Abundance.csv")

All_Communities = All_Communities %>% 
  # filter by year and geographic ranges
  filter(year >= 1970, lat <= 55, long >= -135)

# Read intrinsic interaction potential (IIP)

IIP = read_csv("Data/IP_IIP/Intrinsic_IP_1-AVONET_0.95_TRUE.csv")

# Calculate interaction potential (IP) for all species with abundance estimates

grids = unique(All_Communities$region)

interact_pairs = foreach(g = grids, .combine = rbind, .packages = "tidyverse") %dopar% {
  
  cat("Calculating interaction potential for all species pairs in", g, "\n")
  
  IP_g = data.frame()
  
  # Calculate interaction potential using relative abundances in each year
  
  for (y in unique(All_Communities$year)) {
    
    # Extract all co-occurring pairs in a grid
    
    local_comm_year = All_Communities %>% 
      filter(year == y, region == g)
    
    all_spp = unique(local_comm_year$latin)
    all_spp_pairs = expand_grid(all_spp, all_spp) %>% 
      rename(Spp_1 = `all_spp...1`, Spp_2 = `all_spp...2`) %>% 
      filter(Spp_1 != Spp_2)
    
    # Filter by pairs in the IIP dataset; 
    # order of predator and prey does not matter because expand_grid includes all unique combinations
    
    pairs_subset = all_spp_pairs %>% 
      inner_join(IIP, by = c("Spp_1" = "Predator", "Spp_2" = "Prey")) %>% 
      rename(Predator = Spp_1, Prey = Spp_2) %>% 
      mutate(Recorded = if_else(IIP == 1, 1, 0))
    
    # Calculate the interaction potential as 
    # intrinsic interaction potential * species 1 relative abundance * species 2 relative abundance
    # both relative abundances are scaled over time to remove temporal intercepts
    
    predator_abun = local_comm_year %>% 
      filter(latin %in% unique(pairs_subset$Predator)) %>% 
      select(latin, Relative_Abundance)
    
    prey_abun = local_comm_year %>% 
      filter(latin %in% unique(pairs_subset$Prey)) %>% 
      select(latin, Relative_Abundance)
    
    pairs_abun = pairs_subset %>% 
      left_join(predator_abun, by = c("Predator" = "latin"), relationship = "many-to-one") %>% 
      left_join(prey_abun, by = c("Prey" = "latin"), relationship = "many-to-one") %>% 
      rename(Relative_Abundance_Predator = Relative_Abundance.x, 
             Relative_Abundance_Prey = Relative_Abundance.y) %>% 
      mutate(IP = IIP * Relative_Abundance_Predator * Relative_Abundance_Prey, 
             year = y, grid = g)
    
    # Bind to the total data frame
    
    IP_g = rbind(IP_g, pairs_abun)
    
  }
  
  # Calculate temporal trends of IP for each pair and grid
  
  IP_g_temp = IP_g %>% 
    group_by(Predator, Prey) %>% 
    mutate(year_since_1969 = year - 1969, 
           beta_logA1 = cov(log(Relative_Abundance_Predator), year_since_1969)/var(year_since_1969), 
           beta_logA2 = cov(log(Relative_Abundance_Prey), year_since_1969)/var(year_since_1969), 
           beta_logIP = cov(log(IP), year_since_1969)/var(year_since_1969))
  
  IP_g_temp
  
}

# Write data

write_csv(interact_pairs, "Data/IP_IIP/adj_Pairs_IP_1-AVONET_0.95_TRUE.csv")
