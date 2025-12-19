library(tidyverse)
library(sf)
library(patchwork)
library(lme4)

library(foreach)
library(doSNOW)

cl = makeCluster(4, outfile = "")
registerDoSNOW(cl)

# Read interaction potentials

Interaction_Pairs = read_csv("Data/IP_IIP/adj_Pairs_IP_1-AVONET_0.95_TRUE.csv") %>% 
  separate_wider_delim(grid, delim = "_", names = c("lat", "long"), cols_remove = F) %>% 
  mutate(across(c(lat, long), ~as.numeric(.))) %>% 
  filter(Recorded == 1)

# Read intrinsic interaction potentials

IIP = read_csv("Data/IP_IIP/Intrinsic_IP_1-AVONET_0.95_TRUE.csv")

# Read species groups

species_groups = readxl::read_xlsx("Data/Rosenberg2019_Science_SuppDataS1.xlsx", sheet = 1)

# Check taxonomy of species groups data

setdiff(unique(Interaction_Pairs$Prey), species_groups$sci_name)

# Found 6 name changes (Name in IP | Name in species groups data) and
# two missing species (Junco phaeonotus, Aphelocoma coerulescens)

species_groups_fixed = species_groups %>% 
  mutate(sci_name = if_else(species == "Double-crested Cormorant", "Nannopterum auritum", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Nashville Warbler", "Leiothlypis ruficapilla", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Ruby-crowned Kinglet", "Corthylio calendula", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Orange-crowned Warbler", "Leiothlypis celata", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Tennessee Warbler", "Leiothlypis peregrina", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Double_crested Cormorant", "Nannopterum auritum", sci_name)) %>% 
  mutate(sci_name = if_else(species == "Spruce Grouse", "Canachites canadensis", sci_name)) %>% 
  select(sci_name, Breeding.Biome, Winter.Biome, Family, bird.group, Migrate, AI, native) %>% 
  rename(aerial.insecti = AI)

# Impute the two missing species

species_groups_all = species_groups_fixed %>% 
  filter(sci_name %in% c("Junco hyemalis", "Aphelocoma californica")) %>% 
  mutate(sci_name = if_else(grepl("Junco", sci_name), "Junco phaeonotus", "Aphelocoma coerulescens")) %>% 
  mutate(Breeding.Biome = if_else(grepl("Junco", Breeding.Biome), "Western Forest", "Eastern Forest")) %>% 
  mutate(Winter.Biome = if_else(grepl("Junco", Winter.Biome), "Western Forest", "Eastern Forest")) %>% 
  mutate(Migrate = "R") %>% 
  bind_rows(species_groups_fixed)


########## Tally Interactions by Pair ##########


# Calculate predation risk

IP_prey = Interaction_Pairs %>% 
  group_by(year_since_1969, Prey, grid) %>% 
  mutate(sum_IP = sum(IP), n_pred = n_distinct(Predator)) %>% 
  group_by(Prey, year_since_1969, grid) %>% 
  mutate(sum_pred = sum(IIP * Relative_Abundance_Predator)) %>% 
  group_by(Prey, grid) %>% 
  left_join(species_groups_all, by = c("Prey" = "sci_name"))


########## Model Fitting by Groups ##########


# By bird groups

bird_groups = unique(species_groups_all$bird.group)

m_groups = foreach(b = bird_groups, .packages = "tidyverse") %dopar% {
  
  m_b = lmerTest::lmer(formula = log(sum_IP) ~ year_since_1969 + (1 | grid) + (1 | Prey), 
                       data = IP_prey %>% 
                         select(Prey, grid, sum_IP, year_since_1969, bird.group) %>% 
                         filter(bird.group == b) %>% 
                         distinct())
  
  param_groups = as.data.frame(summary(m_b)$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate(group = b)
  
  conf_int = as.data.frame(confint(m_b)) %>% 
    rownames_to_column("parameters")
  
  R2 = MuMIn::r.squaredGLMM(m_b)
  
  all_diag = conf_int %>% 
    mutate(marginal_r2 = R2[1], conditional_r2 = R2[2]) %>% 
    mutate(AIC = AIC(m_b)) %>% 
    mutate(group = b)
  
  list(param_groups, all_diag)
  
}

names(m_groups) = bird_groups

all_param_groups = do.call(rbind, lapply(m_groups, `[[`, 1)) %>% 
  mutate(p.adjusted = p.adjust(`Pr(>|t|)`, method = "BH"))
all_diag_groups = do.call(rbind, lapply(m_groups, `[[`, 2))

p_trend_group = all_param_groups %>% 
  filter(parameters == "year_since_1969") %>% 
  ggplot(aes(x = group, y = Estimate)) + 
  geom_pointrange(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`)) + 
  geom_hline(yintercept = 0) + 
  xlab("Bird group") + 
  ylab("Temporal trend in predation risk") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

p_trend_group

# By breeding biome

breeding_biomes = unique(species_groups_all$Breeding.Biome)

m_biomes = foreach(b = breeding_biomes, .packages = "tidyverse") %dopar% {
  
  m_b = lmerTest::lmer(formula = log(sum_IP) ~ year_since_1969 + (1 | grid) + (1 | Prey), 
                       data = IP_prey %>% 
                         select(Prey, grid, sum_IP, year_since_1969, Breeding.Biome) %>% 
                         filter(Breeding.Biome == b) %>% 
                         distinct())
  
  param_biomes = as.data.frame(summary(m_b)$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate(group = b)
  
  conf_int = as.data.frame(confint(m_b)) %>% 
    rownames_to_column("parameters")
  
  R2 = MuMIn::r.squaredGLMM(m_b)
  
  all_diag = conf_int %>% 
    mutate(marginal_r2 = R2[1], conditional_r2 = R2[2]) %>% 
    mutate(AIC = AIC(m_b)) %>% 
    mutate(group = b)
  
  list(param_biomes, all_diag)
  
}

names(m_biomes) = breeding_biomes

all_param_biomes = do.call(rbind, lapply(m_biomes, `[[`, 1)) %>% 
  mutate(p.adjusted = p.adjust(`Pr(>|t|)`, method = "BH"))
all_diag_biomes = do.call(rbind, lapply(m_biomes, `[[`, 2))

p_trend_biomes = all_param_biomes %>% 
  mutate(signif = if_else(p.adjusted >= 0.05, "no", "yes")) %>% 
  filter(parameters == "year_since_1969") %>% 
  ggplot(aes(x = group, y = Estimate, alpha = signif)) + 
  geom_pointrange(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`)) + 
  geom_hline(yintercept = 0) + 
  xlab("Breeding biome") + 
  ylab("Temporal trend in predation risk") +
  scale_alpha_manual(values = c(0.15, 1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        legend.position = "none")

p_trend_biomes

# Save Figure 3

ggsave("Figures/Prey_Predation_Risks.png", (p_trend_group | p_trend_biomes), 
       width = 3000, height = 1500, unit = "px")

# Save stats; files starting with Params_ are Table S2

write_csv(all_param_groups, "Data/Stats/Params_lmer_groups.csv")
write_csv(all_diag_groups, "Data/Stats/Diag_lmer_groups.csv")
write_csv(all_param_biomes, "Data/Stats/Params_lmer_biomes.csv")
write_csv(all_diag_biomes, "Data/Stats/Diag_lmer_biomes.csv")
