library(tidyverse)
library(lsa)

source("Code/calculate_similarity.R")


########## Preparation and Empirical Mass Ratios ##########


# Read BBS species list with credible relative abundance estimates (470 species)

all_spp_list = read_csv("Filtered_Spp_List.csv") %>% 
  select(-aou) %>% 
  distinct()

# Read interaction data

all_interactions = read_csv("Data/Interactions/DietDB_Final.csv")

# Filter predators with more than ten records, then filter for only species in the species list

filtered_interactions = all_interactions %>% 
  filter(Scientific_Name %in% all_spp_list$latin, Prey_Scientific_Name %in% all_spp_list$latin) %>% 
  group_by(Scientific_Name) %>% 
  filter(n_distinct(Prey_Scientific_Name) > 10)

# Load functional traits

All_Funct_Data = read_csv("Data/All_Funct_Data.csv")

# Specify set of traits to use; Set 1-AVONET is used for our analyses
# Set 1-AVONET: AVONET base traits + migration, pelagic, nocturnal, parasite, nest
# See below for more information on other sets

trait_set = "1-AVONET"

# Specify threshold for selecting similar species

threshold = 0.95

# Specify strict trait matching

trait_matching = T

# File name based on the above parameters

file_name = paste(trait_set, threshold, trait_matching, sep = "_")

# Calculate empirical range of prey mass

empirical_ratios = all_interactions %>% 
  select(Scientific_Name, Prey_Scientific_Name) %>% 
  left_join(All_Funct_Data %>% 
              select(Species, Mass), by = c("Scientific_Name" = "Species")) %>% 
  rename(Predator_Mass = Mass) %>%
  left_join(All_Funct_Data %>% 
              select(Species, Mass), by = c("Prey_Scientific_Name" = "Species")) %>% 
  rename(Prey_Mass = Mass) %>% 
  group_by(Scientific_Name, Predator_Mass) %>% 
  summarize(min_mass = min(Prey_Mass), max_mass = max(Prey_Mass)) %>% 
  mutate(log_ratio = log(max_mass / min_mass))

# Visualize empirical ratios

p_emp_ratios = empirical_ratios %>% 
  filter(Scientific_Name %in% c("Falco spaverius", "Buteo swainsoni", "Buteo jamaicensis", 
                                "Accipiter cooperii", "Aquila chrysaetos", "Falco mexicanus", 
                                "Buteo regalis", "Buteo lineatus", "Haliaeetus leucocephalus", 
                                "Buteo platypterus", "Accipiter striatus", "Accipiter atricapillus", 
                                "Falco peregrinus", "Falco columbarius")) %>% 
  ggplot(aes(x = log(Predator_Mass), y = log_ratio)) + 
  geom_point(size = 3) + 
  ggrepel::geom_text_repel(aes(label = Scientific_Name), size = 5) + 
  xlab("log(Mass of predator)") + 
  ylab("log(Maximum mass of prey/Minimum mass of prey)") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12))

# Save Figure S4

ggsave("Figures/Mass_Ratios.png", p_emp_ratios, width = 3000, height = 2000, unit = "px")


########## Process Functional Traits ##########


# Construct trait matrix for all species that we have credible abundance estimates of

trait_mat = All_Funct_Data %>% 
  filter(Species %in% all_spp_list$latin) %>% 
  column_to_rownames(var = "Species") %>% 
  select(-starts_with("PC"))

# Spreading the AVONET categorical traits into single columns

cat_spread = trait_mat %>% 
  select(Habitat, Primary.Lifestyle) %>% 
  mutate(row_id = row_number(), value = 1) %>% 
  pivot_wider(names_from = Habitat, values_from = value, values_fill = list(value = 0), names_prefix = "AVONET_") %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = Primary.Lifestyle, values_from = value, values_fill = list(value = 0), names_prefix = "EltonTr_") %>% 
  select(-row_id)

# Combine the spread trait mat

trait_mat_spreaded = trait_mat %>% 
  select(-Habitat, -Primary.Lifestyle, -Diet.5Cat, -Trophic.Niche, -Trophic.Level) %>% 
  cbind(cat_spread)


########## Select Functional Traits ##########


# Categorize traits

morphology = colnames(trait_mat_spreaded)[1:11]
habitat_AVONET = colnames(trait_mat_spreaded)[71:80] # AVONET
habitat_IUCN = colnames(trait_mat_spreaded)[56:68] # IUCN
habitat_density = colnames(trait_mat_spreaded)[12]
lifestyle = colnames(trait_mat_spreaded)[81:85]
diet_perc = colnames(trait_mat_spreaded)[14:23]
foraging_perc = colnames(trait_mat_spreaded)[24:30]
nocturnal = colnames(trait_mat_spreaded)[32]
pelagic = colnames(trait_mat_spreaded)[31]
migration = colnames(trait_mat_spreaded)[13]
parasite = colnames(trait_mat_spreaded)[34]
nest_traits = colnames(trait_mat_spreaded)[35:53]
gen_length = colnames(trait_mat_spreaded)[54]
clutch_size = colnames(trait_mat_spreaded)[33]

# Choose from two sets of habitat traits

if (grepl("AVONET", trait_set)) {habitat = habitat_AVONET} else {habitat = habitat_IUCN}

# Rank traits by potential importance, add data types

trait_meta = tibble(category = c("morphology", "habitat", "habitat_density", "lifestyle", 
                                 "diet_perc", "foraging_perc", "nocturnal", "pelagic", "migration", 
                                 "parasite", "nest_traits", "gen_length", "clutch_size"), 
                    data_type = c("cont", "disc", "disc", "disc", 
                                  "perc", "perc", "disc", "disc", "disc", 
                                  "disc", "disc", "cont", "cont")) %>% 
  mutate(potential_importance = 1:nrow(.))

# Traits ranking 1-6 are most important ones: morphology, habitat, habitat density, lifestyle, diet, foraging

base_traits = trait_meta$category[1:6]

# Combine all trait data

trait_categories = list(morphology, habitat, habitat_density, lifestyle, 
                        diet_perc, foraging_perc, nocturnal, pelagic, migration, 
                        parasite, nest_traits, gen_length, clutch_size)

names(trait_categories) = trait_meta$category

# Discrete traits are in three different sets:

# Set 1-AVONET: AVONET base traits + migration, pelagic, nocturnal, parasite, nest
# Set 2-AVONET: AVONET base traits + pelagic, nocturnal
# Set 3-AVONET: AVONET base traits
# Set 1-IUCN: IUCN base traits + migration, pelagic, nocturnal, parasite, nest
# Set 2-IUCN: IUCN base traits + pelagic, nocturnal
# Set 3-IUCN: IUCN base traits

if (grepl("1", trait_set)) {combination = c(base_traits, "migration", "pelagic", "nocturnal", "parasite", "nest")}
if (grepl("2", trait_set)) {combination = c(base_traits, "pelagic", "nocturnal")}
if (grepl("3", trait_set)) {combination = base_traits}

selected_traits = trait_meta %>% 
  filter(category %in% combination)

cont = selected_traits %>% 
  filter(data_type == "cont") %>% 
  pull(category)
disc = selected_traits %>% 
  filter(data_type == "disc") %>% 
  pull(category)
perc = selected_traits %>% 
  filter(data_type == "perc") %>% 
  pull(category)

# Prepare trait data for calculation: scale continuous traits, transform percentage traits

cont_tr = as.data.frame(scale(trait_mat_spreaded[, unlist(trait_categories[names(trait_categories) %in% cont])]))
disc_tr = trait_mat_spreaded[, unlist(trait_categories[names(trait_categories) %in% disc])]
perc_tr = trait_mat_spreaded[, unlist(trait_categories[names(trait_categories) %in% perc])] / 100


########## Calculate Similarities ##########


# Get consumer and resource species from the database

consumer_sp = unique(filtered_interactions$Scientific_Name) # 14 species
resource_sp = all_spp_list$latin # all possible species; 470 in total

# Calculation based on cosine similarity

similarities_cosine = calculate_similarity(interactions = filtered_interactions, 
                                           consumer_sp = consumer_sp, resource_sp = resource_sp, 
                                           cont_traits = cont_tr, 
                                           disc_traits = disc_tr, 
                                           perc_traits = perc_tr, 
                                           top_n = 10, threshold = 0.95)

if (trait_matching) {
  
# Add extra filtering traits: body mass, habitat density

  similarities_mass = similarities_cosine %>% 
  left_join(trait_mat %>% rownames_to_column("Resource") %>% 
              select(Resource, Mass, Habitat.Density), by = "Resource") %>% 
  rename(Resource_Mass = Mass, Resource_Habitat_Density = Habitat.Density) %>% 
  left_join(trait_mat %>% rownames_to_column("Predicted_Resources") %>% 
              select(Predicted_Resources, Mass, Habitat.Density), by = "Predicted_Resources") %>% 
  rename(Predicted_Resource_Mass = Mass, Predicted_Resource_Habitat_Density = Habitat.Density) %>% 
  mutate(log_Mass_Ratio = log(Predicted_Resource_Mass / Resource_Mass))

# Filter by log mass ratio
# Harder constraint on the upper limit of prey mass than the lower limit 
# And habitat density: resource and predicted resource must has the same density
# These are extremely conservative filters

filtered_similarities = similarities_mass %>% 
  filter(log_Mass_Ratio <= 0.1 & log_Mass_Ratio >= -1) %>%
  filter(Resource_Habitat_Density == Predicted_Resource_Habitat_Density)
  
} else { filtered_similarities = similarities_cosine }

# Write data

write_csv(filtered_similarities, paste0("Data/Similarities/Similarities_", file_name, ".csv"))