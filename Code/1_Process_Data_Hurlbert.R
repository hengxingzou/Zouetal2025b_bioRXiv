library(tidyverse)
library(aviandietdb)


########## Load Data ##########


# Load Hurlbert et al. 2021 data
# See more info at https://aviandiet.unc.edu/
# and https://github.com/ahhurlbert/aviandietdb

data("dietdb")

# Filter for species level only

dietdb_filtered = dietdb %>% 
  # filter for avian preys
  filter(Prey_Class == "Aves") %>% 
  # filter for species-level observations only
  filter(Prey_Scientific_Name != "" & !grepl("/", Scientific_Name))

# Load BBS species list; see Methods on the derivation of this list

All_Spp_List = read_csv("Data/BBS_List_Num_Records.csv")

# Load functional traits

All_Funct_Data = read_csv("Data/All_Funct_Data.csv")


########## Fix Hurlbert et al. Taxonomy ##########


# Read the McTavish et al. taxonomy

McTavish = read_csv("OTT_crosswalk_2023.csv")

# Check predator names

setdiff(dietdb_filtered$Scientific_Name, McTavish$SCI_NAME)

# Found 1 difference (common name | Hurlbert et al. taxonomy | McTavish et al. taxonomy)

# Puerto Rican Owl | Megascops nudipes | Gymnasio nudipes

# Also name change

# American Goshawk | Accipiter gentilis | Accipiter atricapillus

# Fix predator names

dietdb_filtered_predcorr = dietdb_filtered %>% 
  mutate(Scientific_Name = if_else(Scientific_Name == "Megascops nudipes", "Gymnasio nudipes", Scientific_Name)) %>% 
  mutate(Scientific_Name = if_else(Scientific_Name == "Accipiter gentilis", "Accipiter atricapillus", Scientific_Name))

# Check prey names

setdiff(dietdb_filtered_predcorr$Prey_Scientific_Name, McTavish$SCI_NAME)

# Found 30 differences (common name | Hurlbert et al. taxonomy | McTavish et al. taxonomy)

# Neotropic Cormorant | Phalacrocorax brasilianus | Nannopterum brasilianum
# American Wigeon | Anas americana | Mareca americana
# Gadwall | Anas strepera | Mareca strepera
# Bluw-winged Teal | Anas discors | Spatula discors
# Northwestern Crow (American Crow) | Corvus caurinus | Corvus brachyrhynchos
# Snow Goose | Chen caerulescens | Anser caerulescens
# Pelagic Cormorant | Phalacrocorax pelagicus | Urile pelagicus
# Northern Shoveler | Anas clypeata | Spatula clypeata
# Double-crested Cormorant | Phalacrocorax auritus | Nannopterum auritum
# Hairy Woodpecker | Picoides villosus | Dryobates villosus
# Downy Woodpecker | Picoides pubescens | Dryobates pubescens
# Ladder-backed Woodpecker | Picoides scalaris | Dryobates scalaris
# Cinnamon Teal | Anas cyanoptera | Spatula cyanoptera
# Eurasian Jackdaw | Coloeus monedula | Corvus monedula
# Ruby-crowned Kinglet | Regulus calendula | Corthylio calendula
# Evening Grosbeak | Hesperiphona vespertina | Coccothraustes vespertinus
# Olive-throated Parakeet | Aratinga nana | Eupsittula nana
# Eastern Whip-poor-will | Caprimulgus vociferus | Antrostomus vociferus
# Band-winged Nightjar | Caprimulgus longirostris | Systellura longirostris
# Fire-eyed Diucon | Xolmis pyrope | Pyrope pyrope
# Long-tailed Meadowlark | Sturnella loyca | Leistes loyca
# Mourning Sierra Finch | Phrygilus fruticeti | Rhopospina fruticeti
# Band-tailed Sierra Finch | Phrygilus alaudinus | Rhopospina alaudina
# Spruce Grouse | Falcipennis canadensis | Canachites canadensis
# Fork-tailed Storm-Petrel | Oceanodroma furcata | Hydrobates furcatus
# Leach's Storm-Petrel | Oceanodroma leucorhoa | Hydrobates leucorhous
# Wedge-tailed Shearwater | Puffinus pacificus | Ardenna pacifica
# Emperor Goose | Chen canagica | Anser canagicus
# White-headed Woodpecker | Picoides albolarvatus | Dryobates albolarvatus
# Black-faced Grassquit | Tiaris bicolor | Melanospiza bicolor

Hurlbert_name_changes = tibble(
  Prey_Scientific_Name = c("Phalacrocorax brasilianus", 
                           "Anas americana", 
                           "Anas strepera", 
                           "Anas discors",
                           "Corvus caurinus", 
                           "Chen caerulescens",
                           "Phalacrocorax pelagicus", 
                           "Anas clypeata",
                           "Phalacrocorax auritus", 
                           "Picoides villosus", 
                           "Picoides pubescens", 
                           "Picoides scalaris", 
                           "Anas cyanoptera", 
                           "Coloeus monedula", 
                           "Regulus calendula", 
                           "Hesperiphona vespertina", 
                           "Aratinga nana", 
                           "Caprimulgus vociferus", 
                           "Caprimulgus longirostris", 
                           "Xolmis pyrope", 
                           "Sturnella loyca", 
                           "Phrygilus fruticeti", 
                           "Phrygilus alaudinus", 
                           "Falcipennis canadensis", 
                           "Oceanodroma furcata", 
                           "Oceanodroma leucorhoa", 
                           "Puffinus pacificus", 
                           "Chen canagica", 
                           "Picoides albolarvatus", 
                           "Tiaris bicolor"), 
  Prey_Scientific_Corrected = c("Nannopterum brasilianum",
                                "Mareca americana",
                                "Mareca strepera",
                                "Spatula discors",
                                "Corvus brachyrhynchos",
                                "Anser caerulescens",
                                "Urile pelagicus",
                                "Spatula clypeata",
                                "Nannopterum auritum",
                                "Dryobates villosus",
                                "Dryobates pubescens",
                                "Dryobates scalaris",
                                "Spatula cyanoptera",
                                "Corvus monedula", 
                                "Corthylio calendula",
                                "Coccothraustes vespertinus",
                                "Eupsittula nana",
                                "Antrostomus vociferus",
                                "Systellura longirostris",
                                "Pyrope pyrope",
                                "Leistes loyca",
                                "Rhopospina fruticeti",
                                "Rhopospina alaudina",
                                "Canachites canadensis",
                                "Hydrobates furcatus",
                                "Hydrobates leucorhous",
                                "Ardenna pacifica", 
                                "Anser canagicus",
                                "Dryobates albolarvatus", 
                                "Melanospiza bicolor"))

# Correct preys (all species needed to be corrected are preys)

dietdb_corrected = dietdb_filtered_predcorr %>% 
  filter(Prey_Scientific_Name %in% Hurlbert_name_changes$Prey_Scientific_Name) %>% 
  full_join(Hurlbert_name_changes, by = "Prey_Scientific_Name") %>% 
  select(-Prey_Scientific_Name) %>% 
  relocate(Prey_Scientific_Corrected, .after = "Prey_Genus") %>% 
  rename(Prey_Scientific_Name = Prey_Scientific_Corrected)

# Combine the full dataset

dietdb_corrected_full = dietdb_filtered_predcorr %>% 
  filter(Prey_Scientific_Name %in% McTavish$SCI_NAME) %>% 
  bind_rows(dietdb_corrected)

# Construct trait matrix for all species in the diet database

trait_mat = All_Funct_Data %>% 
  filter(Species %in% c(dietdb_corrected_full$Scientific_Name, 
                        dietdb_corrected_full$Prey_Scientific_Name)) %>% 
  column_to_rownames(var = "Species") %>% 
  select(-starts_with("PC")) 

# Total 366 species in the database (based on the BBS list)

# Filter the diet database with this list, also filter out cannibalism (focus only on interspecific interactions)

dietdb_final = dietdb_corrected_full %>% 
  filter(Scientific_Name %in% rownames(trait_mat) & Prey_Scientific_Name %in% rownames(trait_mat)) %>% 
  filter(Scientific_Name != Prey_Scientific_Name)

# Unload data

rm(McTavish, dietdb_filtered, dietdb_filtered_predcorr)

# Write data for future use

write_csv(dietdb_final, "Data/Interactions/DietDB_Final.csv")
