library(tidyverse)

library(foreach)
library(doSNOW)

cl = makeCluster(4)
registerDoSNOW(cl)

# Read data

similarities = read_csv("Data/Similarities/Similarities_1-AVONET_0.95_TRUE.csv")

# Calculate intrinsic interaction potential as the following:
# For each predator, pool all the preys, and for each prey with high similarities but
# undocumented interactions, calculate its average cosine similarities to all preys that have 
# documented interactions with the predator
# For documented predator-prey interactions the interaction potential is 1

# Interaction potential for the documented preys (all 1)

all_doc_potentials = similarities %>% 
  filter(Documented == 1) %>% 
  select(Consumer, Predicted_Resources, Documented) %>% 
  distinct() %>% 
  mutate(IIP = 1, IIP_SD = NA,
         IIP_n = 1) %>% 
  rename(Predator = Consumer, Prey = Predicted_Resources)

# Interaction potential for the undocumented (predicted) preys

predators_undoc = similarities %>% 
  filter(Documented == 0) %>% 
  pull(Consumer) %>% 
  unique()

all_undoc_potentials = foreach(pred = predators_undoc,
                               .combine = rbind, .packages = "tidyverse") %dopar% {
  
  interaction_subset = similarities %>% 
    filter(Consumer == pred)
  
  preys_doc = unique(interaction_subset$Resource)
  preys_undoc = setdiff(unique(interaction_subset$Predicted_Resources), preys_doc)
  
  potentials = data.frame()
  
  for (prey in preys_undoc) {
    
    undoc_similarities = interaction_subset %>% 
      filter(Predicted_Resources == prey) %>% 
      pull(Similarity)
    
    potentials = rbind(potentials, 
                       data.frame(Prey = prey, IIP = mean(undoc_similarities), 
                                  IIP_SD = sd(undoc_similarities), 
                                  IIP_n = length(undoc_similarities)))
    
  }
  
  potentials %>% mutate(Predator = pred, .before = Prey)
  
} %>% 
  mutate(Documented = 0)

# Combine both documented and undocumented preys

all_potentials = rbind(all_doc_potentials, all_undoc_potentials) %>% 
  arrange(Predator, IIP)

# Write data

write_csv(all_potentials, "Data/IP_IIP/Intrinsic_IP_1-AVONET_0.95_TRUE.csv")
