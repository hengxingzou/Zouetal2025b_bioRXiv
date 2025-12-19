# Code modified from Chia Hsieh

calculate_similarity = function(interactions, # interaction database
                                consumer_sp, resource_sp, # species lists
                                cont_traits = NULL, disc_traits = NULL, perc_traits = NULL, # trait values by types of variables
                                top_n = 10, threshold = 0, similarity = "cosine", 
                                return_sim = F, verbose = T) {
  
  output = data.frame()
  
  for (c in consumer_sp) {
    
    if (verbose) {cat("Calculating potential preys for consumer", c, "\n")}
    
    subset_c = interactions %>% 
      filter(Scientific_Name == c)
    subset_r = unique(subset_c$Prey_Scientific_Name) # or whatever the column for the receiving species is
    
    # Normalize continuous variables by the value of the selected consumer
    
    cont_traits_norm = apply(as.matrix(cont_traits), 1, function(x) 
      as.matrix(cont_traits[c, ] - x))
    df_cont_traits_norm = data.frame(t(cont_traits_norm))
    colnames(df_cont_traits_norm) = colnames(cont_traits)
    
    # Combine discrete traits and percentage traits, then normalize

    other_traits = cbind(disc_traits, perc_traits)
    other_traits_norm = apply(as.matrix(other_traits), 1, function(x)
      as.matrix(other_traits[c, ]) - x)
    df_other_traits_norm = data.frame(t(other_traits_norm))
    colnames(df_other_traits_norm) = colnames(other_traits)
    df_other_traits_norm = df_other_traits_norm[rownames(df_cont_traits_norm), ]
    
    # Combine all traits, then remove the row of the consumer
    
    consumer_traits = cbind(df_cont_traits_norm, df_other_traits_norm)
    consumer_traits = consumer_traits[!(rownames(consumer_traits) == c), ]
    
    # Calculate similarities
    
    if (similarity == "gower") {
      
      similarity_mat = 1 - as.matrix(cluster::daisy(as.matrix(consumer_traits), metric = "gower", stand = F))
      
    } else {
      
      similarity_mat = cosine(t(as.matrix(consumer_traits)))
      
    }
    
    similarity_df = as.data.frame(similarity_mat)
    
    colnames(similarity_df) = rownames(consumer_traits)
    rownames(similarity_df) = rownames(consumer_traits)
    
    if (return_sim) { return(similarity_df) }
    
    # Rank and compare resources based on similarity
    
    for (r in subset_r) {
      
      similarity_r_df = similarity_df[order(similarity_df[, r], decreasing = T), ]
      
      if (!is.null(top_n)) {
        
        # Choose top n species with the highest similarity then filter by threshold
        # Top n species in addition to the resource species being looked at, so overall n+1 species
        
        similarity_r_topn = similarity_r_df[1:(top_n+1), ]
        similarity_r_topn = similarity_r_topn[similarity_r_topn[, r] >= threshold, ]
        
        # If no species was identified, then go to the next resource species
        
        if (nrow(similarity_r_topn) == 0) { next }
        
      }
      
      # See if the similar resource species has been documented in consumer-resource interactions
      
      documented = if_else(rownames(similarity_r_topn) %in% subset_r, 1, 0)
      
      # Save results
      
      output = rbind(output, data.frame("Consumer" = c, 
                                        "Resource" = r, 
                                        "Predicted_Resources" = rownames(similarity_r_topn), 
                                        "Similarity" = similarity_r_topn[, r], 
                                        "Documented" = documented, 
                                        check.names = F))
      
    }
    
  }
  
  return(output)
  
}
