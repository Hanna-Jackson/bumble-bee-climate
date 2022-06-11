## Written by: Laura Melissa Guzman

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(stringr)
library(purrr)
library(ggnewscale)

## ~~~~~~~~~~~~ Figure 2 - Heatplot Table of Predictors ~~~~~~~~~~~~

plot_heatplots <- function(res = 100, initial_year = NA){

  ## Here we load the file we need based on the resolution
  ## Add the correct file paths based on your own runs if you're doing that
  if(res == 100){
    
    #### load model with era only
    load("output/modeloutput-era-100.RData", verbose = TRUE)
    
    my.data.m.era <- my.data
    bugs.m.era <- bugs
    summary.m.era <- bugs.m.era$BUGSoutput$summary
    
    ### model with predictors no era
    
    load("output/modeloutput-env-100.RData", verbose = TRUE)
    
    ## set heatplot parameters
    
    # res 100
    limits_era   <- c(-0.5, 0.5)
    limits_tmax   <- c(-5.6, 2.4)
    limits_prec   <- c(-2.4,4.8)
    limits_floral   <- c(-1.4, 0.9)
    limits_perchange <- c(-0.86, 2.2)
    limits_IUCN <- c(-1, 0.9)
    
  } else if(res == 250){
    res <- 250
    
    #### load model with era only
    load("output/modeloutput-era-250.RData", verbose = TRUE)
    
    my.data.m.era <- my.data
    bugs.m.era <- bugs
    
    summary.m.era <- bugs.m.era$BUGSoutput$summary
    
    ### model with predictors no era
    
    load("output/modeloutput-env-250.RData", verbose = TRUE)
    
    ## set the heatplot parameters
    
    ## res 250
    limits_era <- c(-0.4, 0.4)
    limits_tmax <- c(-5.9, 4.3)
    limits_prec <- c(-2.3,3)
    limits_floral <- c(-2.1, 1.39)
    limits_perchange <- c(-0.65, 0.83)
    limits_IUCN <- c(-1, 0.9)
    
  }else if(res == 50){
    
    res <- 50
    
    #### load model with era only
    load("output/modeloutput-era-50.RData", verbose = TRUE)
    
    my.data.m.era <- my.data
    bugs.m.era <- bugs
    
    summary.m.era <- bugs.m.era$BUGSoutput$summary
    
    ### model with predictors no era
    
    load("output/modeloutput-env-50.RData", verbose = TRUE)
    
    ## set the heatplot parameters 
    
    # res 50

    limits_era <- c(-0.46, 0.56)
    limits_tmax <- c(-5, 2.6)
    limits_prec <- c(-1.9,5.6)
    limits_floral <- c(-1.2, 0.71)
    limits_perchange <- c(-0.88, 4.26)
    limits_IUCN <- c(-1, 0.9)
    
  }
  
  if(initial_year == 1960){
    
    res <- 1960
    
    # #### load model with era only
    load("output/modeloutput-start1960-era-100.RData", verbose = TRUE)
    
    my.data.m.era <- my.data
    bugs.m.era <- bugs
    
    summary.m.era <- bugs.m.era$BUGSoutput$summary
    
    ### model with predictors no era
    
    load("output/modeloutput-start1960-env-100.RData", verbose = TRUE)
    
    ## set heatplot parameters
    
    # 1960

    limits_era <- c(-1.26, 1.1)
    limits_tmax <- c(-4.95, 1.9)
    limits_prec <- c(-1.9,5.6)
    limits_floral <- c(-1.37, 1.39)
    limits_perchange <- c(-0.84, 1.4)
    limits_IUCN <- c(-1, 0.9)
    
  }
  
  
  ## assign other name for the predictor model
  
  my.data.m.pred <- my.data
  bugs.m.pred <- bugs
  
  summary.m.pred <- bugs.m.pred$BUGSoutput$summary
  
  ## extact the important information, i.e. psimeanmaxt/floral/meanprec
  
  cols <- c("mean", "2.5%", "97.5%")
  species_names <- rownames(my.data.m.pred$species.ranges)
  species_names_clean <- paste(1:length(species_names), str_remove(species_names, "species."))
  
  rows_temp <- rownames(summary.m.pred)[str_detect(rownames(summary.m.pred), "psi.meanmaxt\\[")]
  rows_floral <- rownames(summary.m.pred)[str_detect(rownames(summary.m.pred), "psi.floral\\[")]
  rows_prec <- rownames(summary.m.pred)[str_detect(rownames(summary.m.pred), "psi.meanprec\\[")]
  
  ## combine all of the useful info into a dataframe
  
  mpred_tmax <- data.frame(summary.m.pred[rows_temp, cols], species = species_names_clean, model = 'pred', term = 'tmax')
  
  mpred_prec <- data.frame(summary.m.pred[rows_prec, cols], species = species_names_clean, model = 'pred', term = 'prec')
  
  mpred_floral <- data.frame(summary.m.pred[rows_floral, cols], species = species_names_clean, model = 'pred', term = 'floral')
  
  
  ##### For era model
  
  ## clean the era summary
  
  cols <- c("mean", "2.5%", "97.5%")
  species_names <- rownames(my.data.m.era$species.ranges)
  
  species_names_order <- str_remove(species_names, "species.")
  
  species_names_clean <- paste(1:length(species_names), str_remove(species_names, "species."))
  
  rows_era <- rownames(summary.m.era)[str_detect(rownames(summary.m.era), "psi.era\\[")]
  m_era_era <- data.frame(summary.m.era[rows_era, cols], species = species_names_clean, model = 'era', term = 'era')
  
  all_summary <- bind_rows(mpred_tmax, mpred_prec, mpred_floral, m_era_era)
  
  rownames(all_summary) <- NULL
  
  ## add winners - losers to the data frame
  
  temp_wl <- data.frame(term = 'tmax', species_n = c(c(2, 6, 7, 30, 32, 34, 41, 42, 43), c(44, 46)), w_l = c(rep("L", 9), rep("W", 2)))
  prec_wl <- data.frame(term = 'prec', species_n = c(c(11)), w_l = c(rep("L", 0), rep("W", 1)))
  floral_wl <- data.frame(term = 'floral', species_n = c(c(28), c(3, 38)), w_l = c(rep("L", 1), rep("W", 2)))
  
  w_l_df <- data.frame(species_n = 1:46, species = paste(1:length(species_names), str_remove(species_names, "species."))) %>% 
    right_join(bind_rows(temp_wl, prec_wl, floral_wl))
  
  ## clean the compiled information
  
  pred_predictors_psi <- all_summary %>% 
    mutate(sig = ifelse(`X2.5.`< 0 & 0 < `X97.5.`, FALSE, TRUE)) %>% 
    mutate(term = factor(term, levels = c('era', 'tmax', 'prec', 'floral'))) %>% 
    group_by(term) %>% 
    mutate(col = ifelse(mean > quantile(mean, 0.8), 'white', 'black')) %>% 
    left_join(w_l_df) %>% 
    mutate(w_l = ifelse(is.na(w_l), "", w_l)) %>% 
    mutate(species = factor(species, levels = rev(species_names_clean)))
  
  ### create supplementary table with values
  
  supp_table <- pred_predictors_psi %>% 
    mutate(values = paste0(round(mean, 2), " (", round(`X2.5.`,2), ", ", round(`X97.5.`,2), ")")) %>% 
    select(values, species, term) %>%
    pivot_wider(names_from = "term", values_from = "values") %>% 
    select(species, era, tmax, prec, floral) %>% 
    mutate(all_cols = paste0(species, " & ", era, " & ", tmax, " & ", prec, " & ", floral, " \\ "))
  
  paste(supp_table$all_cols, collapse = "")
  
  ## add proportional change and IUCN trends
  
  ## expit and logit functions
  expit <- function(x) 1/(1+exp(-x))
  logit <- function(x) log(x/(1-x))
  
  get.mean.bci <- function(chains)
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  
  
  ## load bee status from IUCN
  
  bee_status <- read.csv("Bee status - Sheet6.csv")
  
  ## calculate proportional change with the era model
  
  summary.m.era <- bugs.m.era$BUGSoutput$summary
  
  nsite <- my.data.m.era$nsite
  nvisit <- my.data.m.era$nvisit
  nera <- my.data.m.era$nera
  nspecies <- my.data.m.era$nspecies
  
  summary.m.era['mu.psi.era',]
  
  ## clean the era summary
  
  cols <- c("mean", "2.5%", "97.5%")
  species_names <- rownames(my.data.m.era$species.ranges)
  
  rows_era <- rownames(summary.m.era)[str_detect(rownames(summary.m.era), "psi.era\\[")]
  all_summary <- data.frame(summary.m.era[rows_era, cols], species = species_names, model = 'era', term = 'era')
  
  rownames(all_summary) <- NULL
  
  
  ## Calculate the percent change in occupancy 
  
  sims.mat.m.era <- bugs.m.era$BUGSoutput$sims.matrix
  sims.arr.m.era <- bugs.m.era$BUGSoutput$sims.array
  nsite.m.era <- my.data.m.era$nsite
  nvisit.m.era <- my.data.m.era$nvisit
  nera.m.era <- my.data.m.era$nera
  nspecies.m.era <- my.data.m.era$nspecies
  
  calc.occ.diff.era <- function(species) {
    
    ## figure out which sites are in the species' range
    sites <- which(my.data.m.era$species.ranges[species,])
    
    ## calculate mean site area for species
    sitearea.mean <- mean(my.data.m.era$sitearea[sites])
    
    get.chains <- function(era) {
      
      expit(sims.mat.m.era[,'psi.0'] +
              sims.mat.m.era[,sprintf('psi.sp[%d]', species)] +
              sims.mat.m.era[,sprintf('psi.era[%d]', species)] * era + ## comment out for model without era
              sims.mat.m.era[,'psi.sitearea'] * sitearea.mean)
    }
    chains.era1 <- get.chains(era=min(my.data.m.era$era))
    chains.era2 <- get.chains(era=max(my.data.m.era$era))
    chains.diff <- (chains.era2-chains.era1)/chains.era1
    
    get.mean.bci(chains.diff)
  }
  
  all_actual <- sapply(1:nspecies, function(ss) calc.occ.diff.era(ss))
  
  per_change <- data.frame(t(all_actual), species = str_remove(species_names, 'species.')) %>% 
    dplyr::select(percent_change = mean, species)
  
  
  ### summarise info from all model outputs and IUCN trends
  
  all_change_occupancy <- per_change %>% 
    left_join(bee_status, by = c('species' = 'Species')) %>% 
    left_join(data.frame(species = paste(1:length(species_names), str_remove(species_names, "species.")), species_old = str_remove(species_names, "species.")), by = c("species" = 'species_old')) %>% 
    dplyr::select(`species.y`, percent_change, Mean.Change) %>% 
    mutate(Mean.Change = ifelse(str_detect(species.y, 'distinguendus'), NA, Mean.Change/100)) %>% 
    rename(species = species.y) %>% 
    pivot_longer(names_to = 'term', values_to = 'mean', -c(species)) %>% 
    mutate(species = factor(species, levels = rev(species_names_clean))) %>% 
    mutate(sig = TRUE) %>% 
    group_by(term) %>% 
    mutate(col = ifelse(mean > quantile(mean, 0.8, na.rm = TRUE), 'white', 'black'))
  
  all_model_outputs <- bind_rows(pred_predictors_psi, all_change_occupancy)
  
  
  ## do the heatplot
  
  # species names italizised and change the order of the legends and the columns
  # change the colour scaling to be maxed
  # precipitation and 
  
  mypalette<-brewer.pal(9,"YlGnBu")
  
  heatplot_pred_partial <- ggplot() +
    geom_tile(data = filter(all_model_outputs, term == "era"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2(expression(psi["era"]), limits = limits_era,
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 1))  +
    new_scale("fill") +
    geom_tile(data = filter(all_model_outputs, term == "percent_change"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2("Prop.\nChange", limits = limits_perchange, 
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 2)) +
    new_scale("fill") +
    geom_tile(data = filter(all_model_outputs, term == "Mean.Change"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2("IUCN\nTrend", limits = limits_IUCN, 
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 3)) +
    new_scale("fill") +
    geom_tile(data = filter(all_model_outputs, term == "tmax"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2(expression(psi["temp"]), limits = limits_tmax, 
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 4)) +
    new_scale("fill") +
    geom_tile(data = filter(all_model_outputs, term == "prec"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2(expression(psi["precip"]), limits = limits_prec, 
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 5)) +
    new_scale("fill") +
    geom_tile(data = filter(all_model_outputs, term == "floral"), aes(x = term, y = species, fill = mean)) +
    scale_fill_gradient2(expression(psi["floral"]), limits = limits_floral, 
                         low = mypalette[1], mid = mypalette[round(length(mypalette)/2)], high = mypalette[length(mypalette)], guide = guide_colorbar(order = 6))  +
    
    theme_cowplot() +
    xlab("") + ylab("") +
    scale_x_discrete(limits = c('era', "percent_change", "Mean.Change", 'tmax','prec','floral'),
                     labels = c(expression(psi["era"]), "Prop.\nChange", "IUCN\nTrend", expression(psi["temp"]), 
                                expression(psi["precip"]), expression(psi["floral"]))) +
    geom_text(data = filter(all_model_outputs, term %in% c("percent_change", "Mean.Change")), 
              aes(x = term, y = species, label = round(mean, 2), 
                  colour = col),
              size = 5) +
    scale_colour_manual(values = c('black', 'white'), guide = NULL)+
    theme(axis.text.y = element_text(size = 20, hjust = 0),
          axis.text.x = element_text(size = 22, vjust = -1.2),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 22))
  
  ## add winners and losers only for the 100 resolution
  
  if(res == 100){
    
    heatplot_pred <- heatplot_pred_partial + 
      geom_text(data = filter(all_model_outputs, !(term %in% c("percent_change", "Mean.Change")) & sig == TRUE),
                aes(x = term, y = species, label = paste0("(",round(`X2.5.`,2),",",round(`X97.5.`,2),")", " ", w_l),
                    colour = col),
                size = 5) 
  }else{
    
    heatplot_pred <- heatplot_pred_partial + 
      geom_text(data = filter(all_model_outputs, !(term %in% c("percent_change", "Mean.Change")) & sig == TRUE),
                aes(x = term, y = species, label = paste0("(",round(`X2.5.`,2),",",round(`X97.5.`,2),")"),
                    colour = col),
                size = 5)
    
  }
  
  
  ## save the heatplot
  
  ggsave(heatplot_pred, filename = paste0('heatplot_predictors_',res,'.jpeg'),dpi = 300, height = 15, width = 13)
  
  ggsave(heatplot_pred, filename = paste0('heatplot_predictors_',res,'.pdf'), height = 15, width = 13)
  
  
}







