
# ------------------------------------------------------------------------------
# make a correlation plot based on the average correlation across regions
# ------------------------------------------------------------------------------

dat <- indic # get output 
dat <- subset(dat,dat$Trawling_intensity >=0) # select all trawling gradients

# reverse AMBI and mT
dat$AMBI_abund  <- dat$AMBI_abund * -1 + max(dat$AMBI_abund,na.rm=T)
dat$biom_mT     <- dat$biom_mT * -1 + max(dat$biom_mT,na.rm=T)

# log transform abundance and biomass
dat$Abundance_nc <- log10(dat$Abundance_nc+1) # set abundance to log10(x+1)
dat$Biomass_nc   <- log10(dat$Biomass_nc+1)   # set biomass to log10(x+1)

# change names to match with paper
indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf","SoS","AMBI",
                    "M-AMBI","BENTIX","Dm'","mTDI","TDI","pTDI","mT","DKI")
colnames(dat)[c(15:25,28,29,34:37,46)] <- indicator_list

# ------------------------------------------------------------------------------
# get correlation
# ------------------------------------------------------------------------------
Mtot <- list()
areaid <- unique(dat$Name)

for (j in 1:length(areaid)) {
  #Select the dataset
  td <- subset(dat, dat$Name == areaid[j])
  
  # Select only indicators and remove empty cols
  td <-
    remove_empty(td[, which(colnames(td) %in% indicator_list)], which = "cols")
  
  # Remove indicators with one value
  td <- td %>% select_if( ~ sd(., na.rm = TRUE) != 0)
  
  # Compute the correlation
  M  <- cor(td, use = "pairwise.complete.obs")
  
  # Format correlation into long-format matrix
  M <- M %>%
    as.data.frame() %>%
    rownames_to_column("var1") %>%
    gather(var2, cor,-var1) %>%
    mutate(Name = areaid[j])
  
  Mtot[[j]] <- M
}
Mtot <- map_dfr(Mtot, as.data.frame)

# Summarise mean of correlation across datasets
Mtot_summary <- Mtot %>% 
  group_by(var1, var2) %>%
  summarise(mean_cor = mean(cor, na.rm = TRUE), sd_cor = sd(cor, na.rm = TRUE)) %>%
  ungroup()

mean_cor <- Mtot_summary %>%
  select(-sd_cor)  %>%
  spread(var2, mean_cor, fill = 0) %>%
  column_to_rownames("var1")
mean_cor[mean_cor == 1] <- NA

# plot
pdf(file="Output/Correlation_matrix.pdf",width = 6.5,height = 6.5)
tr <- ggcorrplot(
      mean_cor,
      hc.order = TRUE,
      show.diag = F,
      method = "circle",
      type = "lower",
      outline.color = "white",
      legend.title = "Correlation",
      colors = c("#E46726", "white", "#6D9EC1"),
      ggtheme = ggplot2::theme_minimal())
print(tr)
dev.off()

# get values for supplement
mean_cor <- mean_cor[c(11,5,15,6,16,7,8,4,2,3,1,18,13,12,14,10,9,17),] #need to be re-ordered in excel
write.csv(mean_cor,"Output/appendixS3_mean_correlation.csv")

sd_cor <- Mtot_summary %>%
  select(-mean_cor)  %>%
  spread(var2, sd_cor, fill = 0) %>%
  column_to_rownames("var1") 
sd_cor[sd_cor==0] <- NA
sd_cor <- sd_cor[c(11,5,15,6,16,7,8,4,2,3,1,18,13,12,14,10,9,17),] #need to be re-ordered in excel
write.csv(sd_cor,"Output/appendixS3_sd_correlation.csv")

rm(list=setdiff(ls(), c("indic","stations")))

