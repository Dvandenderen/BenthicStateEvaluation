
# ------------------------------------------------------------------------------
#  estimate mean response across 14 trawl gradient studies
# ------------------------------------------------------------------------------
dat <- indic # get the data

# exclude all indicator estimates not used in main analysis
# biomass/abundance/richness is without cephalopods, i.e. Biomass_nc, Abundance_nc & Richness_nc
exclude <- which(colnames(dat) %in% 
                           c("Abundance","Richness","Biomass","AMBI_bio",
                             "M_AMBI_abund_plusref","Ab_mTDI","Ab_TDI",
                             "Ab_pTDI","Ab_mT","logab_mTDI","logab_TDI",
                             "logab_pTDI","logab_mT","logbiom_mTDI",
                             "logbiom_TDI","logbiom_pTDI","logbiom_mT"))

dat <- dat[,-exclude]

# select all trawl gradients
dat <- subset(dat,dat$Trawling_intensity >=0)

# set the names
indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf","SoS","AMBI","M-AMBI",
                    "BENTIX","Dm'","mTDI","TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list

# re-order the indicators
indicator_reorder    <- c("R","B","A","Dm'","H'","SI","IS","AMBI","M-AMBI","BENTIX","DKI",
                          "TDI","mTDI","mT","Lm","pTDI","SoS","Lf")

#-------------------------------------------------------------------------------
# calculate mean indicator value of reference stations
# ------------------------------------------------------------------------------
ref_stat <- subset(dat, dat$Reference_stations == 'yes')
  
# get mean
Mean_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat[,"Name"]),FUN=mean,na.rm=T)
colnames(Mean_Ctrl)[1] <- "Name"

# get sd
SD_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Ctrl)[1] <- "Name"

# get data count
cou <- function(x) sum(!is.na(x))
n_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=cou)
colnames(n_Ctrl)[1] <- "Name"
n_Ctrl[n_Ctrl == 0] <- NA

#-------------------------------------------------------------------------------
# calculate mean indicator value of high pressure stations
# ------------------------------------------------------------------------------
high_stat <- subset(dat, dat$High_impact == 'yes')
  
Mean_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=mean,na.rm=T)
colnames(Mean_Intr)[1] <- "Name"
Mean_Intr[] <- lapply(Mean_Intr, function(x) ifelse(x<0.01, 0.01, x))

SD_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Intr)[1] <- "Name"

cou <- function(x) sum(!is.na(x))
n_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=cou)
colnames(n_Intr)[1] <- "Name"
n_Intr[n_Intr == 0] <- NA

# ------------------------------------------------------------------------------
# estimate min, mean and max at reference and high pressure stations
# only possible for all with SAR values
# ------------------------------------------------------------------------------
print("low fished stations")
print(ref_stat %>% 
  group_by(Name) %>% 
  summarise(min_trawl    = min(Trawling_intensity),
            mean_trawl = mean(Trawling_intensity),
            max_trawl    = max(Trawling_intensity),
            stations     = length(Trawling_intensity)) %>%
  ungroup())

print("high fished stations")
print(high_stat %>% 
        group_by(Name) %>% 
        summarise(min_trawl    = min(Trawling_intensity),
                  mean_trawl = mean(Trawling_intensity),
                  max_trawl    = max(Trawling_intensity),
                  stations     = length(Trawling_intensity)) %>%
        ungroup())

# ------------------------------------------------------------------------------
# use metafor to estimate the mean response to trawling
# ------------------------------------------------------------------------------
metan <- matrix(data=NA,nrow=18,ncol=3)
for (j in 1:18) {
  parmod    <-
    escalc(
      measure = "ROM",
      m1i = Mean_Intr[, j + 1],
      m2i = Mean_Ctrl[, j + 1],
      sd1i = SD_Intr[, j + 1],
      sd2i = SD_Ctrl[, j + 1],
      n1i = n_Intr[, j + 1],
      n2i = n_Ctrl[, j + 1],
      var.names = c("ES", "Var"),
      digits = 4
    )
  
  mod1       <- rma(parmod[, 1], vi = parmod[, 2], method = "REML")
  metan[j, 1] <- mod1$b
  metan[j, 2] <- mod1$ci.lb
  metan[j, 3] <- mod1$ci.ub
}

# create the plot
  metan <- as.data.frame(metan)
  metan$nam <- indicator_reorder
  metan$nam <- factor(metan$nam, levels = metan$nam)
  
# reverse AMBI and mT
  idx <- which(metan$nam %in% c("AMBI","mT"))
  metan[idx,1:3] <-  metan[idx,1:3]*-1
  
  tt <- SD_Intr[,2:19] + n_Intr[,2:19] + n_Ctrl[,2:19] + SD_Ctrl[,2:19]
  allgear <- ggplot() + geom_point(data=metan,aes(x=nam,y=V1) )+
    geom_errorbar(data=metan, aes(x=1:18,ymin = V3, ymax = V2),width = 0.3,linewidth=0.2) + theme_classic() +
    geom_hline(yintercept = 0,linetype="dashed",linewidth=0.2,col="grey") + #scale_x_discrete(labels = NULL)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.y =  element_text(size=11)) + 
    coord_cartesian(ylim = c(-1.8 , 0.2)) + labs(y="% decline",x="") +
    annotate("text",x=c(1:18),y=rep(-2.2,18),label = as.character(colSums(!(is.na(tt)))),size=3)+
    scale_y_discrete(limits=c(0,log(0.8),log(0.6),log(0.4),log(0.2)),labels = c("0","20","40","60","80"))+
    annotate("text",x=17.5,y=0.6,label="All")
 
# ------------------------------------------------------------------------------
# now do the same for grabs/cores
# ------------------------------------------------------------------------------
dat <- indic # get the data
  
# exclude all indicator estimates not used in main analysis
# biomass/abundance/richness is without cephalopods, i.e. Biomass_nc, Abundance_nc & Richness_nc
exclude <- which(colnames(dat) %in% 
                   c("Abundance","Richness","Biomass","AMBI_bio",
                     "M_AMBI_abund_plusref","Ab_mTDI","Ab_TDI",
                     "Ab_pTDI","Ab_mT","logab_mTDI","logab_TDI",
                     "logab_pTDI","logab_mT","logbiom_mTDI",
                     "logbiom_TDI","logbiom_pTDI","logbiom_mT"))
  
dat <- dat[, -exclude]
  
# select all trawl gradients with grabs/cores
dat <- subset(dat,dat$Trawling_intensity >=0)
dat <- subset(dat,(dat$Name %in% c("CO","DB", "FG","Gotland", "OxyTrawl", "PH", "SEL", "SP", "TH")))
  
# set the names
indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf",
                      "SoS","AMBI","M-AMBI","BENTIX","Dm'","mTDI",
                      "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list
  
indicator_reorder <- c("R","B","A","Dm'","H'","SI","IS","AMBI","M-AMBI","BENTIX","DKI",
                       "TDI","mTDI","mT","Lm","pTDI","SoS","Lf")

#-------------------------------------------------------------------------------
# calculate mean indicator value of reference stations
# ------------------------------------------------------------------------------
ref_stat <- subset(dat, dat$Reference_stations == 'yes')

# get mean
Mean_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat[,"Name"]),FUN=mean,na.rm=T)
colnames(Mean_Ctrl)[1] <- "Name"

# get sd
SD_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Ctrl)[1] <- "Name"

# get data count
cou <- function(x) sum(!is.na(x))
n_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=cou)
colnames(n_Ctrl)[1] <- "Name"
n_Ctrl[n_Ctrl == 0] <- NA

#-------------------------------------------------------------------------------
# calculate mean indicator value of high pressure stations
# ------------------------------------------------------------------------------
high_stat <- subset(dat, dat$High_impact == 'yes')

Mean_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=mean,na.rm=T)
colnames(Mean_Intr)[1] <- "Name"
Mean_Intr[] <- lapply(Mean_Intr, function(x) ifelse(x<0.01, 0.01, x))

SD_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Intr)[1] <- "Name"

cou <- function(x) sum(!is.na(x))
n_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=cou)
colnames(n_Intr)[1] <- "Name"
n_Intr[n_Intr == 0] <- NA

# ------------------------------------------------------------------------------
# now use metafor to estimate the mean response to trawling
# ------------------------------------------------------------------------------
metan <- matrix(data=NA,nrow=18,ncol=3)
for (j in 1:18) {
  parmod    <-
    escalc(
      measure = "ROM",
      m1i = Mean_Intr[, j + 1],
      m2i = Mean_Ctrl[, j + 1],
      sd1i = SD_Intr[, j + 1],
      sd2i = SD_Ctrl[, j + 1],
      n1i = n_Intr[, j + 1],
      n2i = n_Ctrl[, j + 1],
      var.names = c("ES", "Var"),
      digits = 4
    )
  
  mod1       <- rma(parmod[, 1], vi = parmod[, 2], method = "REML")
  metan[j, 1] <- mod1$b
  metan[j, 2] <- mod1$ci.lb
  metan[j, 3] <- mod1$ci.ub
}

# create the plot
metan <- as.data.frame(metan)
metan$nam <- indicator_reorder
metan$nam <- factor(metan$nam, levels = metan$nam)

# reverse AMBI and mT
idx <- which(metan$nam %in% c("AMBI","mT"))
metan[idx,1:3] <-  metan[idx,1:3]*-1

tt <- SD_Intr[,2:19] + n_Intr[,2:19] + n_Ctrl[,2:19] + SD_Ctrl[,2:19]
core <- ggplot() +geom_point(data=metan,aes(x=nam,y=V1) )+
  geom_errorbar(data=metan, aes(x=1:18,ymin = V3, ymax = V2),width = 0.3,linewidth=0.2) + theme_classic() +
  geom_hline(yintercept = 0,linetype="dashed",linewidth=0.2,col="grey") + #scale_x_discrete(labels = NULL)+
  theme(axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y= element_text(size=10),axis.title.y =  element_text(size=11)) + 
  coord_cartesian(ylim = c(-1.8 , 0.2)) + labs(y="% decline",x="")+
  scale_y_discrete(limits=c(0,log(0.8),log(0.6),log(0.4),log(0.2)),labels = c("0","20","40","60","80"))+
  annotate("text",x=c(1:18),y=rep(-2.2,18),label = as.character(colSums(!(is.na(tt)))),size=3)+
  annotate("text",x=16,y=0.6,label="Cores/grabs")

# ------------------------------------------------------------------------------
# now do the same for trawls only
# ------------------------------------------------------------------------------

dat <- indic # get the data

# exclude all indicator estimates not used in main analysis
# biomass/abundance/richness is without cephalopods, i.e. Biomass_nc, Abundance_nc & Richness_nc
exclude <- which(colnames(dat) %in% 
                   c("Abundance","Richness","Biomass","AMBI_bio",
                     "M_AMBI_abund_plusref","Ab_mTDI","Ab_TDI",
                     "Ab_pTDI","Ab_mT","logab_mTDI","logab_TDI",
                     "logab_pTDI","logab_mT","logbiom_mTDI",
                     "logbiom_TDI","logbiom_pTDI","logbiom_mT"))

dat <- dat[, -exclude]

# select all trawl gradients sampled with trawls
dat         <- subset(dat,dat$Trawling_intensity >=0)
dat <- subset(dat,!(dat$Name %in% c("CO","DB", "FG","Gotland", "OxyTrawl", "PH", "SEL", "SP", "TH")))

# set the names
indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf",
                    "SoS","AMBI","M-AMBI","BENTIX","Dm'","mTDI",
                    "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list

indicator_reorder <- c("R","B","A","Dm'","H'","SI","IS","AMBI","M-AMBI","BENTIX","DKI",
                       "TDI","mTDI","mT","Lm","pTDI","SoS","Lf")

#-------------------------------------------------------------------------------
# calculate mean indicator value of reference stations
# ------------------------------------------------------------------------------
ref_stat <- subset(dat, dat$Reference_stations == 'yes')

# get mean
Mean_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat[,"Name"]),FUN=mean,na.rm=T)
colnames(Mean_Ctrl)[1] <- "Name"

# get sd
SD_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Ctrl)[1] <- "Name"

# get data count
cou <- function(x) sum(!is.na(x))
n_Ctrl <- aggregate(list(ref_stat[,indicator_reorder]),by=list(ref_stat$Name),FUN=cou)
colnames(n_Ctrl)[1] <- "Name"
n_Ctrl[n_Ctrl == 0] <- NA

#-------------------------------------------------------------------------------
# calculate mean indicator value of high pressure stations
# ------------------------------------------------------------------------------
high_stat <- subset(dat, dat$High_impact == 'yes')

Mean_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=mean,na.rm=T)
colnames(Mean_Intr)[1] <- "Name"
Mean_Intr[] <- lapply(Mean_Intr, function(x) ifelse(x<0.01, 0.01, x))

SD_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=sd,na.rm=T)
colnames(SD_Intr)[1] <- "Name"

cou <- function(x) sum(!is.na(x))
n_Intr <- aggregate(list(high_stat[,indicator_reorder]),by=list(high_stat$Name),FUN=cou)
colnames(n_Intr)[1] <- "Name"
n_Intr[n_Intr == 0] <- NA

# ------------------------------------------------------------------------------
# now use metafor to estimate the mean response to trawling
# ------------------------------------------------------------------------------
metan <- matrix(data=NA,nrow=18,ncol=3)
for (j in 1:18) {
  if (j != 11) {
    Bycatch    <-
      escalc(
        measure = "ROM",
        m1i = Mean_Intr[, j + 1],
        m2i = Mean_Ctrl[, j + 1],
        sd1i = SD_Intr[, j + 1],
        sd2i = SD_Ctrl[, j + 1],
        n1i = n_Intr[, j + 1],
        n2i = n_Ctrl[, j + 1],
        var.names = c("ES", "Var"),
        digits = 4
      )
    
    mod1       <- rma(Bycatch[, 1], vi = Bycatch[, 2], method = "REML")
    metan[j, 1] <- mod1$b
    metan[j, 2] <- mod1$ci.lb
    metan[j, 3] <- mod1$ci.ub
  }
}

# create the plot
metan <- as.data.frame(metan)
metan$nam <- indicator_reorder
metan$nam <- factor(metan$nam, levels = metan$nam)

# reverse AMBI and mT
idx <- which(metan$nam %in% c("AMBI","mT"))
metan[idx,1:3] <-  metan[idx,1:3]*-1

tt <- SD_Intr[,2:19] + n_Intr[,2:19] + n_Ctrl[,2:19] + SD_Ctrl[,2:19]
trawl <- ggplot() +geom_point(data=metan,aes(x=nam,y=V1) )+
  geom_errorbar(data=metan, aes(x=1:18,ymin = V3, ymax = V2),width = 0.3,linewidth=0.2) + theme_classic() +
  geom_hline(yintercept = 0,linetype="dashed",linewidth=0.2,col="grey")+
  theme(axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y= element_text(size=10),axis.title.y =  element_text(size=11)) + 
  coord_cartesian(ylim = c(-1.8 , 0.2)) + labs(y="% decline",x="")+
  scale_y_discrete(limits=c(0,log(0.8),log(0.6),log(0.4),log(0.2)),labels = c("0","20","40","60","80"))+
  annotate("text",x=c(1:16,17.2,18),y=rep(-2.2,18),label = as.character(colSums(!(is.na(tt)))),size=3)+
  annotate("text",x=17,y=0.6,label="Trawls")

# save the plot
pdf("Output/Metafor.pdf",width=4.5,height = 8)
print(cowplot::plot_grid(
  allgear,core,trawl,
  labels = "AUTO", ncol = 1,
  rel_heights = c(1,1,1)))
dev.off()

rm(list=setdiff(ls(), c("indic","stations")))