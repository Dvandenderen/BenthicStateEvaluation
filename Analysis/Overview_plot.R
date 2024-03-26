# ------------------------------------------------------------------------------
#  examine the effect of the fishing pressure - linear model log10(fishing pressure)
#  response shows the relative response between 0 and 1
#
#  note: not all gradients have values for all indicators   
# ------------------------------------------------------------------------------
dat <- indic # get the data

# get stations and add unique ID number (matching manuscript)
stat <- stations
ord  <- data.frame(loc = c("AS1","AS2","CO","DB","FG","PH","SP","TH","NIC1","NIC2","Gotland","OxyTrawl",
          "FC","SEL","Finland","Saronikos","Vigo"), id= c(1:17))
dat <- cbind(dat, ord[match(dat$Name,ord$loc), c(2)])
colnames(dat)[ncol(dat)] <- "RegName"

# log transform trawling intensity
dat$Ltrawl   <- log10(dat$Trawling_intensity+1)
areaid_Trawl <-  c(1:14)  

# exclude unused indicators and get proper names
exclude <- which(colnames(dat) %in% 
                   c("Abundance","Richness","Biomass","AMBI_bio",
                     "M_AMBI_abund_plusref","Ab_mTDI","Ab_TDI",
                     "Ab_pTDI","Ab_mT","logab_mTDI","logab_TDI",
                     "logab_pTDI","logab_mT","logbiom_mTDI",
                     "logbiom_TDI","logbiom_pTDI","logbiom_mT"))
dat <- dat[,-exclude]

indicator_list <- c("B","A","R","H","SI","IS","Lm","Lf",
                      "SoS","AMBI","M_AMBI","BENTIX","Dm","mTDI",
                      "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list
  
# get table with the indicator's slope coefficient per gradient
#  if AIC is smaller than model without trawling
oview <- data.frame(areaid_Trawl,matrix(data=NA,nrow=14,ncol=18))
  for (j in 1:length(areaid_Trawl)){
    for(p in 1:length(indicator_list)){
      td <- subset(dat,dat$RegName == areaid_Trawl[j])
      Trawling_intensity <- seq(from=min(td$Ltrawl),to=max(td$Ltrawl), length.out =100)
      if(sum(!(is.na(td[,indicator_list[p]]))) > (nrow(td)/2)){
      if(sum(td[,indicator_list[p]],na.rm=T) >0){
        mod1 <- lm(paste(indicator_list[p], "~", "Ltrawl",sep=""),data=td)
        mod2 <- lm(paste(indicator_list[p], "~", "1",sep=""),data=td)
        oview[j,p+1] <- -5
      if(AIC(mod1) < AIC(mod2)){ 
        out <- predict(mod1,newdata = data.frame(Ltrawl = Trawling_intensity))
        out[out<0] <- 0
        oview[j,p+1] <- out[100]/out[1]
    }}}}}
  
colnames(oview)[2:19] <- indicator_list

# get other gradients
# ------------------------------------------------------------------------------
dat <- indic # get the data

# get stations and add unique ID number (matching manuscript)
stat <- stations
ord  <- data.frame(loc = c("AS1","AS2","CO","DB","FG","PH","SP","TH","NIC1","NIC2","Gotland","OxyTrawl",
                           "FC","SEL","Finland","Saronikos","Vigo"), id= c(1:17))
dat <- cbind(dat, ord[match(dat$Name,ord$loc), c(2)])
colnames(dat)[ncol(dat)] <- "RegName"

# select the three other gradients
areaid_other <-  c(15:17)  
dat$other_pres <- NA 
dat$other_pres <- ifelse(dat$RegName %in% c(15,16),dat$Eutrophication_intensity,dat$other_pres)
dat$other_pres <- ifelse(dat$RegName %in% c(17),dat$Pollution_intensity,dat$other_pres)

# exclude unused indicators and get proper names
exclude <- which(colnames(dat) %in% 
                   c("Abundance","Richness","Biomass","AMBI_bio",
                     "M_AMBI_abund_plusref","Ab_mTDI","Ab_TDI",
                     "Ab_pTDI","Ab_mT","logab_mTDI","logab_TDI",
                     "logab_pTDI","logab_mT","logbiom_mTDI",
                     "logbiom_TDI","logbiom_pTDI","logbiom_mT"))
dat <- dat[,-exclude]

indicator_list <- c("B","A","R","H","SI","IS","Lm","Lf",
                    "SoS","AMBI","M_AMBI","BENTIX","Dm","mTDI",
                    "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list

# get table which indicators are showing a response
oview2 <- data.frame(areaid_other,matrix(data=NA,nrow=3,ncol=18))

for (j in 1:length(areaid_other)){
  for(p in 1:length(indicator_list)){
    td <- subset(dat,dat$RegName == areaid_other[j])
    pressure <- seq(from=min(td$other_pres),to=max(td$other_pres),length.out =100)
    if(sum(!(is.na(td[,indicator_list[p]]))) > (nrow(td)/2)){
      if(sum(td[,indicator_list[p]],na.rm=T) >0){
        mod1 <- gam(td[,indicator_list[p]] ~ s(other_pres, k=3),data=td)
        mod2 <- gam(td[,indicator_list[p]] ~ 1,data=td)
        oview2[j,p+1] <- -5
        if(AIC(mod1) < AIC(mod2)){ 
          out <- predict(mod1,newdata = data.frame(other_pres=pressure))
          out[out<0] <- 0
          if(j == 1){oview2[j,p+1] <- out[1]/out[100]}
          if(j %in% c(2,3)){oview2[j,p+1] <- out[100]/out[1]}
        }}}}}

colnames(oview2)[2:19] <- indicator_list

# set plotting
colnames(oview)[1] <- "region"
colnames(oview2)[1] <- "region"
oview <- rbind(oview,oview2)
colnames(oview)[2:19] <- c("B","A","R","H'","SI","IS","Lm","Lf",
                           "SoS","AMBI","M-AMBI","BENTIX","Dm'","mTDI",
                           "TDI","pTDI","mT","DKI")
oview <- oview[,c("region","R","B","A","Dm'","H'","SI","IS",
                  "AMBI","M-AMBI","BENTIX","DKI","TDI","mTDI","mT","Lm",
                  "pTDI","SoS","Lf")]

# reverse AMBI and mT
oview$AMBI <- ifelse(oview$AMBI>0,1/oview$AMBI,oview$AMBI)
oview$mT <- ifelse(oview$mT>0,1/oview$mT,oview$mT)

# group together
oview[oview > 1] <- -1
oview[oview >= 0 & oview < 0.25] <- 4
oview[oview >= 0.25 & oview < 0.5] <- 3
oview[oview >= 0.5 & oview < 0.75] <- 2
oview[oview >= 0.75 & oview < 1] <- 1
oview[oview == -5] <- 0 
oview$region <- 1:17

tot <- data.frame(rep(oview$region,each=18),rep(colnames(oview)[2:19],17))
ind <- c()
for(j in 1:17){
  ind_sub <- as.numeric(oview[j,2:19])
  ind <- c(ind,ind_sub)
}
tot <- cbind(tot,ind)
colnames(tot) <- c("Region","Indicator","Effect")

tot$Region <- factor(tot$Region, levels = unique(tot$Region))
tot$Indicator <- factor(tot$Indicator, levels = rev(unique(tot$Indicator)))
tot$Effect <- factor(tot$Effect, levels = c("-1","0","1","2","3","4"))

# set colors
coldat <- data.frame(cat = as.factor(c(-1,0,1,2,3,4)), 
                     col = c("#fcae91","white","#bdd7e7","#6baed6","#3182bd","#08519c"))

tot <- cbind(tot, coldat[match(tot$Effect,coldat$cat), c(2)])
colnames(tot)[ncol(tot)] <- "coll"
tot$coll[is.na(tot$coll)] <- "grey"
  
rankCV <- data.frame(region= c(1:17),
                     rank = c(13,14,8,4,1,2,3,11,
                              10,12,7,5,6,9,15,16,17))
rankCV <- rankCV[order(rankCV$rank),] 


tot <- cbind(tot, rankCV[match(tot$Region,rankCV$region), c(2)])
colnames(tot)[ncol(tot)] <- "Gradient"

p <- ggplot(tot, aes(x = Gradient, y = Indicator)) + 
  geom_raster(aes(fill=Effect))+ scale_fill_manual(values = c("-1" =  "#fcae91",
                                                              "0" =  "#ffffe5",
                                                              "1" = "#9ecae1",
                                                              "2" ="#6baed6",
                                                              "3" = "#3182bd",
                                                              "4" = "#08519c"),
                                                   na.value="#f0f0f0",
                                                   labels=c('increase', 'no effect',
                                                            '0-25% decline',"25-50% decline",
                                                            "50-75% decline","75-100% decline"))+
  theme (legend.title=element_text(size=10),legend.position = "bottom",panel.background = element_blank())
p <- p + scale_x_continuous(breaks=1:17,labels= rankCV$region)

pdf("Output/Overview_all_gradients.pdf",width = 5,height = 5.5)   
print(p)
dev.off()

rm(list=setdiff(ls(), c("indic","stations")))
  