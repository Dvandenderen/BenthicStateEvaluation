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

indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf",
                      "SoS","AMBI","M-AMBI","BENTIX","Dm'","mTDI",
                      "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list
  
# get table with the indicator's slope coefficient per gradient
#  if AIC is smaller than model without trawling
oview <- data.frame(areaid_Trawl,matrix(data=NA,nrow=14,ncol=18))
  for (j in 1:length(areaid_Trawl)){
    for(p in 1:length(indicator_list)){
      td <- subset(dat,dat$RegName == areaid_Trawl[j])
      if(sum(!(is.na(td[,indicator_list[p]]))) > (nrow(td)/2)){
      if(sum(td[,indicator_list[p]],na.rm=T) >0){
      mod1 <- lm(td[,indicator_list[p]] ~ td$Ltrawl)
      mod2 <- lm(td[,indicator_list[p]] ~1)
      if(AIC(mod1) < AIC(mod2)){ 
        oview[j,p+1] <- round(mod1$coefficients[2],digits = 2)
    }}}}}
  
colnames(oview)[2:19] <- indicator_list
  
# set plotting
pdf(file="Output/Trawl_gradients.pdf",width = 9, height = 7)

par(mfrow=c(3,5),oma = c(2, 1, 2, 3) + 0.1,mar=c(1,1,4,1)+0.1)
colsc <- c("blue","blue","purple","purple","purple","red","red","red","red","black","black")
ltyp  <- c(1,5,1,3,5,1,3,4,5,1,5) 
regname <- c("Italian EEZ sand (1)", "Italian EEZ mud (2)", "Dutch EEZ sand (3)","Dogger Bank (4)",
             "Fladen Ground (5)","Long Forties (6)","Silver Pit (7)","Thames (8)","Iberian Coast sand (9)",
             "Iberian Coast mud (10)","Gotland (11)","Polish EEZ (12)","Flemish Cap (13)","Sellafield (14)")
              
for (j in 1:length(areaid_Trawl)){
  # get one area and gradient
  td <- subset(dat,dat$RegName == areaid_Trawl[j])
  Trawling_intensity <- seq(from=min(td$Ltrawl),to=max(td$Ltrawl), length.out =100)
  TR <- (10^(Trawling_intensity)-1)
    
    # make empty plot
    plot(0,0,ylim=c(0,1.5),las=1,xlab="",ylab="",col="white", xlim=c(0,max(TR)),yaxt="n",xaxt="n",
         main=paste(regname[j]),cex.main=1)
    #text(TR[93],1.3, paste("reg.", areaid_Trawl[j],sep=" "),lwd=2)
    
    axis(1,c(0,round(max(TR))/2,round(max(TR))))
    
    if (j %in% c(1,6,11)){
      axis(2,c(0,0.5,1,1.5),las=1)
    } else{
      axis(2,c(0,0.5,1,1.5),c("","","",""),las=1)
    }
    
    # richness
    mod1 <- lm(R ~ Ltrawl,data=td)
    mod2 <- lm(R ~1, data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(Ltrawl = Trawling_intensity))
      out[out<0] <- 0
      lines(out/out[1]~TR,type="l",col=colsc[1],lty=ltyp[1],lwd=2)
    }
    
    # biomass
    if(sum(!(is.na(td$B))) > (nrow(td)/2)){
    mod1 <- lm((B) ~ Ltrawl,data=td)
    mod2 <- lm((B) ~1, data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl = Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[2],lty=ltyp[2],lwd=2)  
      }}
    
  # diversity indicators
    if(sum(!(is.na(td$`Dm'`))) > (nrow(td)/2)){
      mod1 <- lm(`Dm'` ~ Ltrawl, data=td)
      mod2 <- lm(`Dm'` ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[3],lty=ltyp[3],lwd=2)
      }}
    
    if(sum(!(is.na(td$`H'`))) > (nrow(td)/2)){
      mod1 <- lm(`H'` ~ Ltrawl, data=td)
      mod2 <- lm(`H'` ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[4],lty=ltyp[4],lwd=2)
      }}
   
    if(sum(!(is.na(td$SI))) > (nrow(td)/2)){
      mod1 <- lm(SI ~ Ltrawl, data=td)
      mod2 <- lm(SI ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[5],lty=ltyp[5],lwd=2)
      }} 
    
  # indicators with sensitivity component
    if(sum(!(is.na(td$`M-AMBI`))) > (nrow(td)/2)){
      mod1 <- lm(`M-AMBI` ~ Ltrawl,data=td)
      mod2 <- lm(`M-AMBI`  ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl = Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[6],lty=ltyp[6],lwd=2)
      }}
    
    if(sum(!(is.na(td$BENTIX))) > (nrow(td)/2)){
      mod1 <- lm(BENTIX ~ Ltrawl, data=td)
      mod2 <- lm(BENTIX ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[7],lty=ltyp[7],lwd=2)
      }}
    
    if(sum(!(is.na(td$DKI))) > (nrow(td)/2)){
      mod1 <- lm(DKI ~ Ltrawl, data=td)
      mod2 <- lm(DKI ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[8],lty=ltyp[8],lwd=2)
      }}
    
    if(sum(!(is.na(td$Lm))) > (nrow(td)/2)){
      mod1 <- lm(Lm ~ Ltrawl,data=td)
      mod2 <- lm(Lm ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[9],lty=ltyp[9],lwd=2)
      }}
    
    if(sum(!(is.na(td$SoS))) > (nrow(td)/2)){
      mod1 <- lm(SoS ~ Ltrawl,data=td)
      mod2 <- lm(SoS ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[10],lty=ltyp[10],lwd=2)
      }}
    
    if(sum(!(is.na(td$Lf))) > (nrow(td)/2)){
      mod1 <- lm(Lf ~ Ltrawl,data=td)
      mod2 <- lm(Lf ~ 1,data=td)
      if(AIC(mod1) < AIC(mod2)){
        out <- predict(mod1,newdata = data.frame(Ltrawl =Trawling_intensity))
        out[out<0] <- 0
        lines(out/out[1]~TR,type="l",col=colsc[11],lty=ltyp[11],lwd=2)
      }}
  }

# get legend
  legend(par('usr')[2]+4, par('usr')[4]+0.2, bty='n', xpd=NA,
         legend=c("R","B","Dm'","H'","SI","M-AMBI","BENTIX","DKI","Lm","SoS","Lf"),
         col=colsc,lty=ltyp,y.intersp = 0.9,cex = 1.2,box.col = "white",lwd=1.5)

  dev.off()
  
rm(list=setdiff(ls(), c("indic","stations")))
  