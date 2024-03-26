# ------------------------------------------------------------------------------
# examine the effect of the other pressures - using on GAM
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

indicator_list <- c("B","A","R","H'","SI","IS","Lm","Lf",
                    "SoS","AMBI","M-AMBI","BENTIX","Dm","mTDI",
                    "TDI","pTDI","mT","DKI")
colnames(dat)[12:29] <- indicator_list
  
# get table which indicators are showing a response
oview <- data.frame(areaid_other,matrix(data=NA,nrow=3,ncol=18))
  
  for (j in 1:length(areaid_other)){
    for(p in 1:length(indicator_list)){
      td <- subset(dat,dat$RegName == areaid_other[j])
      if(sum(!(is.na(td[,indicator_list[p]]))) > (nrow(td)/2)){
        if(sum(td[,indicator_list[p]],na.rm=T) >0){
          mod1 <- gam(td[,indicator_list[p]] ~ s(other_pres, k=3),data=td)
          mod2 <- gam(td[,indicator_list[p]] ~ 1,data=td)
          if(AIC(mod1) < AIC(mod2)){ 
            oview[j,p+1] <- 1
          }}}}}
  
  colnames(oview)[2:19] <- indicator_list
  
# create the plot
pdf(file="Output/NonTrawl_gradients.pdf",width = 9, height = 7/2.5)
  
par(mfrow=c(1,5),oma = c(2, 1, 2, 3) + 0.1,mar=c(4,1,1,1)+0.1)
colsc <- c("blue","blue","blue","purple","purple","purple","purple","red","red","red","red",
                   "black","black")
ltyp  <- c(1,5,3,1,3,4,5,1,3,4,5,1,5) 
xnam  <- c("mg O2/L","% total N","Pollution index")
mnam  <- c("Gulf of Finland (15)","Saronikos Gulf (16)", "Vigo Estuary (17)")

for (j in 1:3){
  
  idx <- c(100,1,1) # gradient 1, high pressure value is good
  
  # get one area and gradient
  td <- subset(dat,dat$RegName == areaid_other[j])
  pressure <- seq(from=min(td$other_pres),to=max(td$other_pres),
                  length.out =100)
  
  # make empty plot
  plot(0,0,ylim=c(0,1.5), main=mnam[j],las=1,xlab=xnam[j],ylab="",col="white",
       xlim=c(min(pressure),max(pressure)),yaxt="n",cex.main=0.95)
  
  if (j %in% c(1)){
    axis(2,c(0,0.5,1,1.5),las=1)
  } else{
    axis(2,c(0,0.5,1,1.5),c("","","",""),las=1)
  }
  
  # richness
  mod1 <- gam(R ~ s(other_pres, k=3),data=td)
  mod2 <- gam(R ~ 1,data=td)
  if(AIC(mod1) < AIC(mod2)){
    out <- predict(mod1,newdata = data.frame(other_pres=pressure))
    out[out<0] <- 0
    lines(out/out[idx[j]]~pressure,type="l",col=colsc[1],lty=ltyp[1],lwd=2)
  }
  
  # biomass
  if(sum(!(is.na(td$B))) > (nrow(td)/2)){
    mod1 <- gam(B ~ s(other_pres, k=3),data=td)
    mod2 <- gam(B ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[2],lty=ltyp[2],lwd=2)
    }}
  
  # abundance
  if(sum(!(is.na(td$A))) > (nrow(td)/2)){
    mod1 <- gam(A ~ s(other_pres, k=3),data=td)
    mod2 <- gam(A ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[3],lty=ltyp[3],lwd=2)
    }}
  
  # diversity indicators
  if(sum(!(is.na(td$Dm)))  > (nrow(td)/2)){
    mod1 <- gam(Dm ~ s(other_pres, k=3),data=td)
    mod2 <- gam(Dm ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[4],lty=ltyp[4],lwd=2)
    }}
  
  if(sum(!(is.na(td$H)))  > (nrow(td)/2)){
    mod1 <- gam(H ~ s(other_pres, k=3),data=td)
    mod2 <- gam(H ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[5],lty=ltyp[5],lwd=2)
    }}
  
  if(sum(!(is.na(td$SI)))  > (nrow(td)/2)){
    mod1 <- gam(SI ~ s(other_pres, k=3),data=td)
    mod2 <- gam(SI ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[6],lty=ltyp[6],lwd=2)
    }}
  
  if(sum(!(is.na(td$IS)))  > (nrow(td)/2)){
    mod1 <- gam(IS ~ s(other_pres, k=3),data=td)
    mod2 <- gam(IS ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[7],lty=ltyp[7],lwd=2)
    }}
  
  # indicators with sensitivity component
  if(sum(!(is.na(td$M_AMBI)))  > (nrow(td)/2)){
    mod1 <- gam(M_AMBI ~ s(other_pres, k=3),data=td)
    mod2 <- gam(M_AMBI ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[8],lty=ltyp[8],lwd=2)
    }}
  
  if(sum(!(is.na(td$BENTIX)))  > (nrow(td)/2)){
    mod1 <- gam(BENTIX ~ s(other_pres, k=3),data=td)
    mod2 <- gam(BENTIX ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[9],lty=ltyp[9],lwd=2)
    }}
   
  if(sum(!(is.na(td$DKI)))  > (nrow(td)/2)){
    mod1 <- gam(DKI ~ s(other_pres, k=3),data=td)
    mod2 <- gam(DKI ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[10],lty=ltyp[10],lwd=2)
    }}
  
  if(sum(!(is.na(td$Lm)))  > (nrow(td)/2)){
    mod1 <- gam(Lm ~ s(other_pres, k=3),data=td)
    mod2 <- gam(Lm ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[11],lty=ltyp[11],lwd=2)
    }}
  
  # sensitive fraction
  if(sum(!(is.na(td$SoS)))  > (nrow(td)/2)){
    mod1 <- gam(SoS ~ s(other_pres, k=3),data=td)
    mod2 <- gam(SoS ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      if(out[100] < 0 & j == 1){out <- out - out[100] + 1}
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[12],lty=ltyp[12],lwd=2)
    }}
  
  if(sum(!(is.na(td$Lf)))  > (nrow(td)/2)){
    mod1 <- gam(Lf ~ s(other_pres, k=3),data=td)
    mod2 <- gam(Lf ~ 1,data=td)
    if(AIC(mod1) < AIC(mod2)){
      out <- predict(mod1,newdata = data.frame(other_pres=pressure))
      out[out<0] <- 0
      lines(out/out[idx[j]]~pressure,type="l",col=colsc[13],lty=ltyp[13],lwd=2)
    }}
  
}
  # get legend
  legend(par('usr')[2]+4, par('usr')[4]+0.2, bty='n', xpd=NA,
         legend=c("R","B","A","Dm'","H'","SI","IS","M-AMBI","BENTIX","DKI","Lm","SoS","Lf"),
         col=colsc,lty=ltyp,y.intersp = 0.9,cex = 1.2,box.col = "white",lwd=1.5)
  
  plot.new()
  
  dev.off()
  #save 9 x 7 landscape