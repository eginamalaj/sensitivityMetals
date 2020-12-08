########################################################################################
#     R code for:
#
#     Physiological sensitivity of freshwater macroinvertebrates to heavy
#     metals, Environmental Toxicology and Chemistry; 31: 1754-1764 https://doi.org/10.1002/etc.1868  
#
#     E. Malaj, M. Grote, R.B. Schaefer, W. Brack, and  P.C. von der Ohe
#
#     contact: eginamalaj@gmail.com
#  
#
#######################################################################################
#     Part 1 - Formating
#######################################################################################
#
# rm(list=ls())
#
#
pkg<- c("lme4","reshape","doBy")
#
# lme4: Mixed model,reshape: data manipulation,
# doBy: data manipulation, transformation 
#
# install.packages(pkg)
#
# Call packages
lapply(pkg, require, character.only=T)
#
setwd("C:/Users/Egina/Dropbox/_Work/Github/sensitivityMetals")
#
org_db<-read.csv("Org_File.csv",sep=",", header=T) # Original table
out_db<-read.csv("Outlier.csv",sep=",", header=T) # Original table
#
#
# Restric to HM
db_hm<- org_db[org_db$Type.of.compound== "Inorganic",]
#
# Restrict to Effects, Endpoint, Exposure
db_eff1<- db_hm[db_hm$Effect== "MOR"|
               db_hm$Effect== "ITX",]
#
db_eff2<- db_eff1[db_eff1$Effect.Measurement== "MORT"| 
          db_eff1$Effect.Measurement== "IMBL",]
#
db_eff3<- db_eff2[db_eff2$Endpoint== "LC50"| 
                    db_eff2$Endpoint== "EC50"|
                    db_eff2$Endpoint== "IC50",]
#
# Restrict Exposure Duration
db_exp<- db_eff3[!is.na(db_eff3$Exposure.Duration..Days.),]
db_exp2<- db_exp[db_exp$Exposure.Duration..Days.== 1|
                   db_exp$Exposure.Duration..Days.== 2|
                   db_exp$Exposure.Duration..Days.== 3|
                   db_exp$Exposure.Duration..Days.== 4,]
#               
#
# Restrict it to the media type
db_med<- db_exp2[db_exp2$Media.Type== "FW"|
                   db_exp2$Media.Type== "NR",]
#
# Restrict water conditions (temp and hardness)
#
as.numeric(db_med$Hardness.Mean_calculated)-> db_med$Hardness.Mean_calculated
#
db_con1_t<- db_med[!is.na(db_med$Temperature.Mean_calculated),]
db_con1_t2<- db_con1_t[db_con1_t$Temperature.Mean_calculated>= 5 &
             db_con1_t$Temperature.Mean_calculated<= 35,]
db_con1_t3<- db_med[is.na(db_med$Temperature.Mean_calculated),]
db_con1<- rbind(db_con1_t2, db_con1_t3)
#
db_con2_h<-db_con1[!is.na(db_con1$Hardness.Mean_calculated),]
db_con2_h2<-db_con2_h[db_con2_h$Hardness.Mean_calculated<= 200,]                 
db_con2_h3<-db_con1[is.na(db_con1$Hardness.Mean_calculated),]
db_con2<- rbind(db_con2_h2, db_con2_h3)
#
# Remove duplicates
db_con3<- db_con2[!duplicated(db_con2$Test.Number),]
#
# Table with outliers. Careful, duplicates here.
out_db1<- merge(db_con3, out_db[,c(1,20:21)], by.x="Test.Number", by.y="Test.Number", all.x=T)
out_db2<- out_db1[out_db1$Exclusion=="no",]
#
# Restrict the number of columns
db_rest<- out_db2[,c("Test.Number","HM","Species.Class", "Species.Order",
                    "Species.Family","Species.Genus","Taxon","Reference.Number",
                     "Conc.1..ug.L._free.ion.conver",
                     "Exposure.Duration..Days.","Temperature.Mean_calculated",
                     "Hardness.Mean_calculated")]
#
names(db_rest)<-c("TestNo","HM","Class", "Order","Family","Genus", "Taxon",
                  "ReferenceNo","Conc","Time","Temp","Hard" )
#
# Ratio Time
db_rest$Rt<- ((db_rest$Time)^(-1))/0.5             
#
# Ratio Temperature. If missing replace with 20°C
db_rest$Temp<- replace(db_rest$Temp, is.na(db_rest$Temp),20)
db_rest$RT<- ((3*exp(-0.1099*db_rest$Temp))/0.33)
#
# Ratio Hardness. If missing replace with 50 mg/l CaCO3
db_rest$Hard<- replace(db_rest$Hard, is.na(db_rest$Hard),50)
db_rest$RH<- (12.522*(db_rest$Hard)^(-0.7852))/0.58
#
db_rest$Cc<- (db_rest$Conc)/(db_rest$Rt*db_rest$RH*db_rest$RT)
#
db_rest$log_Cc<-log10(db_rest$Cc)
#
# Restrict to the HM used. 
#
db_rest2<- db_rest[db_rest$HM== "Cd"|
                   db_rest$HM== "Cr"|
                   db_rest$HM== "Cu"|
                   db_rest$HM== "Hg"|
                   db_rest$HM== "Ni"|
                   db_rest$HM== "Pb"|
                   db_rest$HM== "Zn",]
#
hm_tax<- summaryBy(data=db_rest2, Cc~HM+Taxon, FUN=length)
hm_tax2<- summaryBy(data=hm_tax, HM~Taxon, FUN=length)
#
hm_tax3<- merge(db_rest2, hm_tax2, by.x="Taxon",by.y="Taxon",all.x=T)
#
# Restrict it to specific orders
hm_ins<- hm_tax3[hm_tax3$Order=="Anisoptera"|
                   hm_tax3$Order=="Coleoptera"|
                   hm_tax3$Order=="Zygoptera"|
                   hm_tax3$Order=="Trichoptera"|
                   hm_tax3$Order=="Heteroptera"|
                   hm_tax3$Order=="Megaloptera"|
                   hm_tax3$Order== "Odonata",]
#
hm_tax4<- hm_tax3[hm_tax3$HM.length>=3,]
#
hm_all<- rbind(hm_ins,hm_tax4)
#
#
#############################################################################
#             Part 2 - Analysis
#############################################################################
#
spe_v<- summaryBy(data=hm_all, Conc+log_Cc~Taxon+HM, FUN=mean)
#
lmm<- cast(spe_v[,-3], Taxon ~ HM, value=c("log_Cc.mean"))# HM from row to column
#
for (i in 2:8) {
  lmm[,i] <- as.numeric(scale(lmm[,i], scale=T))
}
#
melt.spe<- melt(lmm,variable_name="HM", id.vars=c("Taxon")) # HM from column to row
melt.spe<- rename(melt.spe, c(value="LC50"))
#
#
# Mixed Model
# Are there differences in ranked sensitivities?
#
summary(mmod<- glmer(LC50~ HM + (1|Taxon), data=melt.spe))
anova(mmod)
#
## Model diagnostics
op <- par(mfrow = c(1, 2))
plot(fitted(mmod), resid(mmod), xlab="Fitted", ylab="Residual", main="a") #fix factor
qqnorm(resid(mmod), main="b") #fix factor
qqnorm(ranef(mmod)$Taxon[[1]]) #random factor
#
# No differences - take median across different metals
# 
avg.HM.species<- transform(lmm, HM=apply(lmm, MARGIN= 1, FUN= median, na.rm=T))#standardized data + median HM
#
#
##############################################################################
#              Pairs-Pearson correlation between metals
##############################################################################
#
#
png(file='Pair_Plot.png', 
    width = 7, height = 7, units = "in", res = 300)  
par(mfrow=c(1,1), mar=c(4,4,1,1))
{
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x,breaks=seq(-3.5,3.5,by=1), plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
  }
  
  panel.cor<- function(x, y, method="pearson", digits=2,...) 
  {
    points(x,y,type="n");
    usr<- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1));
    correl<- cor.test(x,y, method=method, use="pairwise.complete.obs");
    sum.row<- rowSums(cbind(x,y));
    count<- length(sum.row[!is.na(sum.row)]);
    n=count;
    r=correl$estimate;
    color="black";
    fontface=2;
    
    txt<- format(r,digits=2)
    txt<- paste("r=", txt, "\nn=", n, sep="")
    text(0.5, 0.5, txt,col=color, cex=1.7)
  } 
}

pairs(avg.HM.species[,(2:8)], diag.panel=panel.hist, cex.labels = 2, font.labels=2,
      upper.panel=panel.cor, xlim = c(-3,3), ylim= c(-3,3)) 
#
dev.off()
#
# End
#
#############################################################################




