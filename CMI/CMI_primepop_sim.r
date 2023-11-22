#########################################################
#################### CMI Working Group ##################
#################### Simulations  ##################
#########################################################


rm(list=ls())

###### Loading packages and settign wd ----
library(car)
library(lme4)
library(msm)
library(ggplot2)
library(dplyr)
library(officer)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)
library(vctrs)



###### Setting up parameters and scenarios ----

n.vec = c(10,50,100,200) # number of subjects per group
ngr = 2 # number of groups
mu.ind.v0 = log10(350)  # V0 Freq_ind
mu.bkg = log10(c(350,350))  # V0, V2 Freq_bkg
mu.ind.v2 = log10(array(c(1300,1300,
                          1300,1500,
                          1300,2000,
                          1300,3000,
                          1300,3500,
                          1300,4000,
                          1300,5000,
                          2500,3500,
                          2500,5000),dim=c(2,9)))  # V2 Freq_ind, for different scenarios
sd.e.vec = c(0.3,0.5,0.7) # error SD
z.corr.vec = c(0.2,0.5,1)   # correlation for Freq_ind V0/V2
## MC = 10   # number of simulations
alpha = 0.05   # significance level


###### Running simulations ----

parallel_sim <- function() {
  res = data.frame(matrix(0,ncol=11))
  colnames(res)<- c("mean","mean1","n","sd.e", "z.corr","mod.val","mod.ll","mod.ul","mod1.val","mod1.ll","mod1.ul")
  na = data.frame(matrix(0, ncol=6))
  colnames(na)<- c("mean","mean1","n","sd.e","z.corr","is_na")
  pv = data.frame(matrix(0, ncol=8))
  colnames(pv)<- c("mean","mean1","n","sd.e","z.corr","pv.mod","pv.mod1","pv.wil")
  impvals.sp = data.frame(matrix(0,ncol=7))
  colnames(impvals.sp) <- c("mean","mean1","n","sd.e","z.corr","group","impvals.sp")
  impvals.sp0 = data.frame(matrix(0,ncol=6))
  colnames(impvals.sp0) <- c("mean","mean1","n","sd.e","z.corr","impvals.sp0")
  corr.prepost = data.frame(matrix(0, ncol=6))
  colnames(corr.prepost) <- c("mean","mean1","n","sd.e","z.corr","corr.prepost")
  
  for(j in 1:dim(mu.ind.v2)[2]){
    for(i.n in 1:length(n.vec)){
      for(i.sd in 1:length(sd.e.vec)){
        for (i.corr in 1:length(z.corr.vec)) {
          # set.seed(j*i.n*i.sd)
          ### Generate the data
          x = NULL
          y0 = NULL
          y = NULL
          sd.subj = sqrt(sd.e.vec[i.sd]^2*z.corr.vec[i.corr])  # to induce correlation V0/V2, (will move to more general unstructured Sigma later on, if LMM simulations will be added)
          re.subj = rnorm(2*n.vec[i.n],0,sd.subj)
          visit = factor(rep(c('V0','V2'),each = n.vec[i.n]*ngr))
          group = factor(rep(rep(c('ctr','vac'),each = n.vec[i.n]),ngr))
          subjid = factor(rep(c(1:(2*n.vec[i.n])), ngr))
          for(i.subj in 1:length(levels(subjid))){
            for(i.vis in 1:length(levels(visit))){
              
              x[subjid == levels(subjid)[i.subj] &
                  visit == levels(visit)[i.vis]] = rnorm(1,mu.bkg[i.vis],sd.e.vec[i.sd])
              if(i.vis == 1){
                
                y0[subjid == levels(subjid)[i.subj] &
                     visit == levels(visit)[i.vis]] = rnorm(1,mu.ind.v0,sd.e.vec[i.sd])+re.subj[i.subj]}
              
              y[subjid == levels(subjid)[i.subj] &
                  visit == levels(visit)[i.vis] &
                  group == 'ctr'] = rnorm(1,mu.ind.v2[1,j],sd.e.vec[i.sd])+re.subj[i.subj]
              
              y[subjid == levels(subjid)[i.subj] &
                  visit == levels(visit)[i.vis] &
                  group == 'vac'] = rnorm(1,mu.ind.v2[2,j],sd.e.vec[i.sd])+re.subj[i.subj]
            }
          }
          
          dat = data.frame('subjid' = subjid[visit == 'V2'],
                           'group' = group[visit == 'V2'],
                           'visit' = visit[visit == 'V2'],
                           'y' = y[visit == 'V2'],
                           'y0' = y0[visit == 'V0'],
                           'x' = x[visit == 'V2'],
                           'x0' = x[visit == 'V0'])
          ##  dat1 = rbind(dat1, data.frame("MC"=mc,"mean"=mu.ind.v2[2,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],dat))
          
          #corr.prepost[mc,j,i.n,i.sd] = cor(dat[,c(4,5)])[1,2]
          xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"corr.prepost"=cor(dat[,c(4,5)])[1,2])
          
          corr.prepost =rbind(corr.prepost, xyz)
          
          ### Fit the models
          ## y = induction freq
          mod = lm(y~group+y0+x,data = dat) # dropping r.e. for only one visit
          summary(mod)
          vcov.mod = vcov(mod) # get var-cov matrix of effects
          mean.mod = coef(mod) # get coefficients
          names(mean.mod) = c('x1','x2','x3','x4')
          rownames(vcov.mod) = c('x1','x2','x3','x4')
          colnames(vcov.mod) = c('x1','x2','x3','x4')
          ## y = spec freq
          dat$sp = log10(10^dat$y-10^dat$x)
          dd <- dat
          dd$na <- 1
          dd$na[is.na(dd$sp)]=0
          if(sum(is.na(dd$sp)) > 0 ){
            
            
            prop <- dd %>% group_by(group, na) %>% summarise(Freq = n())
            
            prop <- data.frame(prop)
            
            prop$prop = prop$Freq/n.vec[i.n]
            
            prop <- prop[which(prop$na==0),]
            
            if (nrow(prop) == 2) {
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "control", "impvals.sp"= prop[prop$group =="ctr", 4])
              impvals.sp = rbind(impvals.sp,xyz)
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "Vaccine", "impvals.sp"= prop[prop$group =="vac", 4])
              impvals.sp = rbind(impvals.sp,xyz)
            } else if (prop$group[1]=="vac") {
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "control", "impvals.sp"= 0)
              impvals.sp = rbind(impvals.sp,xyz)
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "Vaccine", "impvals.sp"= prop[prop$group =="vac", 4])
              impvals.sp = rbind(impvals.sp,xyz)
            } else {
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "control", "impvals.sp"= prop[prop$group =="ctr", 4])
              impvals.sp = rbind(impvals.sp,xyz)
              xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "Vaccine", "impvals.sp"= 0)
              impvals.sp = rbind(impvals.sp,xyz)
            }
          } else {
            xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "control", "impvals.sp"= 0)
            impvals.sp = rbind(impvals.sp,xyz)
            xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr], "group" = "Vaccine", "impvals.sp"= 0)
            impvals.sp = rbind(impvals.sp,xyz)
          }  
          
          dat$sp[is.na(dat$sp)] = 0
          dat$sp0 = log10(10^dat$y0-10^dat$x0)
          xyz <- data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"impvals.sp0"=sum(is.na(dat$sp0))/length(dat$sp0))
          impvals.sp0 = rbind(impvals.sp0,xyz)
          #impvals.sp0[mc,j,i.n,i.sd] = sum(is.na(dat$sp0))/length(dat$sp0)
          dat$sp0[is.na(dat$sp0)] = 0
          mod1 = lm(sp~group+sp0,data = dat)
          ## y = spec freq on original scale for non-parametric test
          dat$sp.orig = 10^dat$y-10^dat$x
          
          ### Compute estimate and SE caculation via delta method for mod
          m.x = mean(dat$x)
          m.y0 = mean(dat$y0)
          
          get.est = function(expr.name,x1,x2,x3,x4,m.y0,m.x){eval(parse(text=sub(".*~ ", "", expr.name)[2]))}
          form=sprintf("~ log10(((10^(x1+x2+x3*%f+x4*%f-%f)-1))/((10^(x1+x3*%f+x4*%f-%f)-1)))",m.y0 ,m.x,m.x,m.y0 ,m.x,m.x)
          # g.eff=formula(~log10(((10^(x1+x2+x3*m.y0+x4*m.x-m.x)-1))/((10^(x1+x3*m.y0+x4*m.x-m.x)-1))))
          g.eff=as.formula(form)
          g.eff.est = get.est(expr.name = g.eff,mean.mod[1],mean.mod[2],mean.mod[3],mean.mod[4],m.y0,m.x)
          #print(c(g.eff.est,m.x,m.y0,mean.mod))
          x1.mc = mean.mod[1]
          x2.mc = mean.mod[2]
          x3.mc = mean.mod[3]
          x4.mc = mean.mod[4]
          
          xx= ifelse (is.na(g.eff.est),1,0)
          xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"is_na"=xx)
          na = rbind(na,xyz) 
          
          ###  na[which(na$is_na==1),]
          
          if (xx==1) {   
            g.eff.est =   
              ifelse((10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)<0,
                     log10(1/((10^(x1.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)/abs(10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)+1)),
                     log10(((10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)/abs(10^(x1.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)+1)))
          }
          
          if (xx==0 & (10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)<0 & (10^(x1.mc +x3.mc*m.y0+x4.mc*m.x-m.x)-1)<0) {   ## FSO :- Here we consider the condition when both D & N will be negative (third condition is not strictly needed but I left it in)
            g.eff.est = log10(((10^(x1.mc +x3.mc*m.y0+x4.mc*m.x-m.x)-1))/((10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1))) 
            g.eff=as.formula( sprintf("~log10(((10^(x1+x3*%f+x4*%f-%f)-1))/((10^(x1+x2+x3*%f+x4*%f-%f)-1)))",m.y0,m.x,m.x,m.y0,m.x,m.x))
            #g.eff = formula(~log10(((10^(x1+x3*m.y0+x4*m.x-m.x)-1))/((10^(x1+x2+x3*m.y0+x4*m.x-m.x)-1))))
          }# Here in some case it may happen both treatment are worsening but then we're are not able to figure out 
          
          
          if (xx==1){ 
            if ((10^(x1.mc+x2.mc+x3.mc*m.y0+x4.mc*m.x-m.x)-1)<0){
              g.eff=as.formula(sprintf("~log10(1/((10^(x1+x3*%f+x4*%f-%f)-1)/sqrt((10^(x1+x2+x3*%f+x4*%f-%f)-1)+1)^2))",m.y0,m.x,m.x,m.y0,m.x,m.x))
              #g.eff = formula(~log10(1/((10^(x1+x3*m.y0+x4*m.x-m.x)-1)/sqrt((10^(x1+x2+x3*m.y0+x4*m.x-m.x)-1)+1)^2))) 
            } else {
              g.eff=as.formula(sprintf("~log10(((10^(x1+x2+x3*%f+x4*%f-%f)-1)/sqrt((10^(x1+x3*%f+x4*%f-%f)-1)+1)^2))",m.y0,m.x,m.x,m.y0,m.x,m.x))
              #g.eff = formula(~log10(((10^(x1+x2+x3*m.y0+x4*m.x-m.x)-1)/sqrt((10^(x1+x3*m.y0+x4*m.x-m.x)-1)+1)^2))) 
            }
          }
          
          g.eff.se = deltamethod(g=g.eff, mean=mean.mod, cov=vcov.mod)
          #g.eff.se=2
          ## Get degrees of freedom (KR approx used in RSV 002 SAS code for LMM) - only different for LMM, not implemented now
          df.mod = summary(mod)$df[2] # for lm (only one visit)
          df.mod1 = summary(mod1)$df[2]
          mod1.eff = summary(mod1)$coefficients[2,1:2]
          
          
          ### Compute CIs and p-values:
          ci.eff.mod = 10^(g.eff.est+c(0,qt(alpha/2,df.mod),qt(1-alpha/2,df.mod))*g.eff.se)
          ci.eff.mod1 = 10^(mod1.eff[1]+c(0,qt(alpha/2,df.mod1),qt(1-alpha/2,df.mod1))*mod1.eff[2])
          
          
          xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"mod.val"=ci.eff.mod[1],"mod.ll"=ci.eff.mod[2],"mod.ul"=ci.eff.mod[3],"mod1.val"=ci.eff.mod1[1],"mod1.ll"=ci.eff.mod1[2],"mod1.ul"=ci.eff.mod1[3])
          res = rbind(res,xyz)
          pv.mod = 1-pt(g.eff.est/g.eff.se,df.mod)
          pv.mod1 = 1-pt(mod1.eff[1]/mod1.eff[2],df.mod1)
          pv.wil = wilcox.test(dat$sp.orig[dat$group == 'ctr'], dat$sp.orig[dat$group == 'vac'],alternative = 'less')$p.value
          
          xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"=mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"pv.mod"=pv.mod,"pv.mod1"=pv.mod1,"pv.wil"=pv.wil)
          pv = rbind(pv,xyz)
          
          list_a<- list(corr.prepost,impvals.sp ,impvals.sp0 , na, res , pv )
        }
        #print(paste('sample size',i.n,'of',length(n.vec)))
      }
      #print(paste('mean scenario',j,'of',dim(mu.ind.v2)[2]))
    }
  }
  return(list_a)
}

#corr.prepost <- data.frame("MC" =1 ,z[[1]] )
#### Parallel programming  :- 
library(foreach)
library(iterators)
library(parallel)
library(doParallel)

library(ranger)

res = data.frame(matrix(0,ncol=12))
colnames(res)<- c("MC","mean","mean1","n","sd.e","z.corr","mod.val","mod.ll","mod.ul","mod1.val","mod1.ll","mod1.ul")
na = data.frame(matrix(0, ncol=7))
colnames(na)<- c("MC","mean","mean1","n","sd.e","z.corr","is_na")
pv = data.frame(matrix(0, ncol=9))
colnames(pv)<- c("MC","mean","mean1","n","sd.e","z.corr","pv.mod","pv.mod1","pv.wil")
impvals.sp = data.frame(matrix(0,ncol=8))
colnames(impvals.sp) <- c("MC","mean","mean1","n","sd.e","z.corr","group","impvals.sp")
impvals.sp0 = data.frame(matrix(0,ncol=7))
colnames(impvals.sp0) <- c("MC","mean","mean1","n","sd.e","z.corr","impvals.sp0")
corr.prepost = data.frame(matrix(0, ncol=7))
colnames(corr.prepost) <- c("MC","mean","mean1","n","sd.e","z.corr","corr.prepost")



n_cores <- 40


registerDoParallel(cores = n_cores)

res2 <- foreach(i = 1:1000, .combine = "rbind") %dopar% {
  set.seed(i)
  parallel_sim() 
}


for(i in 1:1000){
  x <- data.frame("MC" =i ,res2[i,1])
  colnames(x) <- names(corr.prepost)
  corr.prepost <- rbind(corr.prepost,x[-1,] )
  
  x <- data.frame("MC" =i ,res2[i,2])
  colnames(x) <- names(impvals.sp)
  impvals.sp <- rbind(impvals.sp,x[-1,] )
  
  x <- data.frame("MC" =i ,res2[i,3])
  colnames(x) <- names(impvals.sp0)
  impvals.sp0 <- rbind(impvals.sp0,x[-1,] )
  
  x <- data.frame("MC" =i ,res2[i,4])
  colnames(x) <- names(na)
  na <- rbind(na,x[-1,] )
  
  x <- data.frame("MC" =i ,res2[i,5])
  colnames(x) <- names(res)
  res <- rbind(res,x[-1,] )
  
  x <- data.frame("MC" =i ,res2[i,6])
  colnames(x) <- names(pv)
  pv <- rbind(pv,x[-1,] )
  
  
}

###### Summarizing and displaying results ----

#### . Creating power and coverage probability arrays ----

est.res = data.frame(matrix(0,ncol=11))
colnames(est.res)<- c("mean","mean1","n","sd.e","z.corr","mod.val","mod.ll","mod.ul","mod1.val","mod1.ll","mod1.ul")

covpr.res = data.frame(matrix(0,ncol=7))
colnames(covpr.res) <- c("mean","mean1","n","sd.e","z.corr","mod","mod1")


pow.res = data.frame (matrix(0,ncol = 8))
colnames(pow.res) <- c("mean","mean1","n","sd.e","z.corr","pv.mod","pv.mod1","pv.will")

impprop =data.frame(matrix(0,ncol=7))
colnames(impprop) <- c("mean","mean1","n","sd.e","z.corr","group","impprop")

impprop0 =data.frame(matrix(0,ncol=6))
colnames(impprop0) <- c("mean","mean1","n","sd.e","z.corr","impprop0")

corr.prepost.av = data.frame(matrix(0,ncol=6))
colnames(corr.prepost.av) <- c("mean","mean1","n","sd.e","z.corr","corr.prepost.av")
for(j in 1:dim(mu.ind.v2)[2]){
  for(i.n in 1:length(n.vec)){
    for(i.sd in 1:length(sd.e.vec)){
      for(i.corr in 1:length(z.corr.vec)){
        res1 = res[which(res$mean== mu.ind.v2[2,j] & res$mean1== mu.ind.v2[1,j] & res$n==n.vec[i.n] & res$sd.e==sd.e.vec[i.sd] & res$z.corr==z.corr.vec[i.corr]),]
        a1 = 10^apply(log10(res1[,7:12]),2,mean) 
        #a1 = apply(res1[,5:10],2,mean)
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"mod.val"=a1[1],"mod.ll"=a1[2],"mod.ul"=a1[3],"mod1.val"=a1[4],"mod1.ll"=a1[5],"mod1.ul"=a1[6])
        est.res = rbind(est.res,xyz)
        
        true.eff.mod = c((10^mu.ind.v2[2,j]-1)/(10^mu.ind.v2[1,j]-1))
        true.eff.mod1 = c((10^mu.ind.v2[2,j]-10^mu.ind.v0)/(10^mu.ind.v2[1,j]-10^mu.ind.v0))
        
        
        covpr.res.mod = mean(as.numeric(apply(res1[,8:9],1,function(.x){.x[1] < true.eff.mod & .x[2] > true.eff.mod})))
        covpr.res.mod1 = mean(as.numeric(apply(res1[,11:12],1,function(.x){.x[1] < true.eff.mod1 & .x[2] > true.eff.mod1})))
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"mod"=covpr.res.mod,"mod1"=covpr.res.mod1)
        covpr.res = rbind(covpr.res, xyz)
        
        pv1 = pv[which(res$mean== mu.ind.v2[2,j] & res$mean1== mu.ind.v2[1,j] & res$n==n.vec[i.n] & res$sd.e==sd.e.vec[i.sd] & res$z.corr==z.corr.vec[i.corr]),]
        a2 = apply(cbind(apply(pv1[,7:9],2,function(.x){.x<alpha/2})),2,mean)
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"pv.mod"=a2[1],"pv.mod1"=a2[2],"pv.will"=a2[3])
        pow.res = rbind(pow.res, xyz)
        
        impprop.1 = impvals.sp[which( impvals.sp$mean== mu.ind.v2[2,j] &  impvals.sp$mean1== mu.ind.v2[1,j] &  impvals.sp$n==n.vec[i.n] &  impvals.sp$sd.e==sd.e.vec[i.sd] &  impvals.sp$group=="control" & impvals.sp$z.corr==z.corr.vec[i.corr]) ,]
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"group"="control","impprop"=mean(impprop.1$impvals.sp))
        impprop = rbind(impprop,xyz)
        
        impprop.1 = impvals.sp[which( impvals.sp$mean== mu.ind.v2[2,j] &  impvals.sp$mean1== mu.ind.v2[1,j] &  impvals.sp$n==n.vec[i.n] &  impvals.sp$sd.e==sd.e.vec[i.sd] &  impvals.sp$group=="Vaccine" & impvals.sp$z.corr==z.corr.vec[i.corr]) ,]
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"group"="Vaccine","impprop"=mean(impprop.1$impvals.sp))
        impprop = rbind(impprop,xyz)
        
        impprop0.1 = impvals.sp0[which(res$mean== mu.ind.v2[2,j]  & res$mean1== mu.ind.v2[1,j] & res$n==n.vec[i.n] & res$sd.e==sd.e.vec[i.sd] & res$z.corr==z.corr.vec[i.corr]),]
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"impprop0"=mean(impprop0.1$impvals.sp0))
        impprop0 = rbind(impprop0,xyz)
        
        corr =corr.prepost[which(res$mean== mu.ind.v2[2,j]  & res$mean1== mu.ind.v2[1,j] & res$n==n.vec[i.n] & res$sd.e==sd.e.vec[i.sd] & res$z.corr==z.corr.vec[i.corr]),]
        xyz = data.frame("mean"=mu.ind.v2[2,j],"mean1"= mu.ind.v2[1,j],"n"=n.vec[i.n],"sd.e"=sd.e.vec[i.sd],"z.corr"=z.corr.vec[i.corr],"corr.prepost.av"=mean(corr$corr.prepost))
        corr.prepost.av = rbind(corr.prepost.av, xyz)
      }
    }
  }
}


## Checking NA in delta method:-
vv=na%>% group_by(mean,n,sd.e,is_na) %>% summarise(n1 = n()) %>% mutate(Freq = n1/sum(n1))
vv=vv[vv$is_na==0,]
vv$freq1=1-vv$Freq

table(vv$freq1)


## Checking proportion of imputed values
impprop
impprop0

impprop[-1,] %>% group_by(sd.e) %>% summarise(mean= mean(impprop))


## Checking correlation
corr.prepost.av[-1,] %>% group_by(sd.e) %>% summarise(mean = mean(corr.prepost.av), median = median(corr.prepost.av))


# ggplot :-


## Bias and CI width - mod
aa <- data.frame(t(as.data.frame(mu.ind.v2)))
aa$true = ((10^aa$X2)-10^mu.bkg[2])/((10^aa$X1)-10^mu.bkg[1])
#aa <- aa[,2:3]
colnames(aa) <- c("mean1","mean","true")



aa$scenario = c("scenario_1","scenario_2","scenario_3","scenario_4","scenario_5","scenario_6","scenario_7","scenario_8","scenario_9")

impprop <- merge(impprop, aa, by =c("mean", "mean1"))

est.res1 <- merge(est.res, aa, by =c("mean", "mean1"))
est.res1$n_pergroup_sd <- paste0(est.res1$n, "-", est.res1$sd.e)

est.res1$n <- as.factor(est.res1$n)

z.corr.spec <- c(0.17, 0.33, 0.5)

bias_mod <- function(data, Title1,subtitle1){
  
  ggplot(data=data, aes(x=true, y=mod.val, color=n_pergroup_sd)) +
    geom_line()+
    geom_point() + theme(legend.position = "right") + 
    labs(title=Title1,
         subtitle=subtitle1) +
    xlab ("True Vaccine Effect ratio (B/A)") + ylab ("Average point estimated from the model") +
    # geom_ribbon(aes(ymin = mod.ll, ymax = mod.ul), alpha = 0.3) +
    geom_abline()
}

doc <- read_docx()



for (i in 1:length(sd.e.vec)){
  for (j in 1:length(z.corr.vec)){
    p<- bias_mod(data=est.res1[which(est.res1$sd.e==sd.e.vec[i] & est.res1$z.corr==z.corr.vec[j]),], Title1 = paste0("GM ratio of fold increase from background - Vaccine/Control with "),subtitle1 = paste0("SD =",sd.e.vec[i], " & correlation=",z.corr.spec[j]) )
    doc <- body_add_gg(doc, value =p , style = "centered" )
  }
}

print(doc, target  = "graph_mod.docx")

## Bias and CI width - mod1
# aa <- data.frame(t(as.data.frame(mu.ind.v2)))   # FSO added lne, otherwise R cdoes not find the column 'true1'
# aa$true1 = ((10^aa$X2-10^mu.bkg[2])/(10^aa$X1-10^mu.bkg[1]))

# aa <- aa[,2:3]
# colnames(aa) <- c("mean","true1")

#est.res1 <- merge(est.res, aa, by =c("mean", "mean1"))
#est.res1$n_pergroup_sd <- paste0(est.res1$n, "-", est.res1$sd.e)


bias_mod1 <- function(data, Title1,subtitle1){
  
  ggplot(data=data, aes(x=true, y=mod1.val, color=n_pergroup_sd)) +
    geom_line()+
    geom_point() + theme(legend.position = "right") + 
    labs(title=Title1,
         subtitle=subtitle1) +
    xlab ("True Vaccine Effect ratio (B/A)") + ylab ("Average point estimated from the model") +
    # geom_ribbon(aes(ymin = mod.ll, ymax = mod.ul), alpha = 0.3) +
    geom_abline()
}

doc1 <- read_docx()



for (i in 1:length(sd.e.vec)){
  for (j in 1:length(z.corr.vec)){
    g<- bias_mod1(data=est.res1[which(est.res1$sd.e==sd.e.vec[i] & est.res1$z.corr==z.corr.vec[j]),], Title1 = paste0("GM ratio of fold increase from background - Vaccine/Control with "),subtitle1 = paste0("SD =",sd.e.vec[i], " & correlation=",z.corr.spec[j]) )
    doc1 <- body_add_gg(doc1, value =g , style = "centered" )
  }
}

print(doc1, target  = "graph_mod1.docx")


## CI coverage probability
covpr.res

## Tests power
z.corr.spec <- c(0.17, 0.33, 0.5)
pow.res = read_csv("pow.res.csv")


aa <- data.frame(t(as.data.frame(mu.ind.v2)))
aa$true = ((10^aa$X2)-10^mu.bkg[2])/((10^aa$X1)-10^mu.bkg[1])
#aa <- aa[,2:3]
colnames(aa) <- c("mean1","mean","true")

aa$Scenario = "1 to 7"
aa$Scenario[c(8,9)] = "8 and 9"

aa$mean1=round(aa$mean1,5)
aa$mean=round(aa$mean,5)
pow.res$mean=round(pow.res$mean,5)
pow.res$mean1=round(pow.res$mean1,5)
pow.res1 <- merge(pow.res, aa, by =c("mean", "mean1")) %>% rename(MOD=pv.mod, MOD1=pv.mod1,WILCOXON=pv.will)

pow.res2 <- pivot_longer(pow.res1, cols=6:8, names_to = "Method", values_to = "Value")

pow.res2$`n per group` = as.factor(pow.res2$n)

pow_graph <- function(data, Title1,subtitle1){
  
  ggplot(data=data, aes(x=true, y=Value, color=`n per group`,linetype =Method)) +
    geom_line()+
    scale_linetype_manual(values = c("solid","dotted","dashed"))+
    #scale_color_manual(values=c("red","green","blue","black"))
    geom_point(aes(shape=Scenario)) + theme(legend.position = "right") +  
    labs(title=Title1,
         subtitle=subtitle1) +
    xlab ("True vaccine effect ratio (B/A)") + ylab ("Rejection probability") +
    geom_abline(intercept =0.025,slope = 0, col='black')
}

doc3 <- read_docx()



for (i in 1:length(sd.e.vec)){
  for (j in 1:length(z.corr.vec)){
    p3<- pow_graph(data=pow.res2[which(pow.res2$sd.e==sd.e.vec[i] & pow.res2$z.corr==z.corr.vec[j]),],
                   Title1 = paste0("Test for vaccine effect ratio (B/A) - Primed population"),subtitle1 = paste0("SD =",sd.e.vec[i], " & correlation=",z.corr.spec[j]) )
    doc3 <- body_add_gg(doc3, value =p3 , style = "centered" )
  }
}

print(doc3, target  = "power_graph_primed.docx")



## output   ###
write_csv(covpr.res, "covpr.res.csv")
write_csv(pow.res,"pow.res.csv")
write_csv(impprop,"impprop.csv")
write_csv(impprop0,"impprop0.csv")
write_csv(corr.prepost.av,"corr.prepost.av.csv")











