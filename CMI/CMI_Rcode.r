########################################################
##### CMI Data Analysis and Strategy Working group #####
#####                                              #####
##### Inferential model for CMI data analsyis      #####
##### NOV2022                                      #####
########################################################

### Load packages ---

#library(lme4)
library(msm)  # for delta method

### Import data ---

dat = read.csv('data/dat_test.csv')

      
### Fit the model ---

#mod = lmer(y~gr+x+y0+(1|subjid),data = dat) # with r.e. if 2 visits after V0. For > 2 visits the model needs to be changed if unstructured errors variane-covariace matrix needs to be accomodated
mod = lm(y~gr+y0+x,data = dat)
summary(mod)
vcov.mod = vcov(mod) # get var-cov matrix of effects
mean.mod = coef(mod) # get coefficients
names(mean.mod) = c('x1','x2','x3','x4')
rownames(vcov.mod) = c('x1','x2','x3','x4')
colnames(vcov.mod) = c('x1','x2','x3','x4')


### Compute estimate and SE via Delta method ---
m.x = mean(dat$x)
m.y0 = mean(dat$y0)
get.est = function(expr.name,x1,x2,x3,x4,m.y0,m.x){eval(parse(text=sub(".*~ ", "", print(expr.name))[2]))}
g.gmvacc=formula(~(10^(x1+x2+x3*m.y0+x4*m.x)-1)*(10^m.x))
g.gmvacc.est = get.est(expr.name = g.gmvacc,mean.mod[1],mean.mod[2],mean.mod[3],mean.mod[4],m.y0,m.x)
g.gmvacc.se = deltamethod(g=g.gmvacc, mean=mean.mod, cov=vcov.mod)
g.gmctr=formula(~(10^(x1+x3*m.y0+x4*m.x)-1)*(10^m.x))
g.gmctr.est = get.est(expr.name = g.gmctr,mean.mod[1],mean.mod[2],mean.mod[3],mean.mod[4],m.y0,m.x)
g.gmctr.se = deltamethod(g=g.gmvacc, mean=mean.mod, cov=vcov.mod)
g.eff=formula(~log10(((10^(x1+x2+x3*m.y0+x4*m.x)-1)*(10^m.x))/((10^(x1+x3*m.y0+x4*m.x)-1)*(10^m.x))))
g.eff.est = ifelse(!is.na(get.est(expr.name = g.eff,mean.mod[1],mean.mod[2],mean.mod[3],mean.mod[4],m.y0,m.x)),
                   get.est(expr.name = g.eff,mean.mod[1],mean.mod[2],mean.mod[3],mean.mod[4],m.y0,m.x),0.01)
g.eff.se = deltamethod(g=g.eff, mean=mean.mod, cov=vcov.mod)

### Get degrees of freedom (change to KR approximation in case of repeated measures) and compute CIs ---
df.mod = summary(mod)$df[2] # from lm (only one visit)

alpha = 0.05

ci.gmvacc = g.gmvacc.est+c(0,qt(alpha/2,df.mod),qt(1-alpha/2,df.mod))*g.gmvacc.se
ci.gmctr = g.gmctr.est+c(0,qt(alpha/2,df.mod),qt(1-alpha/2,df.mod))*g.gmctr.se
ci.eff = 10^(g.eff.est+c(0,qt(alpha/2,df.mod),qt(1-alpha/2,df.mod))*g.eff.se)
t.value = g.eff.est/g.eff.se
p.value = 1-pt(g.eff.est/g.eff.se,df.mod)
res = array(c(ci.gmvacc,ci.gmctr,ci.eff,t.value,p.value),dim = c(1,11))
colnames(res) = c('GMF_Vac','GMF_Vac_LL','GMF_Vac_UL',
                  'GMF_Ctr','GMF_Ctr_LL','GMF_Ctr_UL',
                  'GMF_Ratio','GMF_Ratio_LL','GMF_Ratio_UL',
                  't-value','p-value')
res
