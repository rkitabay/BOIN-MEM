
library(R2jags)
library(BOIN)

##-----------------Trial Setting--------------------##

pt.target <- 0.30 # Target toxicity probability
pe.null <- 0.25 # Lowest acceptable efficacy probability
assoc <- 0.20 # Association parameter of a Gumbel model
ctox <- 0.95 # Cutoff value of toxicity in equation (14)
ceff <- 0.90 # Cutoff value of efficacy in equation (15)

u00 <- 0.30 # Utility score at no toxicity and no efficacy
u11 <- 0.50 # Utility score at toxicity and efficacy

Ndose <- 5 # Number of dose levels
Np1 <- 30 # Maximum total sample size in Stage 1
s1 <- 12 # Maximum sample size treated at one of dose levels in Stage 1
Np2 <- 24 # Maximum total sample size in Stage 2
s2 <- 24 # Maximum sample size treated at one of dose levels in Stage 2
chsize <- 3 # Cohort size
sint <- 12 # Maximum sample size treated at one of dose levels for re-selection of R

tie <- 0.10 # Significance level in 3+3+CE

# Decision table and boundaries of BOIN
BOIN.design <- get.boundary(target = pt.target, ncohort = Np1, cohortsize = chsize)
BOIN.tab <- BOIN.design$boundary_tab
BOIN.ldes <- BOIN.design$lambda_d

# True toxicity and efficacy in this trial
p1tox <- c(0.06,0.10,0.24,0.40,0.41)
p1eff <- c(0.12,0.14,0.35,0.50,0.60)
p2tox <- c(0.23,0.26,0.28,0.46,0.50)
p2eff <- c(0.41,0.58,0.67,0.69,0.70) # Scenario 11 in the main manuscript

true.et <- list(p1tox = p1tox, p1eff = p1eff, p2tox = p2tox, p2eff = p2eff)
cat_utility(Ndose,true.et,u11,u00)


##-----------------Parameter Specification--------------------##
#---------For BOIN-MEM-----------#

a <- 1 # alpha of a quasi-beta prior distribution and a beta prior distribution to evaluate the admissible dose set
b <- a
w <- 0 # Weight parameter of Section 2.4
gamma <- 0  # Tuning parmaeter of Section 2.4
Narm <- 2 # Number of dose levels in Stage 2 (J_R)

#---------For 3+3+CE-----------#
a <- 1 # alpha of a beta prior distribution to evaluate the admissible dose set
b <- a

##-----------------Analysis--------------------##
#---------For BOIN-MEM-----------#


monitor.init<-list()
for(j in 1:Ndose){
  monitor.init[[j]]<-list(p1=data.frame(n=0,xt=0,xe=0),p2=data.frame(n=0,xt=0,xe=0))
}

# Conduct the BOIN design in Stage 1

p1.res<-BOIN(Ndose, Np1, s1, chsize, assoc, BOIN.tab, true.et, monitor=monitor.init)

act<-p1.res$act[1:Ndose] # Safe dose levels that were not eliminated
left.ch<-Np1/chsize-p1.res$last.ch

nx1<-nx(monitor=p1.res$monitor)
mtd<-select.mtd(target = pt.target, npts = nx1$p1n, ntox = nx1$p1xt)$MTD

if(!act[mtd]&&mtd!=99)stop("MTD is an eliminated dose")
if(mtd==99) act<-rep(F,Ndose)

monitor <- p1.res$monitor

# Select the reccomended Stage 2 dose set
p2dose<-select_dose_BM(Ndose,pt.target,pe.null,ctox,ceff,Narm,BOIN.ldes,act,mtd,monitor.ce=monitor,a,b)

# Steps 1-2 in Stage 2
ce<-enroll_Stage2_BM(chsize,pt.target,pe.null,ctox,ceff,s2,assoc,p2dose,p2ch=sint/chsize,true.et,monitor.ce=monitor,a,b)
monitor<-ce$monitor

# Step 3
p2dose<-select_dose_BM(Ndose,pt.target,pe.null,ctox,ceff,Narm,BOIN.ldes,act,mtd,monitor.ce=monitor,a,b)

# Step 4
ce<-enroll_Stage2_BM(chsize,pt.target,pe.null,ctox,ceff,s2,assoc,p2dose,p2ch=Np2/chsize+left.ch-sum(sint)/chsize,true.et,monitor.ce=monitor,a,b)

stop<-as.numeric(sum(!ce$act2.dose)>=1)

nx2<-nx(monitor=ce$monitor)

# Select tge OBD using MEM (step 5)
if(sum(ce$act2.dose)>1){
  p2act<-ce$p2dose[ce$act2.dose]
  
  post<-list()
  for(j in 1:length(p2act)){
    currj<-p2act[j]
    
    # Define exchangeability models
    if(j==length(p2act)){
      upex<-oneside_ex(Nexdose=Ndose-currj,lowside=F)
    }else upex<-oneside_ex(Nexdose=ceiling((p2act[j+1]-currj-1)/2),lowside=F)
    
    if(j==1){
      lwex<-oneside_ex(Nexdose=currj-1,lowside=T)
    }else lwex<-oneside_ex(Nexdose=ceiling((currj-p2act[j-1]-1)/2),lowside=T)
    
    if(all(is.nan(lwex))){
      Nlwex<-0
    }else Nlwex<-ncol(lwex)
    if(all(is.nan(upex))){
      Nupex<-0
    }else Nupex<-ncol(upex)
    comp<-c((currj-1)-Nlwex,Nlwex,1,Nupex,(Ndose-currj)-Nupex)
    
    if(comp[2]==0&&comp[4]==0){
      omega.234mtx<-matrix(1,ncol=1,nrow=1)
    }else if(comp[4]==0){
      omega.234mtx<-cbind(lwex,1)
    }else if(comp[2]==0){
      omega.234mtx<-cbind(1,upex)
    }else{
      memk<-1
      for(memj in 1:nrow(lwex)){
        for(memi in 1:nrow(upex)){
          if(memk==1)omega.234mtx<-c(lwex[memj,],1,upex[memi,])
          else omega.234mtx<-rbind(omega.234mtx,c(lwex[memj,],1,upex[memi,]))
          memk<-memk+1
        }
      }
    }
    
    omega.all<-omega.234mtx
    if(comp[1]>0) omega.all<-cbind(matrix(0,ncol=comp[1],nrow=nrow(omega.234mtx)),omega.all)
    if(comp[5]>0) omega.all<-cbind(omega.all,matrix(0,ncol=comp[5],nrow=nrow(omega.234mtx)))
    omega.all<-rbind(0,omega.all)
    
    # Estimate utility borrowed information from other dose levels and stages
    sum.p1nu<-omega.all%*%as.matrix(nx2[,c(2,5)])
    sum.p1nu<-matrix(sum.p1nu,ncol=2)
    colnames(sum.p1nu)<-c("n","u")
    
    post.m<-c()
    for(memi in 1:nrow(omega.all)){
      prod.beta<-prod(((beta(a+nx2$p1u,b+(nx2$p1n-nx2$p1u)))/beta(a,b))^(1-omega.all[memi,]))
      post.m[memi]<-(beta(a+sum.p1nu[memi,2]+nx2$p2u[currj],b+(sum.p1nu[memi,1]-sum.p1nu[memi,2])+(nx2$p2n[currj]-nx2$p2u[currj]))/beta(a,b))*prod.beta
    }
    
    if(w==9)pre.omega<-rep(1/(nrow(omega.all)),nrow(omega.all))
    else pre.omega<-c(w,rep((1-w)/(nrow(omega.all)-1),(nrow(omega.all)-1)))
    
    tmp.mtx<-(omega.all%*%matrix(1,nrow=ncol(omega.all),ncol=1))^gamma
    pre.omega[-1]<-tmp.mtx[-1]/sum(tmp.mtx[-1])*sum(pre.omega[-1])
    
    post.omega<-post.m*pre.omega/as.numeric(t(post.m)%*%pre.omega)
    
    post.a<-a+sum.p1nu[,2]+nx2$p2u[currj]
    post.b<-b+(sum.p1nu[,1]-sum.p1nu[,2])+(nx2$p2n[currj]-nx2$p2u[currj])
    
    post[[currj]]<-list(omega=post.omega,a=post.a,b=post.b)
  }
  
  p2indc<-rep(F,Ndose)
  p2indc[p2act]<-T
  obd<-maxu_MEM(indc=p2indc,post=post)
  
}else if(sum(ce$act2.dose)==1){
  obd<-ce$p2dose[ce$act2.dose]
}else if(sum(ce$act2.dose)==0){
  obd<-0
}

obd


#---------For 3+3+CE to select OBD using information from Stage 1 -----------#


monitor.init<-list()
for(j in 1:Ndose){
  monitor.init[[j]]<-list(p1=data.frame(n=0,xt=0,xe=0),p2=data.frame(n=0,xt=0,xe=0))
}

# Conduct the 3+3 design in Stage 1
p1.res<-ThreePlusThree(Ndose,chsize,Np1,assoc,s1,true.et,monitor=monitor.init)

act<-p1.res$act[1:Ndose] # Safe dose levels that were not eliminated
left.ch<-Np1/chsize-p1.res$last.ch

nx1<-nx(monitor=p1.res$monitor)
mtd<-select.mtd(target = pt.target, npts = nx1$p1n, ntox = nx1$p1xt)$MTD

if(!act[mtd]&&mtd!=99)stop("MTD is an eliminated dose")
if(mtd==99) act<-rep(F,Ndose)

# Enroll patients in Stage 2
p2dose<-select_dose_33CE(act,mtd,monitor.ce=p1.res$monitor,a,b)
ce<-enroll_Stage2_33CE(Np2,s2,assoc,p2dose,p2ch=Np2/chsize+left.ch,true.et,monitor.ce=p1.res$monitor,a,b)

stop<-as.numeric(sum(!ce$act2.dose)>=1)

nx2<-nx(monitor=ce$monitor)
nv<-nx2$p1n+nx2$p2n
xev<-nx2$p1xe+nx2$p2xe

# Select the OBD using the exact binomial test for the efficacy
if(sum(ce$act2.dose)>1){
  p2act<-ce$p2dose[ce$act2.dose]
  
  signf<-c()
  for(j in 1:length(p2act)){
    currj<-p2act[j]
    signf[j]<-binom.test(x=xev[currj],n=nv[currj],alternative="greater",p=pe.null)$p.value<tie
  }
  if(sum(signf)>0)obd<-p2act[which(signf==1)][1]
  else obd<-0
}else if(sum(ce$act2.dose)==1){
  obd<-ce$p2dose[ce$act2.dose]
}else if(sum(ce$act2.dose)==0){
  obd<-0
}

obd # 0 means that there is no OBD


#---------For 3+3+CE to select OBD NOT using information from Stage 1 -----------#


monitor.init<-list()
for(j in 1:Ndose){
  monitor.init[[j]]<-list(p1=data.frame(n=0,xt=0,xe=0),p2=data.frame(n=0,xt=0,xe=0))
}

p1.res<-ThreePlusThree(Ndose,chsize,Np1,assoc,s1,true.et,monitor=monitor.init)

act<-p1.res$act[1:Ndose] # Safe dose levels that were not eliminated
left.ch<-Np1/chsize-p1.res$last.ch

nx1<-nx(monitor=p1.res$monitor)
mtd<-select.mtd(target = pt.target, npts = nx1$p1n, ntox = nx1$p1xt)$MTD

if(!act[mtd]&&mtd!=99)stop("MTD is an eliminated dose")
if(mtd==99) act<-rep(F,Ndose)


p2dose<-select_dose_33CE(act,mtd,monitor.ce=p1.res$monitor,a,b)
ce<-enroll_Stage2_33CE(Np2,s2,assoc,p2dose,p2ch=Np2/chsize+left.ch,true.et,monitor.ce=p1.res$monitor,a,b)

stop<-as.numeric(sum(!ce$act2.dose)>=1)

nx2<-nx(monitor=ce$monitor)
nv<-nx2$p2n
xev<-nx2$p2xe

if(sum(ce$act2.dose)>1){
  p2act<-ce$p2dose[ce$act2.dose]
  
  signf<-c()
  for(j in 1:length(p2act)){
    currj<-p2act[j]
    signf[j]<-binom.test(x=xev[currj],n=nv[currj],alternative="greater",p=pe.null)$p.value<tie
  }
  if(sum(signf)>0)obd<-p2act[which(signf==1)][1]
  else obd<-0
}else if(sum(ce$act2.dose)==1){
  obd<-ce$p2dose[ce$act2.dose]
}else if(sum(ce$act2.dose)==0){
  obd<-0
}

obd # 0 means that there is no OBD
