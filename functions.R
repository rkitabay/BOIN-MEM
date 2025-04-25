
BOIN<-function(Ndose, Np1, s1, chsize, assoc, BOIN.tab, true.et, monitor){
  currj<-1
  act.dose<-c(rep(TRUE,Ndose),FALSE)
  cnt.dose<-rep(0,Ndose)
  for(j in 1:(Np1/chsize)){
    cnt.dose[currj]<-cnt.dose[currj]+1
    
    if(cnt.dose[currj]==1) monitor<-monit_update(stage=1,dose=currj,n=chsize,true.et=true.et,c=assoc,monitor=monitor,append=F)
    else monitor<-monit_update(stage=1,dose=currj,n=chsize,true.et=true.et,c=assoc,monitor=monitor,append=T)
    
    
    cjn<-sum(monitor[[currj]]$p1$n)
    cjxt<-sum(monitor[[currj]]$p1$xt)
    
    bound<-BOIN.tab[-1,cjn/chsize]
    if(cjxt<=bound[1]){
      if(act.dose[currj+1])currj<-currj+1
    }else if(bound[2]<=cjxt&&cjxt<bound[3]){
      if(currj!=1)currj<-currj-1
    }else if(cjxt>=bound[3]){
      act.dose[currj:Ndose]<-FALSE
      if(currj!=1)currj<-currj-1
      else break
    }
    
    if(cjn>=s1) break
  }
  list(monitor=monitor,act=act.dose,last.ch=j)
}


cat_utility<-function(sim){
  
  true.et<-create_et(sim$Ndose,sim$scen)
  
  true1.util<-c()
  true2.util<-c()
  for(i in 1:sim$Ndose){
    true1.util[i]<-sim$u00*(1-true.et$p1tox[i])*(1-true.et$p1eff[i])+100*(1-true.et$p1tox[i])*true.et$p1eff[i]+
      0*true.et$p1tox[i]*(1-true.et$p1eff[i])+sim$u11*true.et$p1tox[i]*true.et$p1eff[i]
    true2.util[i]<-sim$u00*(1-true.et$p2tox[i])*(1-true.et$p2eff[i])+100*(1-true.et$p2tox[i])*true.et$p2eff[i]+
      0*true.et$p2tox[i]*(1-true.et$p2eff[i])+sim$u11*true.et$p2tox[i]*true.et$p2eff[i]
  }
  cat(paste(sim$scen[1,],collapse = ""),"\n")
  cat("Phase1: ",true1.util,"\n")
  cat("Phase2: ",true2.util,"\n")
}

select_dose_BM<-function(sim,Ndose,pt.target,pe.null,ctox,ceff,Narm,BOIN.ldes,act,mtd,monitor.ce,a,b){
  
  nx.ce<-nx(monitor = monitor.ce)
  nv<-nx.ce$p1n+nx.ce$p2n
  xtv<-nx.ce$p1xt+nx.ce$p2xt
  xev<-nx.ce$p1xe+nx.ce$p2xe
  uv<-nx.ce$p1u+nx.ce$p2u
  
  if(sum(act)==0){
    p2dose<-0
  }else{
    p2dose<-c()
    
    nzero.dose<-sum(nv!=0)
    indc.n<-c(rep(TRUE,nzero.dose),rep(FALSE,Ndose-nzero.dose))
    #T:dose with not zero patients, F:number of patients is zero
    
    tmp.tox<-(1-pbeta(q=pt.target,shape1=a+xtv,shape2=b+(nv-xtv)))>ctox
    #T:tox rate is higher than the cutoff, F:lower
    if(any(tmp.tox)){
      tox.dose<-which.max(as.numeric(tmp.tox))
      indc.ctox<-c(rep(TRUE,tox.dose-1),rep(FALSE,Ndose-tox.dose+1))
    }else{
      indc.ctox<-rep(TRUE,Ndose)
    }
    
    indc.ceff<-!(pbeta(q=pe.null,shape1=a+xev,shape2=b+(nv-xev))>ceff)
    #T:eff rate is higher than the cutoff, F:lower
    
    indc<-indc.n&indc.ctox&indc.ceff&act
    
    max.nzero<-max(1:Ndose*indc.n)
    if(sum(indc)<Narm&&nzero.dose<Ndose&&sum(!indc.ctox)==0&&sum(!act)==0){
      if(xtv[max.nzero]/nv[max.nzero]<BOIN.ldes)indc[max.nzero+1]<-TRUE
    }
    
    for(p2i in 1:Narm){
      if(sum(indc)==0)break
      else if(sum(indc)==1){
        p2dose[p2i]<-which(indc)
        break
      }else{
        p2dose[p2i]<-maxu(indc=indc,n.vec=nv,a.vec=a+uv,b.vec=b+(nv-uv))
        indc[p2dose[p2i]]<-FALSE
      }
    }
    if(is.null(p2dose))p2dose<-0
    p2dose<-sort(p2dose)
  }
    
  p2dose
}

select_dose_33CE<-function(sim,act,mtd,monitor.ce,a,b){
  
  nx.ce<-nx(monitor = monitor.ce)
  nv<-nx.ce$p1n+nx.ce$p2n
  xtv<-nx.ce$p1xt+nx.ce$p2xt
  xev<-nx.ce$p1xe+nx.ce$p2xe
  uv<-nx.ce$p1u+nx.ce$p2u
  
  if(sum(act)==0){
    p2dose<-0
  }else{
    if(mtd==1){
      p2dose<-1
    }else if(mtd==99){
      p2dose<-0
    }else{
      p2dose<-c(mtd-1,mtd)
    }
  }
  p2dose
}

enroll_Stage2_BM<-function(chsize,pt.target,pe.null,ctox,ceff,s2,assoc,p2dose,p2ch,true.et,monitor.ce,a,b){ 
  nx.ce2<-nx(monitor.ce)
  
  act2.dose<-rep(TRUE,length(p2dose))
  if(sum(p2dose)==0)act2.dose<-FALSE
  
  cnt2.dose<-rep(0,length(p2dose))
  Np2arm<-rep(0.01,length(p2dose))
  
  if(sum(act2.dose)!=0){
    maxu.mtx<-matrix(NaN,nrow=10000,ncol=length(p2dose))
    p2l<-1
    for(p2k in p2dose){
      u<-nx.ce2$p1u[p2k]+nx.ce2$p2u[p2k]
      n<-nx.ce2$p1n[p2k]+nx.ce2$p2n[p2k]
      maxu.mtx[,p2l]<-rbeta(n=10000,shape1=a+u,shape2=b+n-u)
      p2l<-p2l+1
    }
    p.vec<-table(apply(maxu.mtx,1,which.max))/10000
  }
  
  for(p2i in 1:p2ch){
    if(sum(act2.dose)==0)break
    for(p2k in 1:length(p2dose)){
      if(!act2.dose[p2k])p.vec[p2k]<-0
    }
    currch<-rmultinom(n=1,size=chsize,prob=p.vec)
    cnt2.dose<-cnt2.dose+as.numeric(currch>0)+nx.ce2$p2n[p2dose]
    Np2arm<-Np2arm+currch
    for(p2j in 1:length(p2dose)){
      if(currch[p2j]>0){
        if(cnt2.dose[p2j]==1) monitor.ce<-monit_update(stage=2,dose=p2dose[p2j],n=currch[p2j],true.et=true.et,c=assoc,monitor=monitor.ce,append=F)
        else monitor.ce<-monit_update(stage=2,dose=p2dose[p2j],n=currch[p2j],true.et=true.et,c=assoc,monitor=monitor.ce,append=T)
        tmp.nx<-nx(monitor.ce)[p2dose[p2j],]
        if(1-pbeta(q=pt.target,shape1=a+tmp.nx$p1xt+tmp.nx$p2xt,shape2=b+(tmp.nx$p1n-tmp.nx$p1xt)+(tmp.nx$p2n-tmp.nx$p2xt))>ctox){
          Np2arm[p2j:length(p2dose)]<-10000
          act2.dose[p2j:length(p2dose)]<-FALSE
        }
      }
    }
    tmp.nx<-nx(monitor.ce)
    if(any((tmp.nx$p1n+tmp.nx$p2n)>=s2)) break
  }
  for(p2j in 1:length(p2dose)){
    if(sum(act2.dose)==0)break
    tmp.nx<-nx(monitor.ce)[p2dose[p2j],]
    if(pbeta(q=pe.null,shape1=a+tmp.nx$p1xe+tmp.nx$p2xe,shape2=b+(tmp.nx$p1n-tmp.nx$p1xe)+(tmp.nx$p2n-tmp.nx$p2xe))>ceff){
      act2.dose[p2j]<-FALSE
    }
  }
  list(p2dose=p2dose,act2.dose=act2.dose,monitor=monitor.ce)
}

enroll_Stage2<-function(sim,Np2,s2,assoc,p2dose,p2ch,true.et,monitor.ce,a,b){ 
  nx.ce2<-nx(monitor.ce)
  
  act2.dose<-rep(TRUE,length(p2dose))
  if(sum(p2dose)==0)act2.dose<-FALSE
  
  rem.Np2<-Np2
  asg.dose<-act2.dose
  for(p2j in 1:length(p2dose)){
    if(sum(act2.dose)==0)break
    n.each<-min(floor(rem.Np2/sum(asg.dose)),s2-nx.ce2$p1n[p2dose[p2j]])
    rem.Np2<-rem.Np2-n.each
    asg.dose[p2j]<-F
    monitor.ce<-monit_update(stage=2,dose=p2dose[p2j],n=n.each,true.et=true.et,c=assoc,monitor=monitor.ce,append=F)
  }
  list(p2dose=p2dose,act2.dose=act2.dose,monitor=monitor.ce)
}


dgumbel_model<-function(pe,pt,c){
  prob.mtx<-matrix(0,ncol=2,nrow=2)
  k<-1
  rxt.vec<-c()
  rxe.vec<-c()
  for(xt in 0:1){
    for(xe in 0:1){
      prob.mtx[xe+1,xt+1]<-pe^xe*(1-pe)^(1-xe)*pt^xt*(1-pt)^(1-xt)+
        pe*(1-pe)*pt*(1-pt)*(-1)^(xe+xt)*((exp(c)-1)/(exp(c)+1))
      k<-k+1
    }
  }
  rownames(prob.mtx)<-c("xe=0","xe=1")
  colnames(prob.mtx)<-c("xt=0","xt=1")
  prob.mtx
}

maxu<-function(indc,n.vec=rep(1,sim$Ndose),a.vec,b.vec){
  indc.dose<-(1:length(indc))[indc]
  
  for(maxi in 1:length(indc.dose)){
    indci<-indc.dose[maxi]
    if(maxi==1){
      if(n.vec[indci]!=0) u.sample<-rbeta(n=10000,shape1=a.vec[indci],shape2=b.vec[indci])
      else u.sample<-rep(0,10000)
    }else{
      if(n.vec[indci]!=0) u.sample<-cbind(u.sample,rbeta(n=10000,shape1=a.vec[indci],shape2=b.vec[indci]))
      else u.sample<-cbind(u.sample,rep(0,10000))
    }
  }
  maxu.vec<-apply(u.sample,1,which.max)
  max.u<-indc.dose[which.max(table(maxu.vec))]
  
  max.u
}

maxu_MEM<-function(indc,post){
  indc.dose<-(1:length(indc))[indc]
  
  for(maxi in 1:length(indc.dose)){
    indci<-indc.dose[maxi]
    
    model.sample<-matrix(NA,ncol=length(post[[indci]]$omega),nrow=10000)
    for(maxj in 1:ncol(model.sample)){
      model.sample[,maxj]<-rbeta(n=10000,shape1=post[[indci]]$a[maxj],shape2=post[[indci]]$b[maxj])
    }
    
    if(maxi==1)u.sample<-model.sample%*%matrix(post[[indci]]$omega,ncol=1)
    else u.sample<-cbind(u.sample,model.sample%*%matrix(post[[indci]]$omega,ncol=1))
  }
  maxu.vec<-apply(u.sample,1,which.max)
  max.u<-indc.dose[which.max(table(maxu.vec))]
  
  max.u
}

monit_update<-function(stage,dose,n,true.et,c,monitor,append=T){
  if(stage==1){
    n.vec<-rep(1,n)
    rxte<-rgumbel_model(n=n,pe=true.et$p1eff[dose],pt=true.et$p1tox[dose],c)
    xt.vec<-rxte$xt
    xe.vec<-rxte$xe
    
    if(!append) monitor[[dose]]$p1<-data.frame(n=n.vec,xt=xt.vec,xe=xe.vec)
    else monitor[[dose]]$p1<-rbind(monitor[[dose]]$p1,data.frame(n=n.vec,xt=xt.vec,xe=xe.vec))
  }else if(stage==2){
    n.vec<-rep(1,n)
    rxte<-rgumbel_model(n=n,pe=true.et$p2eff[dose],pt=true.et$p2tox[dose],c)
    xt.vec<-rxte$xt
    xe.vec<-rxte$xe
    
    if(!append) monitor[[dose]]$p2<-data.frame(n=n.vec,xt=xt.vec,xe=xe.vec)
    else monitor[[dose]]$p2<-rbind(monitor[[dose]]$p2,data.frame(n=n.vec,xt=xt.vec,xe=xe.vec))
  }
  monitor
}

nx<-function(monitor,Ndose=sim$Ndose,u00=sim$u00,u11=sim$u11){
  p1u<-rep(0,Ndose)
  p2u<-rep(0,Ndose)
  for(ui in 1:Ndose){
    for(uj in 1:nrow(monitor[[ui]]$p1)){
      if(monitor[[ui]]$p1$xt[uj]==0&&monitor[[ui]]$p1$xe[uj]==0)u<-u00
      if(monitor[[ui]]$p1$xt[uj]==0&&monitor[[ui]]$p1$xe[uj]==1)u<-100
      if(monitor[[ui]]$p1$xt[uj]==1&&monitor[[ui]]$p1$xe[uj]==0)u<-0
      if(monitor[[ui]]$p1$xt[uj]==1&&monitor[[ui]]$p1$xe[uj]==1)u<-u11
      if(monitor[[ui]]$p1$n[uj]!=0) p1u[ui]<-p1u[ui]+u/100
    }
    for(uj in 1:nrow(monitor[[ui]]$p2)){
      if(monitor[[ui]]$p2$xt[uj]==0&&monitor[[ui]]$p2$xe[uj]==0)u<-u00
      if(monitor[[ui]]$p2$xt[uj]==0&&monitor[[ui]]$p2$xe[uj]==1)u<-100
      if(monitor[[ui]]$p2$xt[uj]==1&&monitor[[ui]]$p2$xe[uj]==0)u<-0
      if(monitor[[ui]]$p2$xt[uj]==1&&monitor[[ui]]$p2$xe[uj]==1)u<-u11
      if(monitor[[ui]]$p2$n[uj]!=0) p2u[ui]<-p2u[ui]+u/100
    }
  }
  
  p1n<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p1$n))
  p1xt<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p1$xt))
  p1xe<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p1$xe))
  p2n<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p2$n))
  p2xt<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p2$xt))
  p2xe<-sapply(1:Ndose, function(x)sum(monitor[[x]]$p2$xe))
  data.frame(dose=1:Ndose,p1n=p1n,p1xt=p1xt,p1xe=p1xe,p1u=p1u,p2n=p2n,p2xt=p2xt,p2xe=p2xe,p2u=p2u)
}

oneside_ex<-function(Nexdose,lowside){
  if(Nexdose==0) ex<-NaN
  else{
    ex<-matrix(rep(0,Nexdose),nrow=1)
    if(lowside){
      for(memi in 1:Nexdose){
        ex<-rbind(ex,c(rep(0,Nexdose-memi),rep(1,memi)))
      }
    }
    else{
      for(memi in 1:Nexdose){
        ex<-rbind(ex,c(rep(1,memi),rep(0,Nexdose-memi)))
      }
    }
  }
  ex
}

rct_obd<-function(p2dose,act.dose,Ndose,post){
  
  if(sum(act.dose)==2){
    round.post<-rep(0,Ndose)
    round.post[p2dose[act.dose]]<-round(post,digits=3)
    obd<-which.max(round.post) #when there are more than or equal to 2 doses with a same max value, the less dose is selected
  }else if(sum(act.dose)==1){
    obd<-p2dose[act.dose]
  }else if(sum(act.dose)==0){
    obd<-0
  }
  
  obd
}

rgumbel_model<-function(n,pe,pt,c){
  prob.mtx<-matrix(0,ncol=2,nrow=2)
  k<-1
  rxt.vec<-c()
  rxe.vec<-c()
  for(xt in 0:1){
    for(xe in 0:1){
      prob.mtx[xe+1,xt+1]<-pe^xe*(1-pe)^(1-xe)*pt^xt*(1-pt)^(1-xt)+
        pe*(1-pe)*pt*(1-pt)*(-1)^(xe+xt)*((exp(c)-1)/(exp(c)+1))
      rxt.vec[k]<-xt
      rxe.vec[k]<-xe
      k<-k+1
    }
  }
  obs.mtx<-rmultinom(n=n,size=1,prob=c(prob.mtx))
  obs.comp<-sapply(1:n,function(x)which(obs.mtx[,x]==1))
  data.frame(xt=rxt.vec[obs.comp],xe=rxe.vec[obs.comp])
}


ThreePlusThree<-function(sim,true.et,monitor){
  if(sim$chsize!=3)stop("chsize is not 3")
  currj<-1
  act.dose<-c(rep(TRUE,sim$Ndose),FALSE)
  cnt.dose<-rep(0,sim$Ndose)
  add.ch<-F
  for(j in 1:(sim$Np1/sim$chsize)){
    cnt.dose[currj]<-cnt.dose[currj]+1
    
    if(cnt.dose[currj]==1) monitor<-monit_update(stage=1,dose=currj,n=sim$chsize,true.et=true.et,c=sim$assoc,monitor=monitor,append=F)
    else monitor<-monit_update(stage=1,dose=currj,n=sim$chsize,true.et=true.et,c=sim$assoc,monitor=monitor,append=T)
    
    ntotal<-nrow(monitor[[currj]]$p1)
    cjn<-sum(monitor[[currj]]$p1$n[(-3:0)+ntotal])
    cjxt<-sum(monitor[[currj]]$p1$xt[(-3:0)+ntotal])
    
    if(add.ch){
      if(cjxt+add.cjxt<=2)currj<-currj+1
      else if(cjxt+add.cjxt>=3)currj<-currj-1
      add.cjxt<-NaN
      add.ch<-F
    }else{
      if(cjxt==0)currj<-currj+1
      else if(cjxt>=2)currj<-currj-1
      else{
        add.cjxt<-cjxt
        add.ch<-T
      }
    }
    if(currj>sim$Ndose)currj<-sim$Ndose
    if(currj<1)break
    if(cjn>=sim$s1) break
  }
  list(monitor=monitor,act=act.dose,last.ch=j)
}