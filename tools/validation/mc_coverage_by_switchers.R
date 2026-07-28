suppressPackageStartupMessages(library(didgpu))
for(NSW in c(2L,4L,8L,13L)){
  cov1<-0L; est<-numeric(200); ses<-numeric(200)
  for(s in 1:200){
    set.seed(2000+s); nU<-30L; Tn<-20L
    unit<-rep(1:nU,each=Tn); period<-rep(1:Tn,nU)
    D<-as.integer(unit<=NSW & period>=10L)
    Y<-rnorm(nU)[unit]+0.05*period+0.4*D+rnorm(nU*Tn,0,0.5)
    f<-tryCatch(suppressWarnings(didgpu(df=data.frame(unit,period,Y,D),outcome="Y",group="unit",time="period",
        treatment="D",effects=5L,placebo=3L,cluster="unit",bootstrap_reps=300L,seed=s,backend="cpu",verbose=FALSE)),error=function(e)NULL)
    if(is.null(f)){est[s]<-NA;next}
    E<-f$results$Effects; est[s]<-E[1,"Estimate"]; ses[s]<-E[1,"SE"]
    if(is.finite(E[1,"LB.CI"])&&E[1,"LB.CI"]<=0.4&&E[1,"UB.CI"]>=0.4) cov1<-cov1+1L
  }
  ok<-sum(is.finite(est))
  cat(sprintf("switchers=%2d | ok=%d | bias=%+.3f | bootSE=%.3f | empSD=%.3f | coverage=%.1f%%\n",
      NSW, ok, mean(est,na.rm=TRUE)-0.4, mean(ses,na.rm=TRUE), sd(est,na.rm=TRUE), 100*cov1/ok))
}
