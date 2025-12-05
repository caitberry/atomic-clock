oneTrial = function(dat, GGmodel, N,tau,mu,num.iters,alpha,nu,howDatSim,rsessionNum,adaptDeltaVal,maxTreeDepthVal,tdf){ #easy wrapper to catch errors and keep going
  tmp_df = tryCatch(runOneTrial(dat, GGmodel, N,tau,mu,num.iters,alpha,nu,howDatSim,rsessionNum,adaptDeltaVal,maxTreeDepthVal,tdf),
                    error = function(err){print("err"); return(data.frame("method"=NA))})
  return (tmp_df)
}

# FirstDatSim=simulateDataHSSG(N = N,mu = mu,tau = tau,alpha=alpha,nu=nu)
# 
dat = FirstDatSim

# GGmodel = myHGG_noDoF;LGmodel = myHLG_noDoF;
# N = N;
# tau = tau;mu = mu; num.iters = num.iters;alpha = alpha;nu = nu; 
# howDatSim = howDatSim;rsessionNum = rsessionNum;
# adaptDeltaVal=theadaptdelta;maxTreeDepthVal=thetreedepth;tdf=thetdf

runOneTrial = function(dat, GGmodel, N,tau,mu,num.iters,alpha,nu,howDatSim,rsessionNum,adaptDeltaVal,maxTreeDepthVal,tdf){

  # sttime=Sys.time()
  # print("running RMA")
  metaforRes=analyzeDataRMA(data = dat,trueTau = tau,trueMu = mu)
  # print(Sys.time()-sttime)

  
  CSRes=analyzeData_TFmethods(data = dat,trueTau = tau,trueMu = mu)
  
    
  stan_data <- list(N = dim(dat)[1],
                    x = dat$x,
                    u=dat$u,
                    tdf=tdf)
  
  GGres=fitGaussGauss(stanDat = stan_data,trueTau = tau,trueMu = mu, model = GGmodel,n.iter=num.iters,adapt.delta = adaptDeltaVal,tree.depth=maxTreeDepthVal)
  bayesRes=GGres
  
  newOut=bind_rows(metaforRes,CSRes,bayesRes)
  newOut$rsession=rsessionNum
  newOut$simulation=howDatSim
  newOut$N=N
  newOut$mu=mu
  newOut$tau=tau
  newOut$num.iters=num.iters
  newOut$alpha=alpha
  newOut$nu=nu
  newOut$tdf=tdf
  newOut$maxtree=maxTreeDepthVal
  newOut$adaptDelta=adaptDeltaVal
  
  return(newOut)
}
