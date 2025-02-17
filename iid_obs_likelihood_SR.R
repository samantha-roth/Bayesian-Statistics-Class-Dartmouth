

log.lik<- function(params){
  mod.preds<- model(c(param[1],param[2]),t)$mod.obs
  resid<- mod.preds-observations
  -N/2*log(2*pi) - N/2*log(params[3]^2) - 1/2*sum(resid^2)/params[3]^2
}



modelpreds<- model(c(result$par[1],result$par[2]),t)$mod.obs