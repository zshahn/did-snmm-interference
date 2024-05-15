library(nleqslv)
# setwd('C:/Users/zshahn/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
# setwd('~/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
#Generate Data
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}



psi11_t=c(1,.5)
psi12_t=c(2,1)
psi22_t=c(.75,.25,.1)
psi= c(psi11_t,psi12_t,psi22_t)

blip111 = function(a11,a21,psi11=c(1,.5)){
  psi11[1]*a11 + psi11[2]*a21
}

blip211 = function(a11,a21,psi11=c(1,.5)){
  psi11[1]*a21 + psi11[2]*a11
}

blip112 = function(a11,a21,psi12=c(2,1)){
  psi12[1]*a11 + psi12[2]*a21
}

blip212 = function(a11,a21,psi12=c(2,1)){
  psi12[1]*a21 + psi12[2]*a11
}

blip122 = function(a12,a22,a11,a21,psi22=c(.75,.25,.1)){
  psi22[1]*a12 + psi22[2]*a22 + psi22[3]*a12*(a21+a22)
}

blip222 = function(a12,a22,a11,a21,psi22=c(.75,.25,.1)){
  psi22[1]*a22 + psi22[2]*a12 + psi22[3]*a22*(a11+a12)
}


N = 10000
nsims=1000
effects_hat_list = vector(mode='list',length=nsims)
psi_hat_list = vector(mode='list',length=nsims)

#make U correlated across neighboring units
# U = c(rep(1,N/2),rep(0,N/2))
#make U uncorrelated

for(i in 1:nsims){
  U = rbinom(N,1,.5)
  Y10 = rnorm(N,U,.1)
  Y20 = rnorm(N,U,.1)
  A11 = rbinom(N,1,.3+.2*U)
  A21 = rbinom(N,1,.3+.2*U)
  Y11_0 = rnorm(N,U,.1)
  Y21_0 = rnorm(N,U,.1)
  Y11 = rnorm(N,Y11_0 + sapply(1:N,function(i) blip111(a11=A11[i],a21=A21[i],psi11=psi11_t)),.1)
  Y21 = rnorm(N,Y21_0 + sapply(1:N,function(i) blip211(a11=A11[i],a21=A21[i],psi11=psi11_t)),.1)
  
  A12 = ifelse(A11==1,0,rbinom(N,1,.3+ .2*U))
  A22 = ifelse(A21==1,0,rbinom(N,1,.3+ .2*U))
  Y12_0 = rnorm(N,U,.1)
  Y22_0 = rnorm(N,U,.1)
  Y12 = rnorm(N,Y12_0 + sapply(1:N,function(i) blip112(a11=A11[i],a21=A21[i],psi12=psi12_t)) +
               sapply(1:N,function(i) blip122(a12=A12[i],a22=A22[i],a11=A11[i],a21=A21[i],psi22=psi22_t)),.1)
  Y22 = rnorm(N,Y22_0 + sapply(1:N,function(i) blip212(a11=A11[i],a21=A21[i],psi12=psi12_t)) +
                sapply(1:N,function(i) blip222(a12=A12[i],a22=A22[i],a11=A11[i],a21=A12[i],psi22=psi22_t)),.1)
  
  data2 = data.frame(cbind(A11,A21,A12,A22,Y10,Y20,Y11,Y21,Y12,Y22,Y11_0,Y21_0,Y12_0,Y22_0))
  data2$id = 1:N
  ids = unique(data2$id)

  
  
  #nonparametric treatment models
  data2$A11_hat = mean(data2$A11)
  data2$A21_hat = mean(data2$A21)
  A12_mod = glm(A12~A21,data=data2[data2$A11==0,],family='binomial')
  data2$A12_hat = 0
  data2$A12_hat[data2$A11==0] = predict(A12_mod,newdata = data2[data2$A11==0,],type='response')
  A22_mod = lm(A22~A11,data=data2[data2$A21==0,])
  data2$A22_hat = 0
  data2$A22_hat[data2$A21==0] = predict(A22_mod,newdata = data2[data2$A21==0,],type='response')
  
  data2$A12sumA2_hat = 0
  data2$A12sumA2 = data2$A12*(data2$A21+data2$A22)
  A12sumA2_mod = lm(A12sumA2~A21,data=data2[data2$A11==0,])
  data2$A12sumA2_hat[data2$A11==0] = predict(A12sumA2_mod,newdata = data2[data2$A11==0,])
  
  data2$A22sumA1_hat = 0
  data2$A22sumA1 = data2$A22*(data2$A11+data2$A12)
  A22sumA1_mod = lm(A22sumA1~A11,data=data2[data2$A21==0,])
  data2$A22sumA1_hat[data2$A21==0] = predict(A22sumA1_mod,newdata = data2[data2$A21==0,])
  
  # data2$A1_hat = data2$A11_hat
  
  data2$qi1m1k1_1 = data2$A11 
  data2$qi1m1k1_1_hat = data2$A11_hat
  data2$qi1m1k1_2 = data2$A21
  data2$qi1m1k1_2_hat = data2$A21_hat
  
  data2$qi2m1k1_1 = data2$A21 
  data2$qi2m1k1_1_hat = data2$A21_hat
  data2$qi2m1k1_2 = data2$A11
  data2$qi2m1k1_2_hat = data2$A11_hat
  
  data2$qi1m1k2_1 = data2$A11 
  data2$qi1m1k2_1_hat = data2$A11_hat
  data2$qi1m1k2_2 = data2$A21
  data2$qi1m1k2_2_hat = data2$A21_hat
  
  data2$qi2m1k2_1 = data2$A21 
  data2$qi2m1k2_1_hat = data2$A21_hat
  data2$qi2m1k2_2 = data2$A11
  data2$qi2m1k2_2_hat = data2$A11_hat
  
  data2$qi1m2k2_1 = data2$A12 
  data2$qi1m2k2_1_hat = data2$A12_hat
  data2$qi1m2k2_2 = data2$A22
  data2$qi1m2k2_2_hat = data2$A22_hat
  data2$qi1m2k2_3 = data2$A12*(data2$A22+data2$A21)
  data2$qi1m2k2_3_hat = data2$A12sumA2_hat
  
  data2$qi2m2k2_1 = data2$A22 
  data2$qi2m2k2_1_hat = data2$A22_hat
  data2$qi2m2k2_2 = data2$A12
  data2$qi2m2k2_2_hat = data2$A12_hat
  data2$qi2m2k2_3 = data2$A22*(data2$A12+data2$A11)
  data2$qi2m2k2_3_hat = data2$A22sumA1_hat


  compute_est_eq_psi = function(psi,estmat=data2){
    psi11 = psi[1:2]
    psi12 = psi[3:4]
    psi22 = psi[5:7]
    estmat$H110 = estmat$Y10
    estmat$H210 = estmat$Y20
    estmat$H111 = estmat$Y11 - (psi11[1]*estmat$A11 + psi11[2]*estmat$A21)
    estmat$H211 = estmat$Y21 - (psi11[1]*estmat$A21 + psi11[2]*estmat$A11)
    estmat$H112 = estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22)) -
      (psi12[1]*estmat$A11 + psi12[2]*estmat$A21)
    estmat$H212 = estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12)) -
      (psi12[1]*estmat$A21 + psi12[2]*estmat$A11)
    estmat$H122 = estmat$Y12 - (psi22[1]*estmat$A12 + psi22[2]*estmat$A22 + psi22[3]*estmat$A12*(estmat$A21+estmat$A22))
    estmat$H222 = estmat$Y22 - (psi22[1]*estmat$A22 + psi22[2]*estmat$A12 + psi22[3]*estmat$A22*(estmat$A11+estmat$A12))
    estmat$H121 = estmat$Y11
    estmat$H221 = estmat$Y21
    
    estmat$H_diff111 = estmat$H111 - estmat$H110
    estmat$H_diff211 = estmat$H211 - estmat$H210
    estmat$H_diff112 = estmat$H112 - estmat$H111
    estmat$H_diff212 = estmat$H212 - estmat$H211
    estmat$H_diff122 = estmat$H122 - estmat$H121
    estmat$H_diff222 = estmat$H222 - estmat$H221
    
    H_mod111 = lm(H_diff111~1,data=estmat)
    H_mod211 = lm(H_diff211~1,data=estmat)
    H_mod112 = lm(H_diff112~1,data=estmat)
    H_mod212 = lm(H_diff212~1,data=estmat)
    H_mod122 = lm(H_diff122~A11*A21,data=estmat)
    H_mod222 = lm(H_diff222~A11*A21,data=estmat)
    
    estmat$V111 = predict(H_mod111,newdata=estmat)
    estmat$V211 = predict(H_mod211,newdata=estmat)
    estmat$V112 = predict(H_mod112,newdata=estmat)
    estmat$V212 = predict(H_mod212,newdata=estmat)
    estmat$V122 = predict(H_mod122,newdata=estmat)
    estmat$V222 = predict(H_mod222,newdata=estmat)
    
    eqs = cbind(c((estmat$H_diff111-estmat$V111)*(estmat$qi1m1k1_1-estmat$qi1m1k1_1_hat),
                (estmat$H_diff211-estmat$V211)*(estmat$qi2m1k1_1-estmat$qi2m1k1_1_hat)),
                c((estmat$H_diff111-estmat$V111)*(estmat$qi1m1k1_2-estmat$qi1m1k1_2_hat),
                (estmat$H_diff211-estmat$V211)*(estmat$qi2m1k1_2-estmat$qi2m1k1_2_hat)),
                c((estmat$H_diff112-estmat$V112)*(estmat$qi1m1k2_1-estmat$qi1m1k2_1_hat),
                  (estmat$H_diff212-estmat$V212)*(estmat$qi2m1k2_1-estmat$qi2m1k2_1_hat)),
                c((estmat$H_diff112-estmat$V112)*(estmat$qi1m1k2_2-estmat$qi1m1k2_2_hat),
                (estmat$H_diff212-estmat$V212)*(estmat$qi2m1k2_2-estmat$qi2m1k2_2_hat)),
                c((estmat$H_diff122-estmat$V122)*(estmat$qi1m2k2_1-estmat$qi1m2k2_1_hat),
                  (estmat$H_diff222-estmat$V222)*(estmat$qi2m2k2_1-estmat$qi2m2k2_1_hat)),
                c((estmat$H_diff122-estmat$V122)*(estmat$qi1m2k2_2-estmat$qi1m2k2_2_hat),
                  (estmat$H_diff222-estmat$V222)*(estmat$qi2m2k2_2-estmat$qi2m2k2_2_hat)),
                c((estmat$H_diff122-estmat$V122)*(estmat$qi1m2k2_3-estmat$qi1m2k2_3_hat),
                (estmat$H_diff222-estmat$V222)*(estmat$qi2m2k2_3-estmat$qi2m2k2_3_hat)))                                                                                                               
    colSums(eqs)
  }
  
  ss = nleqslv(x=rep(0,7),fn=compute_est_eq_psi)
  ss$termcd
  psi_hat = ss$x
  psi_hat
  psi_hat_list[[i]] = psi_hat
  effects_hat_list[[i]] = c(blip(m=1,k=1,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                            blip(m=1,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                            blip(m=1,k=1,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                            blip(m=1,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                            blip(m=1,k=1,a=c(0,1),past_a=c(0,0),psi=psi_hat),
                            blip(m=1,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat),
                            
                            blip(m=2,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
                            blip(m=2,k=2,a=c(1,0),past_a=c(0,1),psi=psi_hat),
                            blip(m=2,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
                            blip(m=2,k=2,a=c(1,1),past_a=c(0,1),psi=psi_hat),
                            blip(m=2,k=2,a=c(0,1),past_a=c(0,1),psi=psi_hat),
                            blip(m=2,k=2,a=c(0,1),past_a=c(1,0),psi=psi_hat),
                            blip(m=2,k=2,a=c(0,1),past_a=c(1,1),psi=psi_hat),
                            blip(m=2,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat))
}

save(effects_hat_list,psi_hat_list,file="~/snmm/interference/interference_sim_results.RData")

mean(sapply(effects_hat_list,function(x)x[1]))
sd(sapply(effects_hat_list,function(x)x[1]))

mean(sapply(effects_hat_list,function(x)x[2]))
sd(sapply(effects_hat_list,function(x)x[2]))

mean(sapply(effects_hat_list,function(x)x[3]))
sd(sapply(effects_hat_list,function(x)x[3]))

mean(sapply(effects_hat_list,function(x)x[4]))
sd(sapply(effects_hat_list,function(x)x[4]))

mean(sapply(effects_hat_list,function(x)x[5]))
sd(sapply(effects_hat_list,function(x)x[5]))

mean(sapply(effects_hat_list,function(x)x[6]))
sd(sapply(effects_hat_list,function(x)x[6]))

mean(sapply(effects_hat_list,function(x)x[7]))
sd(sapply(effects_hat_list,function(x)x[7]))

mean(sapply(effects_hat_list,function(x)x[8]))
sd(sapply(effects_hat_list,function(x)x[8]))

mean(sapply(effects_hat_list,function(x)x[9]))
sd(sapply(effects_hat_list,function(x)x[9]))

mean(sapply(effects_hat_list,function(x)x[10]))
sd(sapply(effects_hat_list,function(x)x[10]))

mean(sapply(effects_hat_list,function(x)x[11]))
sd(sapply(effects_hat_list,function(x)x[11]))

mean(sapply(effects_hat_list,function(x)x[12]))
sd(sapply(effects_hat_list,function(x)x[12]))

mean(sapply(effects_hat_list,function(x)x[13]))
sd(sapply(effects_hat_list,function(x)x[13]))

mean(sapply(effects_hat_list,function(x)x[14]))
sd(sapply(effects_hat_list,function(x)x[14]))

#all the conditional effects
effect_hats=c(blip(m=1,k=1,a=c(1,0),past_a=c(0,0),psi=psi_hat),
              blip(m=1,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
              blip(m=1,k=1,a=c(1,1),past_a=c(0,0),psi=psi_hat),
              blip(m=1,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
              blip(m=1,k=1,a=c(0,1),past_a=c(0,0),psi=psi_hat),
              blip(m=1,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat),
              
              blip(m=2,k=2,a=c(1,0),past_a=c(0,0),psi=psi_hat),
              blip(m=2,k=2,a=c(1,0),past_a=c(0,1),psi=psi_hat),
              blip(m=2,k=2,a=c(1,1),past_a=c(0,0),psi=psi_hat),
              blip(m=2,k=2,a=c(1,1),past_a=c(0,1),psi=psi_hat),
              blip(m=2,k=2,a=c(0,1),past_a=c(0,1),psi=psi_hat),
              blip(m=2,k=2,a=c(0,1),past_a=c(1,0),psi=psi_hat),
              blip(m=2,k=2,a=c(0,1),past_a=c(1,1),psi=psi_hat),
              blip(m=2,k=2,a=c(0,1),past_a=c(0,0),psi=psi_hat)
)

true_effects=c(blip(m=1,k=1,a=c(1,0),past_a=c(0,0),psi=psi),
               blip(m=1,k=2,a=c(1,0),past_a=c(0,0),psi=psi),
               blip(m=1,k=1,a=c(1,1),past_a=c(0,0),psi=psi),
               blip(m=1,k=2,a=c(1,1),past_a=c(0,0),psi=psi),
               blip(m=1,k=1,a=c(0,1),past_a=c(0,0),psi=psi),
               blip(m=1,k=2,a=c(0,1),past_a=c(0,0),psi=psi),
               
               blip(m=2,k=2,a=c(1,0),past_a=c(0,0),psi=psi),
               blip(m=2,k=2,a=c(1,0),past_a=c(0,1),psi=psi),
               blip(m=2,k=2,a=c(1,1),past_a=c(0,0),psi=psi),
               blip(m=2,k=2,a=c(1,1),past_a=c(0,1),psi=psi),
               blip(m=2,k=2,a=c(0,1),past_a=c(0,1),psi=psi),
               blip(m=2,k=2,a=c(0,1),past_a=c(1,0),psi=psi),
               blip(m=2,k=2,a=c(0,1),past_a=c(1,1),psi=psi),
               blip(m=2,k=2,a=c(0,1),past_a=c(0,0),psi=psi)
)



