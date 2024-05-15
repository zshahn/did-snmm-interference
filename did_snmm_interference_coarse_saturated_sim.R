library(nleqslv)
# setwd('C:/Users/zshahn/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
# setwd('~/Dropbox/Causal/self_controlled/application/favara_imbs_data/')
#Generate Data
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}

#true saturated blip function used to generate data
blip = function(m,k,l=0,a,past_a,psi){
  if(m==1){
    a[1]*psi[1] + a[2]*psi[2] + a[1]*(k-m)*psi[3] + a[2]*(k-m)*psi[4] + a[1]*a[2]*psi[5] +
    a[1]*a[2]*(k-m)*psi[6] 
    }else{
      a[1]*psi[7] + a[2]*psi[8] + a[1]*a[2]*psi[9] + past_a[1]*a[2]*psi[10] + past_a[2]*a[1]*psi[11] +
      past_a[2]*a[2]*psi[12] + past_a[2]*a[1]*a[2]*psi[13]
    }
}

psi=c(1,.5,-.1,-.1,-.2,-.05,
      1,.5,-.1,-.1,-.1,-.05,-0.05)

N = 10000
nsims=1000
effects_hat_list = vector(mode='list',length=nsims)
psi_hat_list = vector(mode='list',length=nsims)

#make U correlated across neighboring units
# U = c(rep(1,N/2),rep(0,N/2))
#make U uncorrelated
for(i in 1:nsims){
U = rbinom(N,1,.5)
Y0 = rnorm(N,U,.1)
A1 = rbinom(N,1,.3+.2*U)
lag1_A1 = c(0,A1[1:(N-1)])
jump1_A1 = c(A1[2:N],0)
phi_A1 = cbind(A1,pmax(lag1_A1,jump1_A1))
colnames(phi_A1) = c('A1','int1')
Y1_0 = rnorm(N,U,.1)
Y1 = rnorm(N,Y1_0 + sapply(1:N,function(i) blip(m=1,k=1,a=phi_A1[i,],past_a=c(0,0),psi=psi)),.1)
A2 = ifelse(A1==1,0,rbinom(N,1,.3+ .2*U))
lag1_A2 = c(0,A2[1:(N-1)])
jump1_A2 = c(A2[2:N],0)
phi_A2 = cbind(A2,pmax(lag1_A2,jump1_A2))
colnames(phi_A2) = c('A2','int2')
Y2_0 = rnorm(N,U,.1)
Y2 = rnorm(N,Y2_0 + sapply(1:N,function(i) blip(m=1,k=2,a=phi_A1[i,],past_a=c(0,0),psi=psi)) +
             sapply(1:N,function(i) blip(m=2,k=2,a=phi_A2[i,],past_a=phi_A1[i,],psi=psi)),.1)

data2 = data.frame(cbind(phi_A1,phi_A2,Y0,Y1,Y2,Y2_0,Y1_0))
data2$id = 1:N
ids = unique(data2$id)
Hmk_mat = data.frame(cbind(id=rep(ids,each=5),m=rep(rep(1:2,3:2),length(ids)),
                           k=rep(c(0:2,1:2),length(ids))))
#add treatment variables
Hmk_mat = merge(Hmk_mat,data2[,c('id','A1','int1','A2','int2')])

#add covariates (we don't have any in simplest simulation)

#add outcomes by k
Hmk_mat = merge(Hmk_mat,data2[,c('id','Y0','Y1','Y2')])
Hmk_mat$Y = ifelse(Hmk_mat$k==0,Hmk_mat$Y0,ifelse(Hmk_mat$k==1,Hmk_mat$Y1,Hmk_mat$Y2))

#add treatments by m
Hmk_mat$A = ifelse(Hmk_mat$m==1,Hmk_mat$A1,Hmk_mat$A2)
Hmk_mat$int = ifelse(Hmk_mat$m==1,Hmk_mat$int1,Hmk_mat$int2)


Hmk_mat = Hmk_mat[order(Hmk_mat$id,Hmk_mat$m,Hmk_mat$k),]


#nonparametric treatment models for A1 and int1
data2$A1_hat = mean(data2$A1)
data2$int1_hat = data2$A1_hat*mean(data2$int1[data2$A1==1]) +
  (1-data2$A1_hat)*mean(data2$int1[data2$A1==0])
#parametric models for A2 and int2
A2_mod = glm(A2~int1,data=data2[data2$A1==0,],family='binomial')
data2$A2_hat = 0
data2$A2_hat[data2$A1==0] = predict(A2_mod,newdata = data2[data2$A1==0,],type='response')
int2_mod = lm(int2~A2*A1*int1,data=data2)
data2$int2_hat = data2$A2_hat*predict(int2_mod,newdata=data.frame(A1=data2$A1,int1=data2$int1,A2=1)) +
  (1-data2$A2_hat)*predict(int2_mod,newdata=data.frame(A1=data2$A1,int1=data2$int1,A2=0)) 

data2$A1int1_hat = mean(data2$A1*data2$int1)

data2$A2int2_hat = 0
data2$A2int2 = data2$A2*data2$int2
A2int2_mod = lm(A2int2~int1,data=data2[data2$A1==0 & data2$int1!=2,])
data2$A2int2_hat[data2$A1==0 & data2$int1!=2] = predict(A2int2_mod,newdata = data2[data2$A1==0 & data2$int1!=2,],type='response')


data2$m = 1
data2$A_hat = data2$A1_hat
Hmk_mat = merge(Hmk_mat,rbind(data.frame(id=data2$id,m=1,A_hat=data2$A1_hat),
                              data.frame(id=data2$id,m=2,A_hat=data2$A2_hat)),all.x=T)

Hmk_mat = merge(Hmk_mat,rbind(data.frame(id=data2$id,m=1,int_hat=data2$int1_hat),
                              data.frame(id=data2$id,m=2,int_hat=data2$int2_hat)),all.x=T)

Hmk_mat = merge(Hmk_mat,rbind(data.frame(id=data2$id,m=1,Aint_hat=data2$A1int1_hat),
                              data.frame(id=data2$id,m=2,Aint_hat=data2$A2int2_hat)),all.x=T)

Hmk_mat = Hmk_mat[order(Hmk_mat$id,Hmk_mat$m,Hmk_mat$k),]

Hmk_mat$past_A = ifelse(Hmk_mat$m==1,0,Hmk_mat$A1)
Hmk_mat$past_A1 = 0
Hmk_mat$past_A2 = Hmk_mat$A1

Hmk_mat$past_int = ifelse(Hmk_mat$m==1,0,Hmk_mat$int1)
Hmk_mat$past_int1 = 0
Hmk_mat$past_int2 = Hmk_mat$int1


Hmk_mat$S1 = Hmk_mat$A*(Hmk_mat$m==1)
Hmk_mat$S2 = Hmk_mat$int*(Hmk_mat$m==1)
Hmk_mat$S3 = Hmk_mat$A*(Hmk_mat$k-Hmk_mat$m)*(Hmk_mat$m==1)
Hmk_mat$S4 = Hmk_mat$int*(Hmk_mat$k-Hmk_mat$m)
Hmk_mat$S5 = Hmk_mat$A*Hmk_mat$int*(Hmk_mat$m==1)
Hmk_mat$S6 = Hmk_mat$A*Hmk_mat$int*(Hmk_mat$k-Hmk_mat$m)

Hmk_mat$S7 = Hmk_mat$A*(Hmk_mat$m==2)
Hmk_mat$S8 = Hmk_mat$int*(Hmk_mat$m==2)
Hmk_mat$S9 = Hmk_mat$A*Hmk_mat$int*(Hmk_mat$m==2)
Hmk_mat$S10 = Hmk_mat$past_A*Hmk_mat$int*(Hmk_mat$m==2)
Hmk_mat$S11 = Hmk_mat$past_int*Hmk_mat$A*(Hmk_mat$m==2)
Hmk_mat$S12 = Hmk_mat$int*Hmk_mat$past_int*(Hmk_mat$m==2)
Hmk_mat$S13 = Hmk_mat$int*Hmk_mat$A*Hmk_mat$past_int*(Hmk_mat$m==2)


years = min(Hmk_mat$m):max(Hmk_mat$m)
compute_est_eq_psi = function(psi,Hmk=Hmk_mat){
  psi_mat = matrix(rep(psi,each=nrow(Hmk)),nrow=nrow(Hmk))
  gammas = matrix(NA,ncol=length(years),nrow=nrow(Hmk))
  # gammas2 = matrix(NA,ncol=length(years),nrow=nrow(Hmk))
  
  colnames(gammas) = paste0('gamma_',years,'_k')
  for(j in 1:length(years)){
    gammas[,paste0('gamma_',years[j],'_k')] = ifelse((Hmk$k<years[j])|(Hmk$m>years[j]),0,rowSums(cbind(Hmk[,paste0("A",years[j])]*(years[j]==1),Hmk[,paste0("int",years[j])]*(years[j]==1),Hmk[,paste0("A",years[j])]*(Hmk$k-years[j])*(years[j]==1),Hmk[,paste0("int",years[j])]*(Hmk$k-years[j])*(years[j]==1),
                                                                                                       Hmk[,paste0("int",years[j])]*Hmk[,paste0("A",years[j])]*(years[j]==1),Hmk[,paste0("int",years[j])]*Hmk[,paste0("A",years[j])]*(Hmk$k-years[j])*(years[j]==1),
                                                                                                       Hmk[,paste0("A",years[j])]*(years[j]==2),Hmk[,paste0("int",years[j])]*(years[j]==2),Hmk[,paste0("int",years[j])]*Hmk[,paste0("A",years[j])]*(years[j]==2),
                                                                                                       Hmk[,paste0("int",years[j])]*Hmk[,paste0("past_A",years[j])]*(years[j]==2),Hmk[,paste0("A",years[j])]*Hmk[,paste0("past_int",years[j])]*(years[j]==2),Hmk[,paste0("int",years[j])]*Hmk[,paste0("past_int",years[j])]*(years[j]==2),Hmk[,paste0("int",years[j])]*Hmk[,paste0("A",years[j])]*Hmk[,paste0("past_int",years[j])]*(years[j]==2))*psi_mat))
  }
  # gammas2[,1] = ifelse(Hmk$k<1|Hmk$m>1,0,sapply(1:nrow(Hmk),function(i)blip(m=1,k=Hmk$k[i],a=c(Hmk$A1[i],Hmk$int1[i]),past_a=c(0,0),psi=psi)))
  # gammas2[,2] = ifelse(Hmk$k<2|Hmk$m>2,0,sapply(1:nrow(Hmk),function(i)blip(m=2,k=Hmk$k[i],a=c(Hmk$A2[i],Hmk$int2[i]),past_a=c(Hmk$A1,Hmk$int1),psi=psi)))
  # # 
  # gammas[,'gamma_1_k'] = ifelse(Hmk$m>1|Hmk$k<1,0,sapply(1:nrow(Hmk),function(i) blip(m=1,k=Hmk$k[i],a=c(Hmk$A1[i],Hmk$int1[i]),past_a=c(0,0),psi=psi[i,])))
  # gammas[,'gamma_2_k'] = ifelse(Hmk$k<2,0,sapply(1:nrow(Hmk),function(i) blip(m=2,k=Hmk$k[i],a=c(Hmk$A2[i],Hmk$int2[i]),past_a=c(Hmk$A1,Hmk$int1),psi=psi[i,])))
  gamma = rowSums(gammas)
  Hmk$H = Hmk$Y - gamma
  Hmk = Hmk %>%
    group_by(id,m) %>%
    dplyr::mutate(H_lag = lag(H, n = 1, default = NA))  
  Hmk$H_diffs = Hmk$H - Hmk$H_lag
  H_mod1 = lm(H_diffs~k,data=Hmk[Hmk$k>=Hmk$m & Hmk$m==1,])
  H_mod2 = lm(H_diffs~A1*int1,data=Hmk[Hmk$k>=Hmk$m & Hmk$m==2,])
  Hmk$V = ifelse(Hmk$m==1,predict(H_mod1,newdata=Hmk),predict(H_mod2,newdata=Hmk))
  est_mat = Hmk[Hmk$k>=Hmk$m & Hmk$m>0,]
  temp = (est_mat$H_diffs-est_mat$V)*(est_mat[,c('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13')]-cbind(est_mat$A_hat*(est_mat$m==1),est_mat$int_hat*(est_mat$m==1),est_mat$A_hat*(est_mat$k-est_mat$m)*(est_mat$m==1),est_mat$int_hat*(est_mat$k-est_mat$m)*(est_mat$m==1),est_mat$Aint_hat*(est_mat$m==1),est_mat$Aint_hat*(est_mat$k-est_mat$m)*(est_mat$m==1),
                                                                                                                        est_mat$A_hat*(est_mat$m==2),est_mat$int_hat*(est_mat$m==2),est_mat$Aint_hat*(est_mat$m==2),est_mat$int_hat*est_mat$past_A*(est_mat$m==2),est_mat$A_hat*est_mat$past_int*(est_mat$m==2),est_mat$int_hat*est_mat$past_int*(est_mat$m==2),est_mat$Aint_hat*est_mat$past_int*(est_mat$m==2)))
  colSums(temp)
}

ss = nleqslv(x=rep(0,13),fn=compute_est_eq_psi)
psi_hat = ss$x
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



