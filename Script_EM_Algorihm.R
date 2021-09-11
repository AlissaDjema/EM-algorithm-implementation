options( "digits"=2)

#  2. Implementation 

## 2.1 Initialisation 

K = runif(1,0,10)
theta = list()
theta$pi = rep(0,K)
theta$lambda = rep(0,K)
for (i in 1:K) {
  theta$pi[i] = c((1/K))
  theta$lambda[i] = runif(1,0,10)}

####Etape E

compute_delta = function(X,theta,K){
  n=length(X)
  num = matrix(nrow = n, ncol = K)
  delta = matrix(nrow = n, ncol = K)
  for (i in 1:n) {
    for (k in 1:K) {
      num[i,k] = (theta$pi[k])*dpois(X[i],theta$lambda[k])}}
  delta = num/rowSums(num)
  return(delta)}
delta = compute_delta(X,theta,K)

####Etape M 
updatepi = function(delta){
  n = nrow(delta)
  return(colSums(delta)/n)}

updatelambda = function(X,delta){
  lambda = rep(0,K)
  for (k in 1:K) {
    lambda[k] = sum(delta[,k] * X)/sum(delta[,k])}
  return(lambda)}

Mstep = function(X, delta) {
  theta = list()
  theta$pi = updatepi(delta)
  theta$lambda = updatelambda(X, delta)
  return(theta)}

compute_llhood = function(X, theta) {
  res = 0
  n = length(X)
  for (i in 1:n) {
    temp = 0
    for (k in 1:K) {
      temp = temp + theta$pi[k]*dpois(X[i],lambda = theta$lambda[k])}
    res = res + log(temp)}
  return(res)}

#Etape 0
llhood0 = compute_llhood(X , theta) 
count = 0
max.iter = 100
#atol = 1e-10
llhood = compute_llhood(X,theta)
logliks = c(llhood)
#conv = Inf
while (count < max.iter) {
  #prev_llhood = llhood
  count = count + 1
  
  #EM
  delta = compute_delta(X, theta,K)
  theta = Mstep(X, delta)
  #print(theta)
  
  llhood = compute_llhood(X, theta)
  logliks = c(logliks, llhood)
  #conv = abs((llhood - prev_llhood)/prev_llhood)
}

BIC = matrix(nrow = 10,ncol = 1) #je fixe K max = 10
n=length(X)
for (K in 1:10) {
  theta = list()
  theta$pi = rep(0,K)
  theta$lambda = rep(0,K)
  for (i in 1:K) {
    theta$pi[i] = c((1/K))
    theta$lambda[i] = runif(1,0,10)
  }
  count = 0
  max.iter = 200
  llhood = compute_llhood(X,theta)
  logliks = c(llhood)
  while (count < max.iter) {
    count = count + 1
    
    #EM
    delta = compute_delta(X, theta,K)
    theta = Mstep(X, delta)
  }
  llhood = compute_llhood(X, theta)
  BIC[K] = -2*llhood + (2*K-1)*log(n)
}
plot(BIC,type = "l",main="Évolution du BIC",col.main="darkblue",col.axis="darkblue",ylab="BIC",xlab="K", col.lab="darkblue",col="#66CCCC", cex.main=0.8,cex.lab=0.8)

K = which.min(BIC)
count = 0
max.iter = 200
atol = 1e-10
llhood = compute_llhood(X,theta)
logliks = c(llhood)
while (count < max.iter) {
  count = count + 1
  
  #EM
  delta = compute_delta(X, theta,K)
  theta = Mstep(X, delta)

  llhood = compute_llhood(X, theta)
}
print(theta)
