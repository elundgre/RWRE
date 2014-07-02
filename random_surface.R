#create random surface
#packages
require(Matrix)
source('~/USC/2014_02_Summer/iitoij.R')
# grid size
n <- 2
# smoothing size
nn <- 3

# white noise
V_iid <- array(dim=c(n+2*nn,n+2*nn),replicate((n+2*nn)^2,rnorm(1,0,9)))
#smoothed matrix
V <- array(dim=c(n,n),0)
for(i in 1:n){
  for(j in 1:n){
    V[i,j] <- mean(V_iid[i:(i+2*nn),j:(j+2*nn)])
  }
}


image(V)

#array for transition rates (q[,,1] is right, q[,,2] is up...  q[,,5] is stay)
q <- array(0,dim=c(n,n,5))
p <- array(0,dim=c(n,n,4)) # transition probabilities

for(i in 1:n){
  for(j in 1:n){
    if(i != n){q[i,j,1] <- exp(V[i,j]-V[i+1,j])} # right
    if(j != n){q[i,j,2] <- exp(V[i,j]-V[i,j+1])} # up
    if(i != 1){q[i,j,3] <- exp(V[i,j]-V[i-1,j])} # left
    if(j != 1){q[i,j,4] <- exp(V[i,j]-V[i,j-1])} # down
    q[i,j,5] <- -sum(q[i,j,1:4])                 # stay
    p[i,j,1:4] <- q[i,j,1:4]/(-q[i,j,5])         # probabilities
  }
}

m = 10  # number of steps
# declare S
S <- array(0,dim=c(m+1,2))
f <- array(0,dim=c(n,n)) # empirical distribution
S[1,] <- c(round(n/2),round(n/2)) # start in center
f[S[1,1],S[1,2]] <- f[S[1,1],S[1,2]] + 1

for(i in 1:m){
  x <- rmultinom(1,1,p[S[i,1],S[i,2],])
  S[i+1,] <- S[i,] + x[1]*c(1,0) + x[2]*c(0,1) +
    x[3]*c(-1,0) + x[4]*c(0,-1)
  f[S[i,1],S[i,2]] <- f[S[i,1],S[i,2]] + 1
}

#image(f/(m+1))

# declare transition/generating matrix
rightii <- 1:(n*(n-1))
rightjj <- rightii + n

upii <- 1:(n-1)
for(i in 2:n){
  upii <- c(upii,((i-1)*(n)+1):(i*n-1))
}
upjj <- upii + 1

leftii <- (n+1):(n^2)
leftjj <- leftii - n

downii <- 2:n
for(i in 2:n){
  downii <- c(downii,((i-1)*(n)+2):(i*n))
}
downjj <- downii - 1

# convert transition values into form that can be put into sparse matrix

rightx = q[iitoi(rightii[1],n),iitoj(rightii[1],n),1]
for(i in 2:length(rightii)){
  rightx = c(rightx,q[iitoi(rightii[i],n),iitoj(rightii[i],n),1])
}

upx = q[iitoi(upii[1],n),iitoj(upii[1],n),2]
for(i in 2:length(upii)){
  upx = c(upx,q[iitoi(upii[i],n),iitoj(upii[i],n),2])
}

leftx = q[iitoi(leftii[1],n),iitoj(leftii[1],n),3]
for(i in 2:length(leftii)){
  leftx = c(leftx,q[iitoi(leftii[i],n),iitoj(leftii[i],n),3])
}

downx = q[iitoi(downii[1],n),iitoj(downii[1],n),4]
for(i in 2:length(downii)){
  downx = c(downx,q[iitoi(downii[i],n),iitoj(downii[i],n),4])
}

stayx = q[iitoi(1,n),iitoj(1,n),5]
for(i in 2:(n^2)){
  stayx = c(stayx,q[iitoi(i,n),iitoj(i,n),5])
}

G <- sparseMatrix(i = c(rightii,upii,leftii,downii,1:(n^2)),
                 j = c(rightjj,upjj,leftjj,downjj,1:(n^2)),
                 x = c(rightx,upx,leftx,downx,stayx))

# convert to a transition matrix
P <- G/(max(abs(G))+1) + sparseMatrix(i=1:(n^2),j=1:(n^2),x=rep(1,n^2))

nu <- rep(1/(n^2),n^2)
#image(matrix(nu,nrow=n,ncol=n))
for(i in 1:1000){
#  nu <- nu%*%(G + sparseMatrix(i=1:(n^2),j=1:(n^2),x=rep(1,n^2)))
  nu <- nu%*%P
  #image(matrix(nu,nrow=n,ncol=n))
  
}
image(t(matrix(nu,nrow=n,ncol=n)))
image(t(matrix(-log(nu),nrow=n,ncol=n)))