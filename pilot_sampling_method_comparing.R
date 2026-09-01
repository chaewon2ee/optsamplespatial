# ============================================================
# Monte Carlo: pilot sampling comparison (Bug Fixed Version)
# ============================================================
rm(list=ls()) # 완벽한 초기화를 위해 맨 처음에 실행

set.seed(42)

n_iter <- 1000
N0 <- 10000; p <- 3; n_N <- 100; n_pil <- 20
beta_true <- c(2,0.5,1)
sigma2 <- 1; nu <- 0.5; rho_star <- 0.3
tau2_ratio <- 0.01
tau2 <- tau2_ratio*sigma2/(1-tau2_ratio)
c1 <- c(0.25,0.25); c2 <- c(0.75,0.75)

matern_corr <- function(u,phi,nu){
  out <- matrix(1,nrow(u),ncol(u))
  idx <- u>1e-10
  a <- sqrt(2*nu)*u[idx]/phi
  out[idx] <- a^nu*besselK(a,nu)/(2^(nu-1)*gamma(nu))
  out
}

find_phi <- function(rho_star,nu,target=0.05){
  uniroot(function(phi)matern_corr(matrix(rho_star),phi,nu)-target,
          c(1e-4,10))$root
}

phi_true <- find_phi(rho_star,nu)

pairwise_dist <- function(A,B)
  sqrt(outer(A[,1],B[,1],"-")^2+outer(A[,2],B[,2],"-")^2)

build_cov <- function(idx1,idx2,phi,nu,sigma2,tau2,coords,same_set=FALSE){
  K <- sigma2*matern_corr(
    pairwise_dist(coords[idx1,,drop=FALSE],coords[idx2,,drop=FALSE]),
    phi,nu
  )
  if(same_set) diag(K) <- diag(K)+tau2
  K
}

# ------------------------------------------------------------
# Pilot sampling
# ------------------------------------------------------------

select_pilot_maximin <- function(n,N,coords,exclude=integer(0)){
  A <- setdiff(seq_len(N),exclude)
  S <- sample(A,1)
  while(length(S)<n){
    R <- setdiff(A,S)
    D <- pairwise_dist(coords[R,,drop=FALSE],coords[S,,drop=FALSE])
    S <- c(S,R[which.max(apply(D,1,min))])
  }
  S
}

select_pilot_grid <- function(n,N,coords,exclude=integer(0)){
  A <- setdiff(seq_len(N),exclude)
  g <- expand.grid(x=(1:5-.5)/5,y=(1:3-.5)/3)
  S <- integer(0)
  for(j in seq_len(min(n,nrow(g)))){
    R <- setdiff(A,S)
    if(!length(R)) break
    d <- (coords[R,1]-g$x[j])^2+(coords[R,2]-g$y[j])^2
    S <- c(S,R[which.min(d)])
  }
  S
}

select_pilot_stratified <- function(n,N,coords,exclude=integer(0)){
  A <- setdiff(seq_len(N),exclude)
  S <- integer(0)
  for(i in 1:5) for(j in 1:3){
    if(length(S)>=n) break
    R <- A[
      coords[A,1]>=(i-1)/5 & coords[A,1]<i/5 &
        coords[A,2]>=(j-1)/3 & coords[A,2]<j/3
    ]
    if(i==5)
      R <- A[
        coords[A,1]>=(i-1)/5 & coords[A,1]<=1 &
          coords[A,2]>=(j-1)/3 & coords[A,2]<j/3
      ]
    R <- setdiff(R,S)
    if(length(R)) S <- c(S,sample(R,1))
  }
  S
}

select_pilot_spatial_balanced <- function(n,N,coords,exclude=integer(0)){
  A <- setdiff(seq_len(N),exclude)
  S <- sample(A,1)
  while(length(S)<n){
    R <- setdiff(A,S)
    D <- pairwise_dist(coords[R,,drop=FALSE],coords[S,,drop=FALSE])
    w <- apply(D,1,min)
    w <- pmax(w,1e-12)
    S <- c(S,sample(R,1,prob=w))
  }
  S
}

select_pilot <- function(method,n,N,coords,close_pair_idx){
  S <- switch(
    method,
    grid=select_pilot_grid(n,N,coords,close_pair_idx),
    stratified=select_pilot_stratified(n,N,coords,close_pair_idx),
    spatial_balanced=select_pilot_spatial_balanced(n,N,coords,close_pair_idx),
    maximin=select_pilot_maximin(n,N,coords,close_pair_idx)
  )
  unique(c(S,close_pair_idx))
}

# ------------------------------------------------------------
# One simulation
# ------------------------------------------------------------

run_one_sim <- function(iter,seed=NULL,method="maximin",stage1=FALSE){
  
  if(!is.null(seed)) set.seed(seed)
  
  N <- N0
  coords <- matrix(runif(2*N),N,2)
  
  n_pairs <- 5
  pair_base_idx <- sample(seq_len(N),n_pairs)
  
  pair_coords <- coords[pair_base_idx,,drop=FALSE] +
    matrix(runif(2*n_pairs,-0.001,0.001),n_pairs,2)
  
  coords <- rbind(coords,pair_coords)
  N <- N+n_pairs
  
  X <- cbind(
    1,
    sqrt(rowSums((coords-matrix(c1,N,2,byrow=TRUE))^2)),
    sqrt(rowSums((coords-matrix(c2,N,2,byrow=TRUE))^2))
  )
  
  close_pair_idx <- c(pair_base_idx,(N-n_pairs+1):N)
  
  Spil <- select_pilot(
    method,n_pil-n_pairs,N,coords,close_pair_idx
  )
  
  if(length(Spil)!=25)
    stop("Pilot size error: ",length(Spil))
  
  Sigma_pil <- build_cov(
    Spil,Spil,phi_true,nu,sigma2,tau2,coords,TRUE
  )
  
  L_pil <- chol(Sigma_pil)
  
  Y_pil <- as.vector(
    X[Spil,,drop=FALSE]%*%beta_true+
      t(L_pil)%*%rnorm(length(Spil))
  )
  
  neg_profile_loglik <- function(par,idx,Xs,Ys){
    phi <- exp(par[1])
    tau <- exp(par[2])
    
    Sigma <- build_cov(idx,idx,phi,nu,sigma2,tau,coords,TRUE)
    
    R <- tryCatch(chol(Sigma),error=function(e)NULL)
    if(is.null(R)) return(1e20)
    
    Si <- chol2inv(R)
    logdet <- 2*sum(log(diag(R)))
    
    XtSiX <- crossprod(Xs,Si%*%Xs)
    XtSiY <- crossprod(Xs,Si%*%Ys)
    
    beta <- tryCatch(
      solve(XtSiX,XtSiY),
      error=function(e)NULL
    )
    
    if(is.null(beta)) return(1e20)
    
    e <- Ys-Xs%*%beta
    quad <- as.numeric(crossprod(e,Si%*%e))
    
    0.5*(logdet+quad)
  }
  
  opt <- optim(
    c(log(phi_true),log(tau2)),
    neg_profile_loglik,
    idx=Spil,
    Xs=X[Spil,,drop=FALSE],
    Ys=Y_pil,
    method="BFGS"
  )
  
  phi_hat <- exp(opt$par[1])
  tau2_hat <- exp(opt$par[2])
  
  if(stage1)
    return(list(
      pilot_phi_hat=phi_hat,
      pilot_tau2_hat=tau2_hat,
      pilot_size=length(Spil)
    ))
  
  # ----------------------------------------------------------
  # Sequential selection
  # ----------------------------------------------------------
  
  S <- Spil
  R_idx <- setdiff(seq_len(N),S)
  total_var <- sigma2+tau2_hat
  
  while(length(S)<n_N){
    
    Sigma_S <- build_cov(
      S,S,phi_hat,nu,sigma2,tau2_hat,coords,TRUE
    )
    
    Sigma_S_inv <- chol2inv(chol(Sigma_S))
    X_S <- X[S,,drop=FALSE]
    
    Ibeta <- crossprod(X_S,Sigma_S_inv%*%X_S)
    Ibeta_inv <- solve(Ibeta)
    
    Sigma_SR <- build_cov(
      S,R_idx,phi_hat,nu,sigma2,tau2_hat,coords
    )
    
    tmp <- Sigma_S_inv%*%Sigma_SR
    
    v_i <- total_var-colSums(Sigma_SR*tmp)
    
    r_i <- t(X[R_idx,,drop=FALSE])-
      t(X_S)%*%tmp
    
    Ibinv_r <- Ibeta_inv%*%r_i
    Ibinv2_r <- Ibeta_inv%*%Ibinv_r
    
    numer <- colSums(r_i*Ibinv2_r)
    denom <- v_i+colSums(r_i*Ibinv_r)
    
    Delta_i <- numer/denom
    
    if(any(!is.finite(Delta_i)))
      stop("Non-finite Delta_i")
    
    i_star <- R_idx[which.max(Delta_i)]
    
    S <- c(S,i_star)
    R_idx <- setdiff(R_idx,i_star)
  }
  
  # ----------------------------------------------------------
  # Conditional response generation
  # ----------------------------------------------------------
  
  S_new <- setdiff(S,Spil)
  
  Sigma_new_new <- build_cov(
    S_new,S_new,phi_true,nu,sigma2,tau2,coords,TRUE
  )
  
  Sigma_new_pil <- build_cov(
    S_new,Spil,phi_true,nu,sigma2,tau2,coords
  )
  
  Sigma_pil_inv <- chol2inv(chol(Sigma_pil))
  
  cond_mean <- X[S_new,,drop=FALSE]%*%beta_true+
    Sigma_new_pil%*%Sigma_pil_inv%*%
    (Y_pil-X[Spil,,drop=FALSE]%*%beta_true)
  
  cond_cov <- Sigma_new_new-
    Sigma_new_pil%*%Sigma_pil_inv%*%t(Sigma_new_pil)
  
  cond_cov <- (cond_cov+t(cond_cov))/2
  
  eig <- eigen(cond_cov,symmetric=TRUE,only.values=TRUE)$values
  
  if(min(eig)<=1e-10)
    cond_cov <- cond_cov+
    diag(1e-8-min(eig),nrow(cond_cov))
  
  Y_new <- as.vector(cond_mean)+
    as.vector(t(chol(cond_cov))%*%rnorm(length(S_new)))
  
  S_all <- c(Spil,S_new)
  Y_all <- c(Y_pil,Y_new)
  
  # ----------------------------------------------------------
  # Final estimation
  # ----------------------------------------------------------
  
  opt_fin <- optim(
    c(log(phi_true),log(tau2)),
    neg_profile_loglik,
    idx=S_all,
    Xs=X[S_all,,drop=FALSE],
    Ys=Y_all,
    method="BFGS"
  )
  
  phi_fin <- exp(opt_fin$par[1])
  tau2_fin <- exp(opt_fin$par[2])
  
  Sigma_fin <- build_cov(
    S_all,S_all,phi_fin,nu,sigma2,tau2_fin,coords,TRUE
  )
  
  Si <- chol2inv(chol(Sigma_fin))
  X_fin <- X[S_all,,drop=FALSE]
  
  beta_hat <- solve(
    crossprod(X_fin,Si%*%X_fin),
    crossprod(X_fin,Si%*%Y_all)
  )
  
  list(
    beta_hat=as.vector(beta_hat),
    tau2_hat=tau2_fin,
    phi_hat=phi_fin,
    pilot_tau2_hat=tau2_hat,
    pilot_phi_hat=phi_hat,
    pilot_size=length(Spil),
    final_size=length(S_all),
    convergence_pilot=opt$convergence,
    convergence_final=opt_fin$convergence
  )
}

# ------------------------------------------------------------
# 완벽하게 방어된 Summary 함수 (핵심 수정 사항)
# ------------------------------------------------------------

summarize <- function(z, method){
  # 에러 방지용 완벽한 타입 필터링
  if(!is.list(z)) return(NULL)
  valid_idx <- vapply(z, function(x) is.list(x) && !is.null(x), logical(1))
  z <- z[valid_idx]
  if(length(z) == 0) return(NULL)
  
  # [[ 연산자를 원천 차단하고 $ 로만 접근하여 get1index 에러 예방
  B <- do.call(rbind, lapply(z, function(x) x$beta_hat))
  tau <- vapply(z, function(x) x$tau2_hat, numeric(1))
  phi <- vapply(z, function(x) x$phi_hat, numeric(1))
  
  d <- sweep(B, 2, beta_true)
  M <- c(colMeans(d^2), mean((tau-tau2)^2), mean((phi-phi_true)^2))
  
  data.frame(
    method=method,
    parameter=c("beta0","beta1","beta2","tau2","phi"),
    true=c(beta_true,tau2,phi_true),
    mean_estimate=c(colMeans(B),mean(tau),mean(phi)),
    bias=c(colMeans(d),mean(tau-tau2),mean(phi-phi_true)),
    variance=c(apply(B,2,var),var(tau),var(phi)),
    MSE=M,
    RMSE=sqrt(M)
  )
}

summarize1 <- function(z, method){
  # 에러 방지용 완벽한 타입 필터링
  if(!is.list(z)) return(data.frame(method=method, successful=0, phi_MSE=NA, tau2_MSE=NA))
  valid_idx <- vapply(z, function(x) is.list(x) && !is.null(x), logical(1))
  z <- z[valid_idx]
  if(length(z) == 0) {
    return(data.frame(method=method, successful=0, phi_MSE=NA, tau2_MSE=NA))
  }
  
  # [[ 연산자를 원천 차단하고 $ 로만 접근
  tau <- vapply(z, function(x) x$pilot_tau2_hat, numeric(1))
  phi <- vapply(z, function(x) x$pilot_phi_hat, numeric(1))
  
  data.frame(
    method=method,
    successful=length(z),
    phi_MSE=mean((phi-phi_true)^2),
    tau2_MSE=mean((tau-tau2)^2)
  )
}

# ------------------------------------------------------------
# Parallel MC
# ------------------------------------------------------------
methods <- c("grid", "stratified", "spatial_balanced", "maximin")
n1 <- 1000
n2 <- 1000

cores <- max(1, parallel::detectCores(logical=TRUE)-1)
cl <- parallel::makeCluster(cores)

run_mc <- function(method,n,stage1,seed0){
  parallel::parLapplyLB(
    cl,seq_len(n),
    function(b){
      tryCatch(
        run_one_sim(b,seed0+b,method,stage1),
        error=function(e){
          message("ERROR | method=",method," | iter=",b," | ",conditionMessage(e))
          NULL
        }
      )
    }
  )
}

parallel::clusterExport(
  cl,
  c(
    "N0","p","n_N","n_pil",
    "beta_true","sigma2","nu","rho_star",
    "tau2_ratio","tau2","c1","c2",
    "phi_true",
    "matern_corr","pairwise_dist","build_cov",
    "select_pilot_maximin",
    "select_pilot_grid",
    "select_pilot_stratified",
    "select_pilot_spatial_balanced",
    "select_pilot",
    "run_one_sim"
  ),
  envir=environment()
)

parallel::clusterSetRNGStream(cl,42)

# ------------------------------------------------------------
# Stage 1
# ------------------------------------------------------------

s1 <- do.call(
  rbind,
  lapply(seq_along(methods),function(k){
    
    m <- methods[k]
    cat("Stage1:",m,"\n")
    
    summarize1(
      run_mc(m,n1,TRUE,100000*k),
      m
    )
  })
)

s1$score <- s1$phi_MSE/phi_true^2 + s1$tau2_MSE/tau2^2
s1 <- s1[order(s1$score),]
print(s1)

# ------------------------------------------------------------
# Stage 2
# ------------------------------------------------------------

top2 <- head(s1$method[s1$method!="maximin"], 2)
stage2_methods <- unique(c(top2,"maximin"))

s2 <- do.call(
  rbind,
  lapply(seq_along(stage2_methods),function(k){
    
    m <- stage2_methods[k]
    cat("Stage2:",m,"\n")
    
    summarize(
      run_mc(m,n2,FALSE,500000+10000*k),
      m
    )
  })
)

print(s2)

write.csv(s1, "Stage1_pilot_MSE.csv", row.names=FALSE)
write.csv(s2, "Stage2_full_MSE.csv", row.names=FALSE)

parallel::stopCluster(cl)