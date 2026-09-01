# ============================================================
# Pilot-only-100  vs  Proposed 20/80
#
# A. Pilot-only 100
#    처음부터 ~100개를 pilot sample로 뽑고, 그 100개만으로
#    beta_hat_pilot-100 을 추정  (sequential selection 없음)
#
# B. Proposed 20/80
#    먼저 ~20~25개의 pilot을 뽑고
#    -> sequential selection -> 100개
#    까지 확장한 뒤 beta_hat_proposed 를 추정
#
# 두 설계를 1000회 반복(MC)하여
#   MSE(beta_hat_pilot-100)  vs  MSE(beta_hat_proposed)
# 를 비교한다.
#
# 구현 방식:
#   원래의 run_one_sim(iter, n_pil, seed) 함수는 이미
#     - pilot_beta_hat / pilot_tau2_hat / pilot_phi_hat  : pilot만으로 추정
#     - beta_hat / tau2_hat / phi_hat                    : pilot -> sequential
#                                                           selection -> n_N개
#                                                           까지 확장 후 추정
#   을 모두 반환한다.
#
#   따라서:
#     n_pil = n_pil_A_target (~95) 로 호출하면
#         actual pilot size ~= n_pil + 5 ~= 100 (>= n_N)
#         => while(length(S) < n_N) 루프가 사실상 실행되지 않음
#         => pilot_beta_hat 이 곧 "Pilot-only-100" 추정치가 된다.
#
#     n_pil = n_pil_B_target (~20) 로 호출하면
#         actual pilot size ~= n_pil + 5 ~= 25 (< n_N)
#         => sequential selection이 25 -> n_N(100) 까지 실제로 진행됨
#         => beta_hat(final) 이 곧 "Proposed 20/80" 추정치가 된다.
#
#   같은 seed 로 두 번 호출하면 후보 모집단(coords, X)이 동일하게
#   생성되므로(첫 random draw들은 n_pil과 무관), 같은 population 위에서
#   두 설계를 비교하는 셈이 되어 더 공정한 비교가 된다.
# ============================================================

# ============================================================
# 0. Settings
# ============================================================

set.seed(42)

n_iter <- 100

N0   <- 10000
p    <- 3
n_N  <- 500     # 최종(full) 표본 크기

# ------------------------------------------------------------
# Pilot target sizes
#   actual pilot size ~= n_pil_target + 5   (원본 코드의 clone/pair 구조 때문)
# ------------------------------------------------------------

n_pil_A_target <- 495   # -> actual pilot ~= 500  (Pilot-only-500)
n_pil_B_target <- 100   # -> actual pilot ~= 205   (Proposed 200/300, pilot 단계)

beta_true <- c(
  2,
  0.5,
  1
)

sigma2   <- 1
nu       <- 0.5
rho_star <- 0.3

tau2_ratio <- 0.01

tau2 <- tau2_ratio * sigma2 /
  (1 - tau2_ratio)

c1 <- c(
  0.25,
  0.25
)

c2 <- c(
  0.75,
  0.75
)


# ============================================================
# 1. Matern correlation
# ============================================================

matern_corr <- function(u, phi, nu) {
  
  out <- matrix(
    1,
    nrow(u),
    ncol(u)
  )
  
  idx <- u > 1e-10
  
  a <- sqrt(2 * nu) *
    u[idx] / phi
  
  out[idx] <-
    (1 /
       (2^(nu - 1) * gamma(nu))) *
    a^nu *
    besselK(a, nu)
  
  out
}


# ============================================================
# 2. Find true phi
# ============================================================

find_phi <- function(
    rho_star,
    nu,
    target = 0.05
) {
  
  f <- function(phi) {
    
    matern_corr(
      matrix(rho_star),
      phi,
      nu
    ) - target
  }
  
  uniroot(
    f,
    c(1e-4, 10)
  )$root
}


phi_true <- find_phi(
  rho_star,
  nu
)


# ============================================================
# 3. Pairwise distance
# ============================================================

pairwise_dist <- function(A, B) {
  
  sqrt(
    outer(
      A[, 1],
      B[, 1],
      "-"
    )^2 +
      outer(
        A[, 2],
        B[, 2],
        "-"
      )^2
  )
}


# ============================================================
# 4. Covariance matrix
# ============================================================

build_cov <- function(
    idx1,
    idx2,
    phi,
    nu,
    sigma2,
    tau2,
    coords,
    same_set = FALSE
) {
  
  D <- pairwise_dist(
    coords[idx1, , drop = FALSE],
    coords[idx2, , drop = FALSE]
  )
  
  K <- sigma2 *
    matern_corr(
      D,
      phi,
      nu
    )
  
  if (same_set) {
    
    diag(K) <-
      diag(K) + tau2
  }
  
  K
}


# ============================================================
# 5. Maximin pilot selection
# ============================================================

select_pilot_maximin <- function(
    n_pil,
    N,
    coords
) {
  
  Spil <- integer(0)
  
  # First point
  Spil <- c(
    Spil,
    sample(
      seq_len(N),
      1
    )
  )
  
  while (
    length(Spil) < n_pil
  ) {
    
    R <- setdiff(
      seq_len(N),
      Spil
    )
    
    D_SR <- pairwise_dist(
      coords[R, , drop = FALSE],
      coords[Spil, , drop = FALSE]
    )
    
    min_dist <- apply(
      D_SR,
      1,
      min
    )
    
    i_star <- R[
      which.max(min_dist)
    ]
    
    Spil <- c(
      Spil,
      i_star
    )
  }
  
  Spil
}


# ============================================================
# 6. One complete simulation (given a target pilot size n_pil)
#
#    원본 코드 그대로 보존.
#    - n_pil이 커서 actual pilot >= n_N 이면 sequential selection은
#      실행되지 않고 final == pilot 상태가 된다.
#    - n_pil이 작으면 sequential selection이 실제로 n_N까지 확장한다.
# ============================================================

run_one_sim <- function(
    iter,
    n_pil,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    
    set.seed(seed)
  }
  
  
  # ==========================================================
  # STEP 1. Candidate locations + clone points
  # ==========================================================
  
  N <- N0
  
  coords <- matrix(
    runif(2 * N),
    N,
    2
  )
  
  # Number of close pairs
  n_pairs <- 5
  
  # Original base points
  pair_base_idx <- sample(
    N,
    n_pairs
  )
  
  # Clone points
  pair_coords <-
    coords[pair_base_idx, ] +
    matrix(
      runif(
        2 * n_pairs,
        -0.001,
        0.001
      ),
      n_pairs,
      2
    )
  
  coords <- rbind(
    coords,
    pair_coords
  )
  
  N <- N + n_pairs
  
  
  # ==========================================================
  # STEP 1-2. Covariates
  # ==========================================================
  
  dist_c1 <- sqrt(
    rowSums(
      (
        coords -
          matrix(
            c1,
            N,
            2,
            byrow = TRUE
          )
      )^2
    )
  )
  
  dist_c2 <- sqrt(
    rowSums(
      (
        coords -
          matrix(
            c2,
            N,
            2,
            byrow = TRUE
          )
      )^2
    )
  )
  
  X <- cbind(
    1,
    dist_c1,
    dist_c2
  )
  
  
  # ==========================================================
  # STEP 2. Pilot sampling
  # ==========================================================
  
  close_pair_idx <- c(
    pair_base_idx,
    (N - n_pairs + 1):N
  )
  
  
  # Preserve original algorithm
  n_pil_maximin <-
    n_pil - n_pairs
  
  
  # Safety check
  if (n_pil_maximin < 1) {
    
    stop(
      "n_pil is too small. ",
      "n_pil - n_pairs must be >= 1."
    )
  }
  
  
  Spil_maximin <-
    select_pilot_maximin(
      n_pil_maximin,
      N,
      coords
    )
  
  
  Spil <- unique(
    c(
      Spil_maximin,
      close_pair_idx
    )
  )
  
  
  # ==========================================================
  # STEP 3. Pilot response
  # ==========================================================
  
  Sigma_pil <- build_cov(
    Spil,
    Spil,
    phi_true,
    nu,
    sigma2,
    tau2,
    coords,
    same_set = TRUE
  )
  
  L_pil <- chol(
    Sigma_pil
  )
  
  Y_pil <- as.vector(
    X[Spil, ] %*%
      beta_true +
      t(L_pil) %*%
      rnorm(
        length(Spil)
      )
  )
  
  
  # ==========================================================
  # STEP 4. Pilot theta estimation
  # ==========================================================
  
  neg_profile_loglik <- function(
    par,
    idx,
    Xs,
    Ys
  ) {
    
    phi <- exp(
      par[1]
    )
    
    tau2_ <- exp(
      par[2]
    )
    
    
    Sigma <- build_cov(
      idx,
      idx,
      phi,
      nu,
      sigma2,
      tau2_,
      coords,
      same_set = TRUE
    )
    
    
    R <- tryCatch(
      chol(Sigma),
      error = function(e)
        NULL
    )
    
    
    if (is.null(R)) {
      
      return(1e20)
    }
    
    
    logdet <-
      2 * sum(
        log(
          diag(R)
        )
      )
    
    
    Sigma_inv <-
      chol2inv(R)
    
    
    XtSiX <-
      t(Xs) %*%
      Sigma_inv %*%
      Xs
    
    XtSiY <-
      t(Xs) %*%
      Sigma_inv %*%
      Ys
    
    
    beta_hat <- tryCatch(
      
      solve(
        XtSiX,
        XtSiY
      ),
      
      error = function(e)
        NULL
    )
    
    
    if (is.null(beta_hat)) {
      
      return(1e20)
    }
    
    
    resid <-
      Ys -
      Xs %*%
      beta_hat
    
    
    quad <- as.numeric(
      t(resid) %*%
        Sigma_inv %*%
        resid
    )
    
    
    0.5 *
      (
        logdet +
          quad
      )
  }
  
  
  # Pilot optimization
  opt <- optim(
    par = c(
      log(phi_true),
      log(tau2)
    ),
    
    fn =
      neg_profile_loglik,
    
    idx =
      Spil,
    
    Xs =
      X[Spil, ],
    
    Ys =
      Y_pil,
    
    method =
      "BFGS"
  )
  
  
  phi_hat <-
    exp(
      opt$par[1]
    )
  
  tau2_hat <-
    exp(
      opt$par[2]
    )
  
  
  # ==========================================================
  # STEP 4-2. Pilot beta estimation
  #
  # This is the pilot-stage beta estimator.
  # ==========================================================
  
  Sigma_pil_hat <-
    build_cov(
      Spil,
      Spil,
      phi_hat,
      nu,
      sigma2,
      tau2_hat,
      coords,
      same_set = TRUE
    )
  
  
  Sigma_pil_hat_inv <-
    chol2inv(
      chol(
        Sigma_pil_hat
      )
    )
  
  
  X_pil <-
    X[Spil, ,
      drop = FALSE
    ]
  
  
  beta_pil_hat <- solve(
    
    t(X_pil) %*%
      Sigma_pil_hat_inv %*%
      X_pil,
    
    t(X_pil) %*%
      Sigma_pil_hat_inv %*%
      Y_pil
  )
  
  
  # ==========================================================
  # STEP 5. Sequential selection
  #   (n_pil이 커서 actual pilot >= n_N 이면 이 루프는 실행되지 않는다)
  # ==========================================================
  
  S <- Spil
  
  R_idx <- setdiff(
    seq_len(N),
    S
  )
  
  
  total_var <-
    sigma2 +
    tau2_hat
  
  
  while (
    length(S) < n_N
  ) {
    
    
    Sigma_S <-
      build_cov(
        S,
        S,
        phi_hat,
        nu,
        sigma2,
        tau2_hat,
        coords,
        same_set = TRUE
      )
    
    
    Sigma_S_inv <-
      chol2inv(
        chol(
          Sigma_S
        )
      )
    
    
    X_S <-
      X[S, ,
        drop = FALSE
      ]
    
    
    Ibeta <-
      t(X_S) %*%
      Sigma_S_inv %*%
      X_S
    
    
    Ibeta_inv <-
      solve(
        Ibeta
      )
    
    
    Sigma_SR <-
      build_cov(
        S,
        R_idx,
        phi_hat,
        nu,
        sigma2,
        tau2_hat,
        coords,
        same_set = FALSE
      )
    
    
    tmp <-
      Sigma_S_inv %*%
      Sigma_SR
    
    
    v_i <-
      total_var -
      colSums(
        Sigma_SR * tmp
      )
    
    
    r_i <-
      t(
        X[R_idx, ,
          drop = FALSE
        ]
      ) -
      t(X_S) %*%
      tmp
    
    
    Ibinv_r <-
      Ibeta_inv %*%
      r_i
    
    
    Ibinv2_r <-
      Ibeta_inv %*%
      Ibinv_r
    
    
    numer <-
      colSums(
        r_i * Ibinv2_r
      )
    
    
    denom <-
      v_i +
      colSums(
        r_i * Ibinv_r
      )
    
    
    Delta_i <-
      numer /
      denom
    
    
    i_star <-
      R_idx[
        which.max(
          Delta_i
        )
      ]
    
    
    S <-
      c(
        S,
        i_star
      )
    
    
    R_idx <-
      setdiff(
        R_idx,
        i_star
      )
  }
  
  
  # ==========================================================
  # STEP 6. Generate final responses conditionally
  #   (S_new이 비어있으면, 즉 sequential selection이 실행되지
  #    않았으면 이 블록은 사실상 아무 것도 하지 않는다)
  # ==========================================================
  
  S_new <-
    setdiff(
      S,
      Spil
    )
  
  
  if (length(S_new) > 0) {
    
    Sigma_new_new <-
      build_cov(
        S_new,
        S_new,
        phi_true,
        nu,
        sigma2,
        tau2,
        coords,
        same_set = TRUE
      )
    
    
    Sigma_new_pil <-
      build_cov(
        S_new,
        Spil,
        phi_true,
        nu,
        sigma2,
        tau2,
        coords,
        same_set = FALSE
      )
    
    
    Sigma_pil_inv <-
      solve(
        Sigma_pil
      )
    
    
    cond_mean <-
      X[S_new, ] %*%
      beta_true +
      
      Sigma_new_pil %*%
      Sigma_pil_inv %*%
      (
        Y_pil -
          X[Spil, ] %*%
          beta_true
      )
    
    
    cond_cov <-
      Sigma_new_new -
      
      Sigma_new_pil %*%
      Sigma_pil_inv %*%
      t(
        Sigma_new_pil
      )
    
    
    cond_cov <-
      (
        cond_cov +
          t(cond_cov)
      ) / 2
    
    
    # Numerical stabilization
    eig <- eigen(
      cond_cov,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    
    if (
      min(eig) <= 1e-10
    ) {
      
      cond_cov <-
        cond_cov +
        
        diag(
          1e-8 -
            min(eig),
          nrow(cond_cov)
        )
    }
    
    
    L_new <-
      chol(
        cond_cov
      )
    
    
    Y_new <-
      as.vector(
        cond_mean
      ) +
      
      as.vector(
        t(L_new) %*%
          rnorm(
            length(S_new)
          )
      )
    
  } else {
    
    Y_new <- numeric(0)
  }
  
  
  S_all <-
    c(
      Spil,
      S_new
    )
  
  
  Y_all <-
    c(
      Y_pil,
      Y_new
    )
  
  
  # ==========================================================
  # STEP 7. Final theta estimation
  #   (S_new이 비어있으면 S_all == Spil 이므로 이 값은
  #    pilot 추정치와 사실상 동일하다)
  # ==========================================================
  
  opt_fin <- optim(
    
    par = c(
      log(phi_true),
      log(tau2)
    ),
    
    fn =
      neg_profile_loglik,
    
    idx =
      S_all,
    
    Xs =
      X[S_all, ],
    
    Ys =
      Y_all,
    
    method =
      "BFGS"
  )
  
  
  phi_fin <-
    exp(
      opt_fin$par[1]
    )
  
  
  tau2_fin <-
    exp(
      opt_fin$par[2]
    )
  
  
  # ==========================================================
  # STEP 8. Final GLS beta
  # ==========================================================
  
  Sigma_fin <-
    build_cov(
      S_all,
      S_all,
      phi_fin,
      nu,
      sigma2,
      tau2_fin,
      coords,
      same_set = TRUE
    )
  
  
  Sigma_fin_inv <-
    chol2inv(
      chol(
        Sigma_fin
      )
    )
  
  
  X_fin <-
    X[S_all, ,
      drop = FALSE
    ]
  
  
  beta_hat <- solve(
    
    t(X_fin) %*%
      Sigma_fin_inv %*%
      X_fin,
    
    t(X_fin) %*%
      Sigma_fin_inv %*%
      Y_all
  )
  
  
  # ==========================================================
  # STEP 9. Return results
  # ==========================================================
  
  list(
    
    # --------------------------
    # Final estimators (pilot -> sequential selection 이후)
    # --------------------------
    
    beta_hat =
      as.vector(
        beta_hat
      ),
    
    tau2_hat =
      tau2_fin,
    
    phi_hat =
      phi_fin,
    
    
    # --------------------------
    # Pilot-only estimators
    # --------------------------
    
    pilot_beta_hat =
      as.vector(
        beta_pil_hat
      ),
    
    pilot_tau2_hat =
      tau2_hat,
    
    pilot_phi_hat =
      phi_hat,
    
    
    # --------------------------
    # Sample sizes
    # --------------------------
    
    pilot_size =
      length(Spil),
    
    final_size =
      length(S_all)
  )
}


# ============================================================
# 7. Comparison wrapper: A(Pilot-only-100) vs B(Proposed 20/80)
#
#    동일한 seed로 run_one_sim을 두 번 호출한다.
#    -> coords/population 생성 부분(STEP 1, 1-2)은 n_pil과 무관하게
#       완전히 동일하게 재현되므로, 같은 후보 모집단 위에서
#       두 설계를 비교하는 셈이 된다.
# ============================================================

run_comparison_sim <- function(
    iter,
    seed,
    n_pil_A = n_pil_A_target,
    n_pil_B = n_pil_B_target
) {
  
  # ---- A. Pilot-only-100 ----
  # n_pil_A(~95)를 주면 actual pilot ~= 100 >= n_N 이 되어
  # sequential selection이 실행되지 않으므로
  # pilot_beta_hat 이 그대로 "pilot-only-100" 추정치가 된다.
  
  res_A <- run_one_sim(
    iter = iter,
    n_pil = n_pil_A,
    seed = seed
  )
  
  
  # ---- B. Proposed 20/80 ----
  # n_pil_B(~20)를 주면 actual pilot ~= 25 < n_N 이 되어
  # sequential selection이 25 -> n_N(100) 까지 실제로 진행되고
  # beta_hat(final) 이 "proposed 20/80" 추정치가 된다.
  
  res_B <- run_one_sim(
    iter = iter,
    n_pil = n_pil_B,
    seed = seed
  )
  
  
  list(
    
    # A: pilot-only-100
    beta_A = res_A$pilot_beta_hat,
    tau2_A = res_A$pilot_tau2_hat,
    phi_A  = res_A$pilot_phi_hat,
    size_A = res_A$pilot_size,
    
    # B: proposed 20/80
    beta_B = res_B$beta_hat,
    tau2_B = res_B$tau2_hat,
    phi_B  = res_B$phi_hat,
    pilot_size_B = res_B$pilot_size,
    final_size_B = res_B$final_size
  )
}


# ============================================================
# 8. Parallel Monte Carlo simulation
# ============================================================

library(parallel)


n_cores_total <- detectCores()

# SLURM이 이 job에 실제로 할당해준 코어 수를 우선 사용.
# (detectCores()는 노드 전체 코어 수를 반환하므로, sbatch에서 요청한
#  --cpus-per-task보다 훨씬 많은 프로세스를 띄워 다른 job과 충돌/과부하가
#  날 수 있음. SLURM_CPUS_PER_TASK 환경변수가 있으면 그 값을 우선 사용.)
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)

n_cores <- if (!is.na(slurm_cpus)) {
  as.integer(slurm_cpus)
} else {
  max(1, n_cores_total - 1)
}


cat("\n")
cat("============================================================\n")
cat("Pilot-only-100 vs Proposed 20/80 : Parallel Monte Carlo\n")
cat("============================================================\n")

cat(sprintf("Total CPU cores detected : %d\n", n_cores_total))
cat(sprintf("Cores used               : %d\n", n_cores))
cat(sprintf("Number of iterations     : %d\n", n_iter))
cat(sprintf("n_pil_A_target (~actual) : %d (~%d)\n", n_pil_A_target, n_pil_A_target + 5))
cat(sprintf("n_pil_B_target (~actual) : %d (~%d)\n", n_pil_B_target, n_pil_B_target + 5))
cat(sprintf("n_N (final size)         : %d\n", n_N))
cat("============================================================\n\n")


cl <- makeCluster(
  n_cores,
  type = "PSOCK"
)


tryCatch({
  
  clusterExport(
    cl,
    varlist = c(
      "N0",
      "p",
      "n_N",
      "n_pil_A_target",
      "n_pil_B_target",
      "beta_true",
      "sigma2",
      "nu",
      "rho_star",
      "tau2_ratio",
      "tau2",
      "c1",
      "c2",
      "phi_true",
      
      "matern_corr",
      "find_phi",
      "pairwise_dist",
      "build_cov",
      "select_pilot_maximin",
      "run_one_sim",
      "run_comparison_sim"
    ),
    envir = .GlobalEnv
  )
  
  
  cat("Simulation started...\n\n")
  
  start_time <- Sys.time()
  
  
  results <- parLapply(
    cl,
    seq_len(n_iter),
    
    function(b) {
      
      tryCatch(
        
        run_comparison_sim(
          iter = b,
          seed = 1000 + b
        ),
        
        error = function(e) {
          
          message(
            "Iteration ",
            b,
            " failed: ",
            conditionMessage(e)
          )
          
          NULL
        }
      )
    }
  )
  
  
  end_time <- Sys.time()
  
  stopCluster(cl)
  
  cat("\n")
  cat("Parallel simulation finished.\n")
  
  cat(
    sprintf(
      "Elapsed time: %.2f minutes\n",
      as.numeric(
        difftime(
          end_time,
          start_time,
          units = "mins"
        )
      )
    )
  )
  
}, error = function(e) {
  
  stopCluster(cl)
  
  stop(e)
})


# ============================================================
# 9. Successful iterations
# ============================================================

ok <-
  !vapply(
    results,
    is.null,
    logical(1)
  )


cat(
  sprintf(
    "Successful iterations: %d / %d\n",
    sum(ok),
    n_iter
  )
)


if (sum(ok) == 0) {
  
  stop("No successful simulation.")
}


results_ok <-
  results[ok]


# ============================================================
# 10. Extract beta estimates
# ============================================================

beta_A_mat <-
  do.call(
    rbind,
    lapply(
      results_ok,
      function(z) z$beta_A
    )
  )

colnames(beta_A_mat) <- c("beta0", "beta1", "beta2")


beta_B_mat <-
  do.call(
    rbind,
    lapply(
      results_ok,
      function(z) z$beta_B
    )
  )

colnames(beta_B_mat) <- c("beta0", "beta1", "beta2")


# ============================================================
# 11. Extract theta estimates
# ============================================================

tau2_A_vec <-
  vapply(results_ok, function(z) z$tau2_A, numeric(1))

tau2_B_vec <-
  vapply(results_ok, function(z) z$tau2_B, numeric(1))

phi_A_vec <-
  vapply(results_ok, function(z) z$phi_A, numeric(1))

phi_B_vec <-
  vapply(results_ok, function(z) z$phi_B, numeric(1))


# ============================================================
# 12. Helper function
# ============================================================

calc_summary <- function(
    estimates,
    truth
) {
  
  bias <- mean(estimates - truth)
  
  variance <- var(estimates)
  
  mse <- mean((estimates - truth)^2)
  
  rmse <- sqrt(mse)
  
  c(
    mean_estimate = mean(estimates),
    bias = bias,
    variance = variance,
    MSE = mse,
    RMSE = rmse
  )
}


# ============================================================
# 13. Regression coefficient comparison
# ============================================================

beta_summary_list <- list()


for (j in seq_along(beta_true)) {
  
  A_result <- calc_summary(beta_A_mat[, j], beta_true[j])
  B_result <- calc_summary(beta_B_mat[, j], beta_true[j])
  
  beta_summary_list[[j]] <-
    data.frame(
      
      parameter = paste0("beta", j - 1),
      true = beta_true[j],
      
      pilot100_mean     = A_result["mean_estimate"],
      pilot100_bias     = A_result["bias"],
      pilot100_variance = A_result["variance"],
      pilot100_MSE      = A_result["MSE"],
      pilot100_RMSE     = A_result["RMSE"],
      
      proposed_mean     = B_result["mean_estimate"],
      proposed_bias     = B_result["bias"],
      proposed_variance = B_result["variance"],
      proposed_MSE      = B_result["MSE"],
      proposed_RMSE     = B_result["RMSE"]
    )
}


beta_comparison <-
  do.call(rbind, beta_summary_list)

rownames(beta_comparison) <- NULL


# ============================================================
# 14. Tau2 comparison
# ============================================================

tau2_A_summary <- calc_summary(tau2_A_vec, tau2)
tau2_B_summary <- calc_summary(tau2_B_vec, tau2)

tau2_comparison <-
  data.frame(
    
    parameter = "tau2",
    true = tau2,
    
    pilot100_mean     = tau2_A_summary["mean_estimate"],
    pilot100_bias     = tau2_A_summary["bias"],
    pilot100_variance = tau2_A_summary["variance"],
    pilot100_MSE      = tau2_A_summary["MSE"],
    pilot100_RMSE     = tau2_A_summary["RMSE"],
    
    proposed_mean     = tau2_B_summary["mean_estimate"],
    proposed_bias     = tau2_B_summary["bias"],
    proposed_variance = tau2_B_summary["variance"],
    proposed_MSE      = tau2_B_summary["MSE"],
    proposed_RMSE     = tau2_B_summary["RMSE"]
  )


# ============================================================
# 15. Phi comparison
# ============================================================

phi_A_summary <- calc_summary(phi_A_vec, phi_true)
phi_B_summary <- calc_summary(phi_B_vec, phi_true)

phi_comparison <-
  data.frame(
    
    parameter = "phi",
    true = phi_true,
    
    pilot100_mean     = phi_A_summary["mean_estimate"],
    pilot100_bias     = phi_A_summary["bias"],
    pilot100_variance = phi_A_summary["variance"],
    pilot100_MSE      = phi_A_summary["MSE"],
    pilot100_RMSE     = phi_A_summary["RMSE"],
    
    proposed_mean     = phi_B_summary["mean_estimate"],
    proposed_bias     = phi_B_summary["bias"],
    proposed_variance = phi_B_summary["variance"],
    proposed_MSE      = phi_B_summary["MSE"],
    proposed_RMSE     = phi_B_summary["RMSE"]
  )


# ============================================================
# 16. Combined summary
# ============================================================

comparison_summary <-
  rbind(
    beta_comparison,
    tau2_comparison,
    phi_comparison
  )


# ============================================================
# 17. Print results
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PILOT-ONLY-100 vs PROPOSED 20/80\n")
cat("============================================================\n\n")


comparison_summary_print <- comparison_summary

numeric_cols <- vapply(comparison_summary_print, is.numeric, logical(1))

comparison_summary_print[numeric_cols] <-
  lapply(comparison_summary_print[numeric_cols], round, digits = 6)


print(comparison_summary_print, row.names = FALSE)


# ============================================================
# 18. KEY RESULT: beta1
# ============================================================

cat("\n")
cat("============================================================\n")
cat("KEY RESULT: beta1\n")
cat("============================================================\n")

beta1_result <-
  comparison_summary[
    comparison_summary$parameter == "beta1", ,
    drop = FALSE
  ]

beta1_numeric <- vapply(beta1_result, is.numeric, logical(1))

beta1_result[beta1_numeric] <-
  lapply(beta1_result[beta1_numeric], round, digits = 6)

print(beta1_result, row.names = FALSE)


# ============================================================
# 19. KEY RESULT: tau2
# ============================================================

cat("\n")
cat("============================================================\n")
cat("KEY RESULT: tau2\n")
cat("============================================================\n")

tau2_result <-
  comparison_summary[
    comparison_summary$parameter == "tau2", ,
    drop = FALSE
  ]

tau2_numeric <- vapply(tau2_result, is.numeric, logical(1))

tau2_result[tau2_numeric] <-
  lapply(tau2_result[tau2_numeric], round, digits = 6)

print(tau2_result, row.names = FALSE)


# ============================================================
# 20. MSE improvement (Proposed vs Pilot-only-100)
# ============================================================

beta_mse_improvement <-
  data.frame(
    
    parameter = beta_comparison$parameter,
    
    pilot100_MSE = beta_comparison$pilot100_MSE,
    proposed_MSE = beta_comparison$proposed_MSE,
    
    MSE_reduction =
      beta_comparison$pilot100_MSE -
      beta_comparison$proposed_MSE,
    
    relative_MSE_reduction =
      1 - (beta_comparison$proposed_MSE / beta_comparison$pilot100_MSE)
  )


cat("\n")
cat("============================================================\n")
cat("MSE IMPROVEMENT OF PROPOSED (20/80) vs PILOT-ONLY-100\n")
cat("============================================================\n\n")


beta_mse_improvement_print <- beta_mse_improvement

numeric_cols <- vapply(beta_mse_improvement_print, is.numeric, logical(1))

beta_mse_improvement_print[numeric_cols] <-
  lapply(beta_mse_improvement_print[numeric_cols], round, digits = 6)


print(beta_mse_improvement_print, row.names = FALSE)


# ============================================================
# 21. Save results
# ============================================================

write.csv(
  beta_A_mat,
  "pilot_only_100_beta_estimates.csv",
  row.names = FALSE
)

write.csv(
  beta_B_mat,
  "proposed_20_80_beta_estimates.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    pilot100_tau2 = tau2_A_vec,
    proposed_tau2 = tau2_B_vec,
    pilot100_phi  = phi_A_vec,
    proposed_phi  = phi_B_vec
  ),
  "pilot100_vs_proposed_theta_estimates.csv",
  row.names = FALSE
)

write.csv(
  beta_comparison,
  "pilot100_vs_proposed_beta_summary.csv",
  row.names = FALSE
)

write.csv(
  tau2_comparison,
  "pilot100_vs_proposed_tau2_summary.csv",
  row.names = FALSE
)

write.csv(
  phi_comparison,
  "pilot100_vs_proposed_phi_summary.csv",
  row.names = FALSE
)

write.csv(
  comparison_summary,
  "pilot100_vs_proposed_complete_summary.csv",
  row.names = FALSE
)

write.csv(
  beta_mse_improvement,
  "proposed_beta_MSE_improvement.csv",
  row.names = FALSE
)


# ============================================================
# 22. Save complete RDS
# ============================================================

saveRDS(
  
  list(
    
    settings = list(
      n_iter = n_iter,
      N0 = N0,
      n_N = n_N,
      n_pil_A_target = n_pil_A_target,
      n_pil_B_target = n_pil_B_target,
      beta_true = beta_true,
      sigma2 = sigma2,
      nu = nu,
      rho_star = rho_star,
      tau2 = tau2,
      phi_true = phi_true,
      c1 = c1,
      c2 = c2
    ),
    
    beta_A = beta_A_mat,
    beta_B = beta_B_mat,
    
    tau2_A = tau2_A_vec,
    tau2_B = tau2_B_vec,
    
    phi_A = phi_A_vec,
    phi_B = phi_B_vec,
    
    beta_comparison = beta_comparison,
    tau2_comparison = tau2_comparison,
    phi_comparison = phi_comparison,
    
    comparison_summary = comparison_summary,
    beta_mse_improvement = beta_mse_improvement
  ),
  
  "pilot_only_500_vs_proposed_200_300_results.rds"
)


cat("\n")
cat("============================================================\n")
cat("ALL DONE\n")
cat("============================================================\n")