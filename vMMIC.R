vMMIC <- function(inputdata, max_iter, epsilon_conv){
  
  converged <- 0
  
  n<-length(inputdata[[1]])
  m<-dim(inputdata[[1]][[1]])[1]
  d0<-dim(inputdata[[1]][[1]])[2]
  d1<-dim(inputdata[[2]][[1]])[2]
  
  # specify hyperparameters of alpha, beta, gamma, a, b, c and d
  
  y_vec <- inputdata$y
  
  # sample initials alpha, beta, gamma, a, b, c and d
  
  E_delta <- matrix(0.5, nrow = m, ncol = n)
  E_U <- matrix(0.5, nrow = m, ncol = n)
  E_y <- matrix(0.5, nrow = n, ncol=1)

  V_beta_inv <- solve(Sigma_beta)
  m_beta <- solve(Sigma_beta)%*%mu_beta
  for (i in 1:n){
    w <- inputdata$w[[i]]
    Xi <- inputdata$X[[i]]
    Zi <- inputdata$Z[[i]]
    E_li <- matrix(E_delta[,i]*(1-w), ncol=1)
    E_li_prime <- matrix(E_delta[,i]*w, ncol=1)
    E_y[i,] <- pnorm(beta_true[1]+t(E_li)%*%Xi[,-1]%*%beta_true[-1] + 
                       gamma_true[1] + t(E_li_prime)%*%Zi[,-1]%*%gamma_true[-1])
    E_delta[, i] <- pnorm((1-w)*(b[1] + Xi[,-1]%*%b[-1]) + w*(d[1]+Zi[,-1]%*%d[-1]))
    temp_E_delta <- matrix(E_delta[,i], ncol=1)
    E_lili <- (E_li) %*% t(E_li) + diag(0.25, m)
    E_lili_prime <- (E_li) %*% t(E_li_prime)
    V_beta_inv <- V_beta_inv + t(Xi)%*%E_lili%*%Xi
    m_beta <- m_beta + E_y[i]*t(Xi)%*%E_li - t(Xi)%*%E_lili_prime%*%Zi%*%gamma_true
  }
  V_beta <- solve(V_beta_inv)
  m_beta <- V_beta %*% m_beta
  E_beta <- m_beta
  E_beta2 <- matrix(outer(E_beta, E_beta), nrow = length(E_beta), ncol = length(E_beta)) + V_beta
  
  V_gamma_inv <- solve(Sigma_gamma)
  m_gamma <- solve(Sigma_gamma)%*%mu_gamma
  for (i in 1: n){
    w <- inputdata$w[[i]]
    Xi <- inputdata$X[[i]]
    Zi <- inputdata$Z[[i]]
    E_li <- matrix(E_delta[,i]*(1-w), ncol=1)
    E_li_prime <- matrix(E_delta[,i]*w, ncol=1)
    E_li_prime_li <- (E_li_prime) %*% t(E_li)
    E_li_prime_li_prime <- E_li_prime %*% t(E_li_prime) + diag(0.25, m)
    V_gamma_inv <- V_gamma_inv + t(Zi)%*%E_li_prime_li_prime%*%Zi
    m_gamma <- m_gamma + E_y[i]*t(Zi)%*%E_li_prime - t(Zi)%*%E_li_prime_li%*%Xi%*%E_beta
  }
  V_gamma <- solve(V_gamma_inv)
  m_gamma <- V_gamma %*% m_gamma
  E_gamma <- m_gamma
  E_gamma2 <- matrix(outer(E_gamma, E_gamma), nrow = length(E_gamma), ncol = length(E_gamma)) + V_gamma
  
  V_b_inv <- solve(Sigma_b)
  m_b <- solve(Sigma_b)%*%mu_b
  for (i in 1:n){
    w <- inputdata$w[[i]]
    Xi <- inputdata$X[[i]]
    Zi <- inputdata$Z[[i]]
    E_U[, i] <- (1-w)*(b[1]+Xi[,-1]%*%b[-1]) + w*(d[1]+Zi[,-1]%*%d[-1])
    V_b_inv <- V_b_inv +t((1-w)*Xi)%*%Xi
    m_b <- m_b + t(matrix((1-w)*E_U[, i], nrow = 1) %*% Xi)
  }
  V_b <- solve(V_b_inv)
  m_b <- V_b %*% m_b
  E_b <- m_b
  E_b2 <- matrix(outer(E_b, E_b), nrow = length(E_b), ncol = length(E_b)) + V_b
  
  V_d_inv <- solve(Sigma_d)
  m_d <- solve(Sigma_d)%*%mu_d
  for (i in 1:n){
    w <- inputdata$w[[i]]
    Xi <- inputdata$X[[i]]
    Zi <- inputdata$Z[[i]]
    E_U[, i] <- (1-w)*(E_b[1]+Xi[,-1]%*%E_b[-1]) +w*(d[1]+Zi[,-1]%*%d[-1])
    V_d_inv <- V_d_inv +t(w*Zi)%*%Zi
    m_d <- m_d + t(matrix(w*E_U[, i], nrow = 1) %*% Zi)
  }
  V_d <- solve(V_d_inv)
  m_d <- V_d %*% m_d
  E_d <- m_d
  E_d2 <- matrix(outer(E_d, E_d), nrow = length(E_d), ncol = length(E_d)) + V_d
  
  w <- inputdata$w
  X <- inputdata$X
  Z <- inputdata$Z
  
  E_y1 <- rep(NA, n)
  E_y0 <- rep(NA, n)
  m_y <- rep(NA, n)
  m_y <- sapply(1:n, function(i) E_beta[1]+t(E_delta[,i]*(1-w[[i]]))%*%X[[i]][,-1]%*%E_beta[-1] + 
                  E_gamma[1] + t(E_delta[,i]*w[[i]])%*%Z[[i]][,-1]%*%E_gamma[-1])
  m_y[m_y > 5] <- 5
  m_y[m_y < -5] <- -5
  E_y1 <- m_y + dnorm(-m_y)/(1-pnorm(-m_y))
  E_y0 <- m_y - dnorm(-m_y)/pnorm(-m_y)
  E_y[y_vec == 0] <- E_y0[y_vec ==0]
  E_y[y_vec == 1] <- E_y1[y_vec ==1]
  
  m_U <- sapply(1:n, function(i) E_b[1]+(1-w[[i]])*(X[[i]][,-1]%*%E_b[-1]) + 
                  E_d[1]+w[[i]]*(Z[[i]][,-1]%*%E_d[-1]))
  E_U <- m_U
  E_U2 <- m_U^2+1
  
  ad1_ini <- sapply(1:n, function(i) E_beta[1]+ E_gamma[1]+ E_delta[,i]*((as.numeric(1-w[[i]]) * X[[i]][,-1]) %*% E_beta[-1]+
                                                                           (as.numeric(w[[i]]) * Z[[i]][,-1]) %*% E_gamma[-1]))
  ad1 <- matrix(rep(colSums(ad1_ini), m), byrow = TRUE, nrow = m) - ad1_ini 
  ad1 <- (matrix(rep(E_y, m), byrow = TRUE, nrow = m)-ad1)*
    sapply(1:n, function(i) E_beta[1]+ (as.numeric(1-w[[i]]) * X[[i]][,-1]) %*% E_beta[-1] + 
             E_gamma[1] + (as.numeric(w[[i]]) * Z[[i]][,-1]) %*% E_gamma[-1])
  ad2 <- sapply(1:n, function(i) 
    -diag((as.numeric(1-w[[i]]) * X[[i]]) %*% E_beta2 %*% t(X[[i]]) +
            (as.numeric(w[[i]]) * Z[[i]]) %*% E_gamma2  %*% t(Z[[i]]))/2)
  A <- (1-pnorm(-m_U)) * exp(ad1) * exp(ad2)
  B <- pnorm(-m_U)
  AA <- A/(A+B)
  BB <- B/(A+B)
  AA[is.nan(AA)] <- 0.5
  BB[is.nan(BB)] <- 0.5
  E_delta <- AA
  V_delta <- AA*BB
  E_delta2 <- E_delta^2 + V_delta
  
  E_U2 <- (BB-AA)*(m_U^2+1)*pnorm(-m_U) + (AA-BB)*m_U*exp(-1/2*m_U^2)/sqrt(2*pi) +
    AA*(m_U^2 + 1)/(AA + (BB-AA*pnorm(-m_U)))
  p <- c()
  p <- sapply(1:n, function(i) pnorm(E_beta[1] + t(as.numeric(1-w[[i]]) * E_delta[, i]) %*% X[[i]][,-1] %*% E_beta[-1] +
                                       E_gamma[1] + t(as.numeric(w[[i]]) * E_delta[, i]) %*% Z[[i]][,-1] %*% E_gamma[-1]))
  P <- matrix(p, nrow = 1)
  
  for (h in 2: max_iter){
    V_beta_inv <- solve(Sigma_beta)
    m_beta <- solve(Sigma_beta)%*%mu_beta
    
    for (i in 1:n){
      w <- inputdata$w[[i]]
      Xi <- inputdata$X[[i]]
      Zi <- inputdata$Z[[i]]
      E_li <- matrix(E_delta[,i]*(1-w), ncol=1)
      E_li_prime <- matrix(E_delta[,i]*w, ncol=1)
      temp_E_delta <- matrix(E_delta[,i], ncol=1)
      E_lili <- (temp_E_delta*(1-w)) %*% t(temp_E_delta*(1-w)) + diag(V_delta[, i])
      E_lili_prime <- (temp_E_delta*(1-w)) %*% t(temp_E_delta*w)
      V_beta_inv <- V_beta_inv + t(Xi)%*%E_lili%*%Xi
      m_beta <- m_beta + E_y[i]*t(Xi)%*%(temp_E_delta*(1-w)) - t(Xi)%*%E_lili_prime%*%Zi%*%E_gamma
    }
    
    V_beta <- solve(V_beta_inv)
    m_beta <- V_beta %*% m_beta
    E_beta <- m_beta
    E_beta2 <- matrix(outer(E_beta, E_beta), nrow = length(E_beta), ncol = length(E_beta)) + V_beta
    
    V_gamma_inv <- solve(Sigma_gamma)
    m_gamma <- solve(Sigma_gamma)%*%mu_gamma
    for (i in 1: n){
      w <- inputdata$w[[i]]
      Xi <- inputdata$X[[i]]
      Zi <- inputdata$Z[[i]]
      temp_E_delta <- matrix(E_delta[,i], ncol=1)
      E_li_prime_li <- (temp_E_delta*w) %*% t(temp_E_delta*(1-w))
      E_li_prime_li_prime <- (temp_E_delta*w) %*% t(temp_E_delta*w) + diag(V_delta[, i])
      V_gamma_inv <- V_gamma_inv + t(Zi)%*%E_li_prime_li_prime%*%Zi
      m_gamma <- m_gamma + E_y[i]*t(Zi)%*%(temp_E_delta*w) - t(Zi)%*%E_li_prime_li%*%Xi%*%E_beta
    }
    V_gamma <- solve(V_gamma_inv)
    m_gamma <- V_gamma %*% m_gamma
    E_gamma <- m_gamma
    E_gamma2 <- matrix(outer(E_gamma, E_gamma), nrow = length(E_gamma), ncol = length(E_gamma)) + V_gamma
    
    V_b_inv <- solve(Sigma_b)
    m_b <- solve(Sigma_b)%*%mu_b
    
    for (i in 1:n){
      w <- inputdata$w[[i]]
      Xi <- inputdata$X[[i]]
      Zi <- inputdata$Z[[i]]
      E_U[, i] <- (1-w)*(E_b[1] + Xi[,-1]%*%b[-1]) + w*(E_d[1] + Zi[,-1]%*%E_d[-1])
      V_b_inv <- V_b_inv +t((1-w)*Xi)%*%Xi
      m_b <- m_b + t(matrix((1-w)*E_U[, i], nrow = 1) %*% Xi)
    }
    V_b <- solve(V_b_inv)
    m_b <- V_b %*% m_b
    E_b <- m_b
    E_b2 <- matrix(outer(E_b, E_b), nrow = length(E_b), ncol = length(E_b)) + V_b
    
    V_d_inv <- solve(Sigma_d)
    m_d <- solve(Sigma_d)%*%mu_d
    for (i in 1:n){
      w <- inputdata$w[[i]]
      Xi <- inputdata$X[[i]]
      Zi <- inputdata$Z[[i]]
      E_U[, i] <- (1-w)*(E_b[1] + Xi[,-1]%*%b[-1]) + w*(E_d[1] + Zi[,-1]%*%E_d[-1])
      V_d_inv <- V_d_inv +t(w*Zi)%*%Zi
      m_d <- m_d + t(matrix(w*E_U[, i], nrow = 1) %*% Zi)
    }
    V_d <- solve(V_d_inv)
    m_d <- V_d %*% m_d
    E_d <- m_d
    E_d2 <- matrix(outer(E_d, E_d), nrow = length(E_d), ncol = length(E_d)) + V_d
    
    w <- inputdata$w
    X <- inputdata$X
    Z <- inputdata$Z
    
    E_y1 <- rep(NA, n)
    E_y0 <- rep(NA, n)
    m_y <- rep(NA, n)
    m_y <- sapply(1:n, function(i) E_beta[1]+t(E_delta[,i]*(1-w[[i]]))%*%X[[i]][,-1]%*%E_beta[-1] + 
                    E_gamma[1]+t(E_delta[,i]*w[[i]])%*%Z[[i]][,-1]%*%E_gamma[-1])
    m_y[m_y > 5] <- 5
    m_y[m_y < -5] <- -5
    E_y1 <- m_y + dnorm(-m_y)/(1-pnorm(-m_y))
    E_y0 <- m_y - dnorm(-m_y)/pnorm(-m_y)
    E_y[y_vec == 0] <- E_y0[y_vec ==0]
    E_y[y_vec == 1] <- E_y1[y_vec ==1]
    
    m_U <- sapply(1:n, function(i) (1-w[[i]])*(E_b[1] + X[[i]][,-1]%*%E_b[-1])
                  +w[[i]]*(E_d[1] + Z[[i]][,-1]%*%E_d[-1]))
    E_U <- ((BB-AA)*m_U*pnorm(-m_U)+AA*m_U+(AA-BB)/sqrt(2*pi)*exp(-1/2*m_U^2))/
      (AA+(BB-AA)*pnorm(-m_U))
    E_U2 <- ((BB-AA)*(m_U^2+1)*pnorm(-m_U)+(AA-BB)*m_U/sqrt(2*pi)*exp(-1/2*m_U^2)+
               AA*(m_U^2+1))/(AA+(BB-AA)*pnorm(-m_U))
    
    ad1_ini <- sapply(1:n, function(i) E_beta[1]+E_gamma[1]+E_delta[,i]*((as.numeric(1-w[[i]]) * X[[i]][,-1]) %*% E_beta[-1] +
                                                                           (as.numeric(w[[i]]) * Z[[i]][,-1]) %*% E_gamma[-1]))
    ad1 <- matrix(rep(colSums(ad1_ini), m), byrow = TRUE, nrow = m) - ad1_ini 
    ad1 <- (matrix(rep(E_y, m), byrow = TRUE, nrow = m)-ad1)*
      sapply(1:n, function(i) E_beta[1]+ (as.numeric(1-w[[i]]) * X[[i]][,-1]) %*% E_beta[-1] + 
               E_gamma[1] + (as.numeric(w[[i]]) * Z[[i]][,-1]) %*% E_gamma[-1])
    ad2 <- sapply(1:n, function(i) 
      -diag((as.numeric(1-w[[i]]) * X[[i]]) %*% E_beta2 %*% t(X[[i]]) +
              (as.numeric(w[[i]]) * Z[[i]]) %*% E_gamma2  %*% t(Z[[i]]))/2)
    A <- AA*(1-pnorm(-m_U)) * exp(ad1) * exp(ad2)
    B <- BB*pnorm(-m_U)
    AA <- A/(A+B)
    BB <- B/(A+B)
    AA[is.nan(AA)] <- 0.5
    BB[is.nan(BB)] <- 0.5
    E_delta <- AA
    V_delta <- AA*BB
    
    p <- sapply(1:n, function(i) pnorm(E_beta[1]+ E_gamma[1] + t(as.numeric(1-w[[i]]) * E_delta[, i]) %*% X[[i]][,-1] %*% E_beta[-1] +
                                         t(as.numeric(w[[i]]) * E_delta[, i]) %*% Z[[i]][,-1] %*% E_gamma[-1]))
    P <- rbind(P, p)
    
    if (mean(abs(P[h, ] - P[h-1, ])) < epsilon_conv) { 
      print(paste("this dataset has converged at", h, "iteration"))
      converged <- 1
      break
    }
    
    if (h == max_iter) {warning("VMIL did not converge!\n")}
  }
  
  results <- list(E_beta = E_beta, E_gamma = E_gamma, E_b = E_b, E_d = E_d, 
                  V_beta = V_beta, V_gamma = V_gamma, V_b= V_b, V_d = V_d, 
                  P = as.matrix(P), converged=converged)
  return(results)
}
