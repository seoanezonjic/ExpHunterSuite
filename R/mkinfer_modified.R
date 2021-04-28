boot.t.test.2 <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                        mu = 0, paired = FALSE, var.equal = FALSE, 
                        conf.level = 0.95, R = 9999, symmetric = FALSE, vectorial = TRUE, chunk_size = 2,...){
  alternative <- match.arg(alternative)
  if(!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                               conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.null(y)){
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else{
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }else{
    dname <- deparse(substitute(x))
    if (paired) 
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if(paired){
    x <- x - y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- if (paired) "Bootstrapped Paired t-test" else "Bootstrapped One Sample t-test"
    estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    x.cent <- x - mx
    X <- matrix(sample(x.cent, size = nx*R, replace = TRUE), nrow = R)
    MX <- rowMeans(X)
    VX <- rowSums((X-MX)^2)/(nx-1)
    STDERR <- sqrt(VX/nx)
    TSTAT <- MX/STDERR
    EFF <- MX+mx
  }else{
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2)) 
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2)) 
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3) 
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste("Bootstrapped", paste(if (!var.equal) "Welch", "Two Sample t-test"))
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if(var.equal){
      df <- nx + ny - 2
      v <- 0
      if (nx > 1) 
        v <- v + (nx - 1) * vx
      if (ny > 1) 
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
      z <- c(x, y)
      Z <- matrix(sample(z, size = (nx+ny)*R, replace = TRUE), nrow = R)
      X <- Z[,1:nx]
      Y <- Z[,(nx+1):(nx+ny)]
      MX <- rowMeans(X)
      MY <- rowMeans(Y)
      V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
      STDERR <- sqrt(V*(1/nx + 1/ny))
      EFF <- (MX+mx) - (MY+my)
    }else{
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
      z <- c(x, y)
      mz <- mean(z)
      x.cent <- x - mx + mz
      y.cent <- y - my + mz

      if (vectorial) {
         X <- matrix(sample(x.cent, size = nx*R, replace = TRUE), nrow = R)
         Y <- matrix(sample(y.cent, size = ny*R, replace = TRUE), nrow = R)
         MX <- rowMeans(X)
         MY <- rowMeans(Y)
         VX <- rowSums((X-MX)^2)/(nx-1)
         VY <- rowSums((Y-MY)^2)/(ny-1)

      } else  {
          iterations <- round(R / chunk_size)
          modulus <- R%%chunk_size
          X_dat <- data.table::rbindlist(
                    lapply(seq(iterations), function(i){
                        perm_mean_sd(x.cent, sample_size = chunk_size)}))
          MX <- X_dat$M
          VX <- X_dat$V
          X_dat <- NULL
          Y_dat <- data.table::rbindlist(
                    lapply(seq(iterations), function(i){
                        perm_mean_sd(y.cent, sample_size = chunk_size)}))
          MY <- Y_dat$M
          VY <- Y_dat$V
          Y_dat <- NULL
      }
      STDERR <- sqrt(VX/nx + VY/ny)
      EFF <- (MX+mx) - (MY+my)
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
    TSTAT <- (MX - MY)/STDERR
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    boot.pval <- mean(TSTAT < tstat)
    cint <- c(-Inf, tstat + qt(conf.level, df))
    boot.cint <- c(-Inf, quantile(EFF, conf.level))
  }else if(alternative == "greater") {
    boot.pval <- mean(TSTAT > tstat)
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
    boot.cint <- c(quantile(EFF, 1-conf.level), Inf)
  }else{
    pval <- 2 * pt(-abs(tstat), df)
    if(symmetric)
      boot.pval <- mean(abs(TSTAT) > abs(tstat))
    else
      boot.pval <- 2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat))
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
    boot.cint <- quantile(EFF, c(alpha/2, 1-alpha/2))
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(y)) "difference in means" else "mean"
  attr(cint, "conf.level") <- conf.level
  attr(boot.cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval, 
               boot.p.value = boot.pval,
               conf.int = cint, boot.conf.int = boot.cint,
               estimate = estimate, null.value = mu, 
               stderr = stderr, alternative = alternative, method = method, 
               data.name = dname)
  class(rval) <- c("boot.htest", "htest")
  rval
}
perm_mean_sd <- function(all_values, sample_size){
         X_length <- length(all_values)
         X <- matrix(sample(all_values, size = X_length*sample_size, replace = TRUE), nrow = sample_size)
         MX <- rowMeans(X)
         VX <- rowSums((X-MX)^2)/(X_length-1)
         data.frame(M = MX,V = VX)
      }
