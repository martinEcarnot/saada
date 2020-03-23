saada2019_dalm_mec <-  function (Xr, Yr, Xu, Yu = NULL, weights = NULL) 
{
  Xr <- as.matrix(Xr)
  n <- nrow(Xr)
  Xu <- as.matrix(Xu)
  m <- nrow(Xu)
  rownam.Xu <- row.names(Xu)
  colnam.Yr <- colnames(Yr)
  if (is.null(colnam.Yr)) 
    colnam.Yr <- "y1"
  Yr <- as.factor(Yr)
  ni <- c(table(Yr))
  nclas <- length(ni)
  lev <- levels(Yr)
  if (!is.null(Yu)) 
    Yu <- as.character(Yu)
  else Yu <- rep(NA, m)
  if (is.null(weights)) 
    d <- rep(1/n, n)
  else d <- weights/sum(weights)
  if (nclas == 1) {
    fit <- rep(lev, m)
    dummyfit <- NULL
  }
  else {
    Xr <- cbind(rep(1, n), Xr)
    Xu <- cbind(rep(1, m), Xu)
    Xr.d <- d * Xr
    beta <- solve(crossprod(Xr.d, Xr)) %*% crossprod(Xr.d, 
                                                     dummy(Yr))
    
    # correction de la moyenne du jeu de cal
    diffm=colMeans(Xr)-colMeans(Xu)
    Xu=Xu + rep(1, nrow(Xu)) %*% t(diffm)
    
    z <- Xu %*% beta
    row.names(z) <- rownam.Xu
    colnames(z) <- lev
    dummyfit <- z
    z <- apply(dummyfit, FUN = function(x) which.max(x), 
               MARGIN = 1)
    fit <- sapply(z, FUN = function(x) lev[x])
  }
  y <- Yu
  r <- as.numeric(y != fit)
  y <- data.frame(rownum = 1:m, rownam = rownam.Xu, y, stringsAsFactors = FALSE)
  fit <- data.frame(rownum = 1:m, rownam = rownam.Xu, fit, 
                    stringsAsFactors = FALSE)
  r <- data.frame(rownum = 1:m, rownam = rownam.Xu, r)
  names(r)[ncol(r)] <- names(fit)[ncol(fit)] <- names(y)[ncol(y)] <- colnam.Yr
  list(y = y, fit = fit, r = r, dummyfit = dummyfit, ni = ni)
}