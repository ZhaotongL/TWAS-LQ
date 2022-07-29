stage2 <- function(betahat,G,corY,n)
{
  if(abs(det(t(betahat)%*%t(G)%*%G%*%betahat/711)) < 1e-10)
  {
    return(0)
  }
  theta = solve(t(betahat)%*%t(G)%*%G%*%betahat/711)%*%t(betahat)%*%corY
  sigma2 = 
    (1 - 2*t(corY)%*%betahat%*%theta + crossprod(G%*%betahat%*%theta)/711)
  sigma2 = as.numeric(sigma2)
  vartheta = solve(crossprod(G%*%betahat)/711*n)*sigma2
  return(list(theta = theta,
              sigma2 = sigma2,
              vartheta = vartheta))
}
