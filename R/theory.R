
#' Function to compute a theoretical z-score given source split times and admixture proportions.
#'
#' This function computes the theoretical z-score of an f4-statistic of the form f4(PO,P2,PX,P1) 
#' as a function of the admixture time a, source (P1 & P2) split time s, and admixture proportion alpha. 
#' It assumes a single lineage is sampled from each population. 
#'
#' @param t Twigstats cutoff time in units of 2Ne generations
#' @param a Admixture time in units of 2Ne generations
#' @param s Source split time in units of 2Ne generations
#' @param alpha Admixture proportions. Proportion of P2 in PX.
#' @return Returns the theoretical z-score value.
#' @export
theoretical_zscore <- function(t,a,s,alpha){
  x <- -alpha
  if(t <= s){
    x <- x * (exp(a-t) - a + t - 1)
  }else{
    foo <- ( exp(-t) * (a * (exp(s) - exp(t)) + exp(a)*(t-s) + s*exp(t) - t*exp(s)) )
    foo <- foo + exp(-(t-s)) * (exp(a-s)*(s - t + 1) - a + t - 1)
    x <- x * foo
  }

  y <- alpha
  if(t <= s){
    y <- y * (a^2 - 2*a*(t-1) - 2*exp(a-t) + t^2 - 2*t + 2)
  }else{
    foo <- ( exp(-t)*(1-exp(a-s))*( (s*(s+2)+2)*exp(t) - exp(s)*(t*(t+2)+2) ) +
            ( (a^2 + 2*a + 2)*exp(s) - exp(a)*(s^2+2*s+2) ) * exp(-s-t)*(exp(t)-exp(s)) -
            2*(exp(a)*(s+1) - exp(s)*(a+1))*exp(-s-t)*(exp(s)*(t+1)-exp(t)*(s+1)))

    foo <- foo + exp(-(t-s)) * (-exp(a-s)*(-2*(s+1)*t+s*(s+2)+t^2+2) - 2*(a+1)*t+a*(a+2)+t^2+2)
    y <- y * foo
  }

  z <- 2/3 * exp(-s+a)
  if(t > s){
    foo <- (1/18 * (18*s^2 - 27*(t^2+2*t+2)*exp(s-t) + (9*t^2+6*t+2) *exp(3*s-3*t) + 48*s + 52) +
            1/36 * (36*s^2 - 27*(2*s*(s+1)+1)*exp(s-t) + (6*t*(3*t+5) +19)*exp(3*s-3*t) + 24*s + 8) +
            1/18 * (-4*(9*s^2 + 15*s + 5)- (18*t^2 + 21*t + 7) * exp(3*s - 3*t) + 27*(2*s+1)*(t+1)*exp(s-t)))

    foo <- foo + exp(-(t-s)) * (s^2 + s*(2/3-2*t)-2/9*exp(3*s-3*t) + t^2 - 2*t/3 + 2/9)
    z <- z*foo

  }else{
    z <- 0
  }

  res <- -1/sqrt((y+z)/(x*x) - 1)

  return(res)
}

theoretical_zscore<- Vectorize(theoretical_zscore)

