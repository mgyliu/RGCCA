# The function soft.threshold() soft-thresholds a vector such that the L1-norm constraint is satisfied.
# @param x A numeric vector.
# @param sumabs A numeric constraint on x's L1 norm.
#
# @return Returns a vector resulting from the soft thresholding of \eqn{x} given sumabs
# @keywords manip
#
soft.threshold <- function(x,sumabs=1){
  proj = proj_l1_l2(x,sumabs)
  if (proj$l2_SAT){
    x_proj = soft(x, proj$lambda)
    return(x_proj/norm2(x_proj))
  }else{
    return(proj$sol)
  }
}
