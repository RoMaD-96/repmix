#' @title Posterior density of effect size and weight parameter
#'
#' @description This function computes the posterior density of effect size and
#'     weight parameter based on the data from original and replication study.
#'     See the details section for details regarding data model and prior
#'     distributions.
#'
#' @details A normal likelihood around the underlying effect size \eqn{\theta}
#'     is assumed for the effect estimate from the replication study
#'     \deqn{\code{tr} \mid \theta \sim \text{N}(\theta, \code{sr}^2)}{tr |
#'     theta ~ N(theta, sr^2)}
#'
#'     A mixture prior is specified for the effect size \eqn{\theta} with one
#'     component being a normal distribution with mean equal to the original
#'     effect estimate \code{to}, variance equal to the squared original
#'     standard error \code{so^2}, and mixture weight \code{w}, and the other
#'     component being a normal distribution with mean \code{m}, variance
#'     \code{v}, and mixture weight 1 - \code{w}
#' \deqn{\theta \mid w \sim w
#'     \times \text{N}(\code{to}, \code{so}^2) + (1 - w) \times
#'     \text{N}(\code{m}, \code{v})}{ theta | w ~ w * N(to, so^2) + (1 - w) *
#'     N(m, v)}
#'
#'     A Beta prior is specified for the weight parameter
#' \deqn{w \sim \text{Beta}(x, y)}{w ~ Beta(x, y)}
#'
#' @param theta Effect size. Has to be of length one or the same length as
#'     \code{w}.
#' @param w Weight parameter. Has to be of length one or the same length as
#'     \code{theta}.
#' @param tr Effect estimate from the replication study.
#' @param sr Standard error of the replication effect estimate.
#' @param to Effect estimate from the original study.
#' @param so Standard error of the original effect estimate.
#' @param x Number of successes parameter of the beta prior for \code{w}.
#'     Defaults to \code{1}.
#' @param y Number of failures parameter of the beta prior for \code{w}.
#'     Defaults to \code{1}.
#' @param m Mean parameter of the normal prior component. Defaults to \code{0}.
#' @param v Variance parameter of the normal prior component. Defaults to
#'     \code{1}.
#'
#' @return The posterior density for a pair of effect size and weight parameter.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' ## 2D plot of posterior density
#' w <- seq(0, 1, length.out = 200)
#' theta <- seq(0, 0.3, length.out = 200)
#' parGrid <- expand.grid(w = w, theta = theta)
#' postdens <- posteriormix(theta = parGrid$theta, w = parGrid$w, tr = 0.09,
#'                     sr = 0.05, to = 0.2, so = 0.05, m = 0, v = 2)
#' postdensMat <- matrix(data = postdens, ncol = 200, byrow = TRUE)
#' filled.contour(x = theta, y = w, z = postdensMat,
#'                xlab = bquote("Effect size" ~ theta),
#'                ylab = "Weight parameter w", nlevels = 15,
#'                color.palette = function(n) hcl.colors(n = n, palette = "Blues 3", rev = TRUE))
#'
#' @export
posteriormix <- function(theta, w, tr, sr, to, so, x = 1, y = 1, m = 0, v = 1) {
    ## TODO is this a good function name?

    ## input checks
    stopifnot(
        1 <= length(theta),
        1 <= length(w),
        (length(theta) == length(w) ||
         (length(theta) == 1 & length(w) > 1) ||
         (length(theta) > 1 & length(w) == 1)
        ),

        all(is.numeric(theta)),
        all(is.finite(theta)),

        all(is.numeric(w)),
        all(is.finite(w)),
        all(0 <= w),
        all(w <= 1),

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        is.finite(v),
        0 < v
    )

    ## normalizing constant
    nc <- marglik(tr = tr, sr = sr, to = to, so = so, x = x, y = y, m = m, v = v)

    ## posterior density
    res <- stats::dnorm(x = tr, mean = theta, sd =  sr) *
        (w * stats::dnorm(x = theta, mean = to, sd = so) +
         (1 - w) * stats::dnorm(x = theta, mean = m, sd = sqrt(v))) *
        stats::dbeta(x = w, shape1 = x, shape2 = y) /
        nc

    return(res)
}
