#' @title Joint posterior density of effect size and weight parameter
#'
#' @description This function computes the posterior density of effect size and
#'     weight parameter based on the data from original and replication study
#'     using a mixture prior to incorporate the original data into the analysis
#'     of the replication data. See the details section for details regarding
#'     data model and prior distributions.
#'
#' @details A normal likelihood around the underlying effect size \eqn{\theta}
#'     is assumed for the effect estimate from the replication study
#' \deqn{\code{tr} \mid \theta \sim \text{N}(\theta, \code{sr}^2)}{tr | theta ~
#'     N(theta, sr^2)}
#'
#'     A mixture prior is assumed for the effect size \eqn{\theta} with one
#'     component being a normal distribution with mean equal to the original
#'     effect estimate \code{to}, variance equal to the squared original
#'     standard error \code{so^2}, and mixture weight \eqn{w}, and the other
#'     component being a normal distribution with mean \code{m}, variance
#'     \code{v}, and mixture weight \eqn{1 - w}
#' \deqn{\theta \mid w \sim w \times \text{N}(\code{to}, \code{so}^2) + (1 - w)
#'     \times \text{N}(\code{m}, \code{v})}{ theta | w ~ w * N(to, so^2) + (1 -
#'     w) * N(m, v)}
#'
#'     A Beta prior is assumed for the weight parameter
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
#' @return The joint posterior density for a pair of effect size and weight
#'     parameter.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' ## 2D plot of posterior density
#' w <- seq(0, 1, length.out = 200)
#' theta <- seq(0, 0.3, length.out = 200)
#' parGrid <- expand.grid(w = w, theta = theta)
#' postdens <- posteriormix(theta = parGrid$theta, w = parGrid$w, tr = 0.09,
#'                          sr = 0.05, to = 0.2, so = 0.05, m = 0, v = 4)
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
        0 < x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 < y,

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


#' @title Marginal posterior density of effect size
#'
#' @description This function computes the marginal posterior density of the
#'     effect size based on the data from original and replication study using a
#'     mixture prior to incorporate the original data into the analysis of the
#'     replication data. See the details section for details regarding data
#'     model and prior distributions.
#'
#' @details A normal likelihood around the underlying effect size \eqn{\theta}
#'     is assumed for the effect estimate from the replication study
#' \deqn{\code{tr} \mid \theta \sim \text{N}(\theta, \code{sr}^2)}{tr | theta ~
#'     N(theta, sr^2)}
#'
#'     A mixture prior is assumed for the effect size \eqn{\theta} with one
#'     component being a normal distribution with mean equal to the original
#'     effect estimate \code{to}, variance equal to the squared original
#'     standard error \code{so^2}, and mixture weight \eqn{w}, and the other
#'     component being a normal distribution with mean \code{m}, variance
#'     \code{v}, and mixture weight \eqn{1 - w}
#' \deqn{\theta \mid w \sim w \times \text{N}(\code{to}, \code{so}^2) + (1 - w)
#'     \times \text{N}(\code{m}, \code{v})}{ theta | w ~ w * N(to, so^2) + (1 -
#'     w) * N(m, v)}
#'
#'     The weight parameter can either be fixed to a specified value or a Beta
#' prior can assumed for it \deqn{w \sim \text{Beta}(x, y)}{w ~ Beta(x, y)}
#'
#' @param theta Effect size.
#' @param tr Effect estimate from the replication study.
#' @param sr Standard error of the replication effect estimate.
#' @param to Effect estimate from the original study.
#' @param so Standard error of the original effect estimate.
#' @param w Fixed weight parameter. Is only taken into account when not
#'     \code{NULL}. Defaults to \code{NULL}.
#' @param x Number of successes parameter of the beta prior for \code{w}. Only
#'     taken into account when \code{w} is \code{NULL}. Defaults to \code{1}.
#' @param y Number of failures parameter of the beta prior for \code{w}. Only
#'     taken into account when \code{w} is \code{NULL}. Defaults to \code{1}.
#' @param m Mean parameter of the normal prior component. Defaults to \code{0}.
#' @param v Variance parameter of the normal prior component. Defaults to
#'     \code{1}.
#'
#' @return The marginal posterior density for the effect size.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' theta <- seq(0, 0.3, length.out = 100)
#' post <- thetaposteriormix(theta = theta, tr = 0.09, sr = 0.05, to = 0.21, so = 0.05,
#'                           x = 1, y = 1, m = 0, v = 4)
#' plot(theta, post, type = "l", xlab = bquote("Effect size" ~ theta),
#'      ylab = "Marginal posterior density", las = 1)
#'
#' @export
thetaposteriormix <- function(theta, tr, sr, to, so, w = NULL, x = 1, y = 1, m = 0, v = 1) {
    ## TODO this is a bad name, find something better

    ## input checks
    stopifnot(
        1 <= length(theta),
        all(is.numeric(theta)),
        all(is.finite(theta)),

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

        is.null(w) || (length(w) == 1 &
                       is.numeric(w) &
                       is.finite(w) &
                       (0 <= w) &
                       (w <= 1)),

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 < x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 < y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        is.finite(v),
        0 < v
    )

    ## compute likelihood, prior components, and normalizing constant
    lik <- stats::dnorm(x = tr, mean = theta, sd = sr)
    p1 <- stats::dnorm(x = theta, mean = to, sd = so)
    p2 <- stats::dnorm(x = theta, mean = m, sd = sqrt(v))
    nc <- marglik(tr = tr, sr = sr, to = to, so = so, w = w, x = x, y = y, m = m, v = v)

    ## a random weight is the same as a fixed weight with w = x / (x + y)
    if (is.null(w)) {
        w <- x / (x + y)
    }

    res <- lik * (w * (p1 - p2) + p2) / nc
    return(res)
}


#' @title Marginal posterior density of mixture weight
#'
#' @description This function computes the marginal posterior density of the
#'     mixture weight based on the data from original and replication study
#'     using a mixture prior to incorporate the original data into the analysis
#'     of the replication data. See the details section for details regarding
#'     data model and prior distributions.
#'
#' @details A normal likelihood around the underlying effect size \eqn{\theta}
#'     is assumed for the effect estimate from the replication study
#' \deqn{\code{tr} \mid \theta \sim \text{N}(\theta, \code{sr}^2)}{tr | theta ~
#'     N(theta, sr^2)}
#'
#'     A mixture prior is assumed for the effect size \eqn{\theta} with one
#'     component being a normal distribution with mean equal to the original
#'     effect estimate \code{to}, variance equal to the squared original
#'     standard error \code{so^2}, and mixture weight \eqn{w}, and the other
#'     component being a normal distribution with mean \code{m}, variance
#'     \code{v}, and mixture weight \eqn{1 - w}
#' \deqn{\theta \mid w \sim w \times \text{N}(\code{to}, \code{so}^2) + (1 - w)
#'     \times \text{N}(\code{m}, \code{v})}{ theta | w ~ w * N(to, so^2) + (1 -
#'     w) * N(m, v)}
#'
#'     A Beta prior is assumed for the weight parameter \deqn{w \sim
#' \text{Beta}(x, y)}{w ~ Beta(x, y)}
#'
#' @param w Weight parameter.
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
#' @return The marginal posterior density for the mixture weight.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' w <- seq(0, 1, length.out = 100)
#' post <- wposteriormix(w = w, tr = 0.44, sr = 0.04, to = 0.21, so = 0.05,
#'                       x = 1, y = 1, m = 0, v = 4)
#' plot(w, post, type = "l", xlab = "Mixture weight",
#'      ylab = "Marginal posterior density", las = 1)
#'
#' @export
wposteriormix <- function(w, tr, sr, to, so, x = 1, y = 1, m = 0, v = 1) {
    ## TODO this is a bad name, find something better

    ## input checks
    stopifnot(
        1 <= length(w),
        all(is.numeric(w)),
        all(is.finite(w)),
        all((0 <= w) & (w <= 1)),

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
        0 < x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 < y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        is.finite(v),
        0 < v
    )

    ## compute likelihood, prior components, and normalizing constant
    lik1 <- stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2))
    lik2 <- stats::dnorm(x = tr, mean = m, sd = sqrt(sr^2 + v))
    prior <- stats::dbeta(x = w, shape1 = x, shape2 = y)
    nc <- marglik(tr = tr, sr = sr, to = to, so = so, x = x, y = y, m = m, v = v)

    ## compute marginal posterior density
    res <- prior * (w * lik1 + (1 - w) * lik2) / nc
    return(res)
}
