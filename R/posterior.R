#' @title Posterior based on two component normal mixture prior
#'
#' @description This function computes the posterior means, variances, and the
#'     weight parameter based on data with normal likelihood and a two component
#'     mixture prior for the underlying mean.
#'
#' @param y Data mean.
#' @param se Standard error of the data mean.
#' @param m1 Mean parameter of the first normal prior component.
#' @param v1 Variance parameter of first the normal prior component.
#' @param m2 Mean parameter of the second normal prior component.
#' @param v2 Variance parameter of second the normal prior component.
#' @param w Mixture weight.
#'
#' @return A list of the posterior means, variances, and weight.
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## normal2mixposterior(y = 0.09, se = 0.05, m1 = 0.21, v1 = 0.05^2,
#' ##                     m2 = 0, v2 = 1, w = 0.5)
#'
#' @keywords internal
normal2mixposterior <- function(y, se, m1, v1, m2, v2, w) {
    ## updated weight
    wpost <- 1 /
        (1 +
         ## priors odds
         (1 - w) / w *
         ## Bayes factor
         stats::dnorm(x = y, mean = m2, sd = sqrt(se^2 + v2))/   
         stats::dnorm(x = y, mean = m1, sd = sqrt(se^2 + v1))
         
        )
    ## updated component means and variances
    v1post <- 1/(1/v1 + 1/se^2)
    m1post <- (m1/v1 + y/se^2)*v1post
    v2post <- 1/(1/v2 + 1/se^2)
    m2post <- (m2/v2 + y/se^2)*v2post
    res <- list("w" = wpost, "m1" = m1post, "v1" = v1post, "m2" = m2post,
                "v2" = v2post)
    return(res)
}

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
thetaposteriormix <- function(theta, tr, sr, to, so, w = NULL, x = 1, y = 1,
                              m = 0, v = 1) {
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

    ## a random weight is the same as a fixed weight with w = x / (x + y)
    if (is.null(w)) {
        w <- x / (x + y)
    }

    ## compute posterior density
    postpars <- normal2mixposterior(y = tr, se = sr, m1 = to, v1 = so^2,
                                    m2 = m, v2 = v, w = w)
    res <- postpars$w * stats::dnorm(x = theta, mean = postpars$m1,
                                     sd = sqrt(postpars$v1)) +
        (1 - postpars$w) * stats::dnorm(x = theta, mean = postpars$m2,
                                     sd = sqrt(postpars$v2))
    return(res)
}

#' @title Highest posterior density interval for effect size
#'
#' @description This function computes the highest posterior density interval
#'     for the effect size based on the data from original and replication study
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
#'     The weight parameter can either be fixed to a specified value or a Beta
#' prior can assumed for it \deqn{w \sim \text{Beta}(x, y)}{w ~ Beta(x, y)}
#'
#' @param level Credibility level. Defaults to \code{0.95}.
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
#' @return The highest posterior density interval at the specified credibility
#'     level and posterior median for the effect size.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' ## tipping point analysis (computing HPDs for range of weights from 0 to 1)
#' to <- 0.21
#' tr <- 0.09
#' so <- 0.05
#' sr <- 0.045
#' w <- seq(0, 1, length.out = 20)
#' HPDs <- t(sapply(X = w, FUN = function(w) {
#'                thetaHPD(level = 0.95, tr = tr, sr = sr, to = to,
#'                         so = so, w = w, m = 0, v = 4)
#' }))
#' plot(NA, NA, type = "n", ylim = c(0, 0.3), xlim = c(-0.25, 1.25),
#'      xlab = "Weight of original study",
#'      ylab = "Posterior median with 95% HPDI", las = 1)
#' arrows(x0 = w, y0 = HPDs[,1], y1 = HPDs[,3], code = 3, angle = 90, length = 0.05,
#'        col = 4)
#' points(x = w, y = HPDs[,2], pch = 20, col = 4)
#' arrows(x0 = c(-0.2, 1.2), y0 = c(tr - 1.96*sr, to - 1.96*so),
#'        y1 = c(tr + 1.96*sr, to + 1.96*so), code = 3, angle = 90, length = 0.05,
#'        col = c(1, 2))
#' points(x = c(-0.2, 1.2), y = c(tr, to), pch = 20, col = c(1, 2))
#' axis(side = 1, at = c(-0.2, 1.2), labels = c("Replication", "Original"))
#'
#' @export
thetaHPD <- function(level = 0.95, tr, sr, to, so, w = NULL, x = 1, y = 1,
                     m = 0, v = 1) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

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
    ## a random weight is the same as a fixed weight with w = x / (x + y)
    if (is.null(w)) {
        w <- x / (x + y)
    }

    ## compute parameters of posterior
    postpars <- normal2mixposterior(y = tr, se = sr, m1 = to, v1 = so^2,
                                    m2 = m, v2 = v, w = w)

    ## HPD is known when weight is zero or one
    if (w == 1) {
        CI <- stats::qnorm(p = c((1 - level)*0.5, 0.5, (1 + level)*0.5),
                           mean = postpars$m1, sd = sqrt(postpars$v1))
    } else if (w == 0) {
        CI <- stats::qnorm(p = c((1 - level)*0.5, 0.5, (1 + level)*0.5),
                           mean = postpars$m2, sd = sqrt(postpars$v2))
    } else {
        ## determine HPD numerically

        ## CDF of the posterior
        cdf <- function(theta) {
            postpars$w * stats::pnorm(q = theta, mean = postpars$m1,
                                      sd = sqrt(postpars$v1)) +
                (1 - postpars$w) * stats::pnorm(q = theta, mean = postpars$m2,
                                                sd = sqrt(postpars$v2))
        }

        ## quantile function of the posterior
        quantileFun. <- function(p) {
            if (p == 0) {
                res <- -Inf
            } else if (p == 1) {
                res <- Inf
            } else {
                q1 <- stats::qnorm(p = p, mean = postpars$m1, sd = sqrt(postpars$v1))
                q2 <- stats::qnorm(p = p, mean = postpars$m2, sd = sqrt(postpars$v2))
                ## quantile must be somewhere between the quantiles from the
                ## components
                rootFun <- function(theta) cdf(theta) - p
                res <- stats::uniroot(f = rootFun,
                                      interval = c(min(c(q1, q2)), max(c(q1, q2)))
                                      )$root
            }
            return(res)
        }
        quantileFun <- Vectorize(FUN = quantileFun.)

        ## search for the shortest interval
        optFun <- function(lowerq) {
            width <- quantileFun(p = lowerq + level) - quantileFun(p = lowerq)
            return(width)
        }
        res <- try(stats::optim(par = (1 - level)*0.5, fn = optFun,
                                method = "L-BFGS-B",
                                ## TODO make argument/interface to control search range
                                lower = (1 - level)*0.001,
                                upper = (1 - level)*0.999)$par)
        if (inherits(res, "try-error")) {
            CI <- c(NaN, NaN, NaN)
        } else {
            CI <- c(quantileFun(p = res),
                    quantileFun(p = 0.5),
                    quantileFun(p = res + level))
        }
    }

    names(CI) <- c("lower", "median", "upper")
    return(CI)
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


#' @title Highest posterior density interval for the mixture weight
#'
#' @description This function computes the highest posterior density interval
#'     for the mixture weight based on the data from original and replication study
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
#' @param level Credibility level. Defaults to \code{0.95}.
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
#' @return The highest posterior density interval at the specified credibility
#'     level and posterior median for the mixture weight.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' wHPD(level = 0.95, tr = 0.09, sr = 0.05, to = 0.21, so = 0.05, x = 1, y = 1,
#'       m = 0, v = 4)
#' wHPD(level = 0.95, tr = 0.21, sr = 0.06, to = 0.21, so = 0.05, x = 1, y = 1,
#'       m = 0, v = 4)
#' wHPD(level = 0.95, tr = 0.44, sr = 0.04, to = 0.21, so = 0.05, x = 1, y = 1,
#'       m = 0, v = 4)
#'
#' @export
wHPD <- function(level = 0.95, tr, sr, to, so, x = 1, y = 1, m = 0, v = 1) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

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

    ## CDF of the posterior
    cdf <- function(w) {
        res <- try(stats::integrate(f = wposteriormix, lower = 0, upper = w,
                                    tr = tr, sr = sr, to = to, so = so, x = x,
                                    y = y, m = m, v = v)$value)
        if (inherits(res, "try-error")) return(NaN)
        else return(res)
    }

    ## quantile function of the posterior
    quantileFun. <- function(p) {
        if (p == 0) {
            res <- 0
        } else if (p == 1) {
            res <- 1
        } else {
            rootFun <- function(theta) cdf(theta) - p
            res <- stats::uniroot(f = rootFun, interval = c(0, 1))$root
        }
        return(res)
    }
    quantileFun <- Vectorize(FUN = quantileFun.)

    ## search for the shortest interval
    optFun <- function(lowerq) {
        width <- quantileFun(p = lowerq + level) - quantileFun(p = lowerq)
        return(width)
    }
    res <- try(stats::optim(par = (1 - level)*0.5, fn = optFun, method = "L-BFGS-B",
                            lower = 0, upper = 1 - level)$par)
    if (inherits(res, "try-error")) {
        CI <- c(NaN, NaN, NaN)
    } else {
        CI <- c(quantileFun(p = res),
                quantileFun(p = 0.5),
                quantileFun(p = res + level))
    }


    names(CI) <- c("lower", "median", "upper")
    return(CI)
}
