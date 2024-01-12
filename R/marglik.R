#' @title Marginal likelihood of of the replication effect estimate
#'
#' @description This function computes the marginal likelihood of the
#'     replication effect estimate \code{tr} under a mixture prior for the
#'     underlying effect size with one component being a normal distribution
#'     with mean equal to the original effect estimate \code{to}, variance equal
#'     to the squared original standard error \code{so^2}, and mixture weight
#'     \code{w}, and the other component being a normal distribution with mean
#'     \code{m}, variance \code{v}, and mixture weight 1 - \code{w}. The weight
#'     \code{w} parameter can either be fixed to a value between zero and one,
#'     or a beta prior \eqn{\code{w} \sim \text{Beta}(x,y)}{\code{w} ~ Beta(x,
#'     y)} can be assigned to it.
#'
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
#' @return The marginal likelihood of the replication effect estimate.
#'
#' @author Samuel Pawel, Roberto Macri Demmartino
#'
#' @examples
#' ## uniform prior specified for weight w
#' marglik(tr = 0.09, sr = 0.05, to = 0.21, so = 0.05, x = 1, y = 1, m = 0, v = 2)
#' marglik(tr = 0.21, sr = 0.06, to = 0.21, so = 0.05, x = 1, y = 1, m = 0, v = 2)
#' marglik(tr = 0.44, sr = 0.04, to = 0.21, so = 0.05, x = 1, y = 1, m = 0, v = 2)
#'
#' ## fixed weight w = 0.5
#' marglik(tr = 0.09, sr = 0.05, to = 0.21, so = 0.05, w = 0.2, m = 0, v = 2)
#' marglik(tr = 0.21, sr = 0.06, to = 0.21, so = 0.05, w = 0.2, m = 0, v = 2)
#' marglik(tr = 0.44, sr = 0.04, to = 0.21, so = 0.05, w = 0.2, m = 0, v = 2)
#'
#' @export
marglik <- function(tr, sr, to, so, w = NULL, x = 1, y = 1, m = 0, v = 1) {
    ## input checks
    stopifnot(
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

    ## compute marginal likelihood under both models
    marg1 <- stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2))
    marg2 <- stats::dnorm(x = tr, mean = m, sd = sqrt(sr^2 + v))

    ## weight marginal likelihood with fixed weight
    if (!is.null(w)) {
        res <- w * marg1 + (1 - w) * marg2
    } else {
        ## weight marginal likelihood with random weight
        res <- x / (x + y) * (marg1 - marg2) + marg2
    }

    return(res)
}
