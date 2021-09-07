#' @name bodyfat
#' @aliases bodyfat
#'
#' @title Percentage of body fat data set
#'
#' @description The body fat percentage of individuals assisted in a public hospital in Curitiba, Paraná, Brasil.
#'
#' @format A \code{\link{data.frame}} with 298 observations and 9 columns:
#'
#' \itemize{
#' \item \code{arms}: Arms fat percentage.
#' \item \code{legs}: Legs fat percentage.
#' \item \code{body}: Body fat percentage.
#' \item \code{android}: Android fat percentage.
#' \item \code{gynecoid}: Ginecoid fat percentage.
#' \item \code{bmi}: Body mass index - 24.71577.
#' \item \code{age}: Age - 46.00.
#' \item \code{sex}: Sex of individual. Female or male.
#' \item \code{ipaq}: Factor variable indicating the sedentary, insufficiently active or active.
#' }
#'
#' @author André Felipe Menezes \email{andrefelipemaringa@gmail.com}
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Bruna Alves \email{pg402900@uem.br}
#'
#' @usage data(bodyfat, package = "unitquantreg")
#'
#' @references
#'
#' Petterle, R. R., Bonat, W. H., Scarpin, C. T., Jonasson, T., and Borba, V. Z. C., (2020). Multivariate quasi-beta regression models for continuous bounded data. \emph{The International Journal of Biostatistics}, 1--15, (preprint).
#'
#' Mazucheli, J., Leiva, V., Alves, B., and Menezes A. F. B., (2021). A new quantile regression for modeling bounded data under a unit Birnbaum-Saunders distribution with applications in medicine and politics. \emph{Symmetry}, \bold{13}(4) 1--21.
#'
#' @source \url{http://www.leg.ufpr.br/doku.php/publications:papercompanions:multquasibeta}
#'
"bodyfat"


#' @name water
#' @aliases water
#'
#' @title Access to piped water supply data set
#'
#' @description the access of people in households with piped water supply in the cities of Brazil from
#' the Southeast and Northeast regions. Information obtained during the census of 2010.
#'
#' @format \code{\link{data.frame}} with 3457 observations and 5 columns:
#'
#' \itemize{
#' \item \code{phpws}: the proportion of households with piped water supply.
#' \item \code{mhdi}: municipal human development index.
#' \item \code{incpc}: per capita income.
#' \item \code{region}: 0 for Southeast, 1 for Northeast.
#' \item \code{pop}: total population.
#' }
#'
#' @author André Felipe Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @usage data(water, package = "unitquantreg")
#'
#' @references Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., Oliveira, R. P. and Ghitany, M. E., (2020).
#' The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the
#' modeling of quantiles conditional on covariates. \emph{Jounal of Applied Statistics}, \bold{47}(6), 954--974.
#'
"water"


