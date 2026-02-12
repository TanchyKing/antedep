#' Cochlear implant categorical outcomes (derived)
#'
#' Categorical version of the cochlear implant speech recognition data, created
#' by mapping percentage scores to ordered 10-point bins (1..11). This dataset
#' is intended for categorical AD examples related to the cochlear application
#' discussed by Xie and Zimmerman (2013).
#'
#' @format A list with seven components:
#'   \describe{
#'     \item{y}{integer matrix of dimension N by n_time containing all subjects}
#'     \item{y_A}{integer matrix for Group A subjects}
#'     \item{y_B}{integer matrix for Group B subjects}
#'     \item{blocks}{integer vector of length N (1 = Group A, 2 = Group B)}
#'     \item{group}{character vector of group labels ("A"/"B")}
#'     \item{n_categories}{number of categories (11)}
#'     \item{category_breaks}{numeric cut points used for categorization}
#'   }
#'
#' @source Raw percentages from
#'   \url{https://www.stat.uiowa.edu/~dzimmer/Data-for-AD/speech_recognition_data.txt};
#'   categorical mapping performed in \code{data-raw/cochlear_implant_cat.R}.
"cochlear_implant_cat"
