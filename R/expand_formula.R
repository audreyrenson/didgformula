#' Create formulas by expanding out string vector arguments
#'
#' @param formula An R formula containing terms referencing objects defined in the parent environment or in ... Terms must each consist of a single lowercase or uppercase letter.
#' @param ... Variables defined in `formula` (see examples).
#'
#' @return chr.
#' @export
#'
#' @examples
#' expand_formula(~a:b + d, a=c('a1','a2'), b=c('b1','b2'), d='d1')
#'
#' #we could have equivalently defined any of the variables in the formula in the calling environment:
#' a=c('a1','a2')
#' b=c('b1','b2')
#' d='d1'
#' expand_formula(~a:b + d)
expand_formula <- function(formula, ...) {

  term_maker = function(letter) {
    if(! letter %in% c(letters, LETTERS)) stop('formula terms must be a-z or A-Z')
    letter_object = eval(parse(text = letter), envir = list(...))
    if (length(letter_object) > 1) {
      return( paste0('(', paste(letter_object, collapse= " + "), ")") )
    } else {
      return ( letter_object )
    }
  }

  all_terms = attr(terms.formula(formula), 'term.labels')
  return (paste(stringr::str_replace_all(all_terms, '[a-z]', term_maker), collapse=" + "))
}
