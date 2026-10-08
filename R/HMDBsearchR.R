#' Submit 1H NMR peaks to HMDB for compound search 
#'
#' This function allows to search HMDB from within R by simply submitting the peaks you want to search for. 
#' The function will open a webpage with the query results or provide a link to the HMDB page with the results.
#'
#' @param peakVector A vector with ppm values of peaks
#' @param ppmTol The ppm tolerance for the HMDB search (default = 0.02).
#' @param returnURL Return the URL instead of opening a webpage.
#' 
#' @return Opens a webpage or returns a URL with the HMDB results
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' \dontrun{
#' HMDBsearchR(peakVector = c(3.2, 3.38), ppmTol = 0.2, returnURL = TRUE)
#' }
#' 
#' @importFrom rvest session html_form html_form_set session_submit
#' @importFrom xml2 read_html
#' @importFrom utils browseURL      
#'         
#' @export
#' 
HMDBsearchR <- function(peakVector, ppmTol = 0.02, returnURL = FALSE){
    
    peakSubmit <- paste(as.character(peakVector), collapse = " \n ")
    
    session <- rvest::session("http://www.hmdb.ca/spectra/nmr/one_d/search/new")
    form <- rvest::html_form(xml2::read_html("http://www.hmdb.ca/spectra/nmr/one_d/search/new"))[[2]]
    filled_form <- rvest::html_form_set(form, peaks = peakSubmit, cs_tolerance = ppmTol)
    result <- rvest::session_submit(session, filled_form, submit = 'peaks')
    
    HMDBResultsURL = result$url
    if(returnURL){
        return(HMDBResultsURL)
    }else{
        utils::browseURL(HMDBResultsURL, browser = getOption("browser"))
    }
}