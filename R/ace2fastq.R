#' ace_to_fastq
#' 
#' Converts a sequence in .ace file format to .fastq format.
#'
#' @param filename .ace file
#' @param target_dir target directory
#' @param name2id use the file name as primary id or not. Default is TRUE.
#' @importFrom stringr str_ends str_trim str_replace
#' @return target file name
#' @author Reinhard Simon
#' @export
#'
#' @examples
#' 
#'   library(ace2fastq)
#'   filename <- system.file("sampledat/1.seq.ace", package = "ace2fastq")
#'  
#'   ace_to_fastq(filename, target_dir = tempdir())
#' 
ace_to_fastq <- function(filename,
                         target_dir = dirname(filename),
                         name2id = TRUE) {
  stopifnot(stringr::str_ends(filename, ".ace"))
  stopifnot(file.exists(filename))
  stopifnot(dir.exists(target_dir))
  stopifnot(is.logical(name2id))
  

  lines <- readLines(filename)

  # read and combine sequences
  id <- stringr::str_trim(lines[3])
  eofs <- which(lines == "")[2] # get start of sequence lines
  seqs <- paste(lines[4:eofs], collapse = "") # get all sequence lines

  # read, combine, and transform quality values
  svls <- eofs + 2 # get start of quality value lines
  evls <- which(lines == "")[3] # get end of sequence lines
  qvls <- paste(lines[svls:evls], collapse = "") # get all quality lines
  qvls <- stringr::str_split(stringr::str_trim(qvls), " ") # separate
  qvls <- as.integer(qvls[[1]]) + 33
  qvls <- paste(vapply(qvls, intToUtf8, ""), collapse = "")

  # prepare id line
  filebase <- stringr::str_replace(basename(filename), ".ace", "")
  if (name2id) {
    id <- paste0("@", filebase, " ", id)
  } else {
    id <- paste0("@", id)
  }

  # prepare resulting final fastq lines
  txt <- character(4)
  txt[1] <- id
  txt[2] <- seqs
  txt[3] <- "+"
  txt[4] <- qvls

  # write fastq
  target_name <- file.path(target_dir, paste0(filebase, ".fastq"))
  writeLines(text = txt, con = target_name)

  return(target_name)
}
