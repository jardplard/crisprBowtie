#' @importFrom stats complete.cases
NULL

N_MAX_MISMATCHES <- 3

utils::globalVariables(c('.', "SpCas9", "AsCas12a"))


#' @importFrom methods is
.validateCrisprNuclease <- function(crisprNuclease){
    if (is.null(crisprNuclease)){
        crisprNuclease <- .getDefaultCrisprNuclease()
    } else {
        if (!is(crisprNuclease, "CrisprNuclease")){
            stop("Provided nuclease must be a 'CrisprNuclease' object. ")
        }
    }
    return(crisprNuclease)
}

#' @importFrom utils data
.getDefaultCrisprNuclease <- function(type=c("Cas9", "Cas12a")){
    type <- match.arg(type)
    if (type=="Cas9"){
        data("SpCas9",
             package="crisprBase",
             envir=environment())
        nuc <- SpCas9
    } else {
        data("AsCas12a",
             package="crisprBase",
             envir=environment())
        nuc <- AsCas12a
    }
    return(nuc)
}


#' @importFrom methods is
.validateBSGenome <- function(bsgenome){
    if (is.null(bsgenome)){
        stop("Provided bsgenome argument cannot be NULL.")
    } else {
        if (!is(bsgenome, "BSgenome")){
            stop("Provided bsgenome argument must be a 'BSgenome' object. ")
        }
    }
    return(bsgenome)
}




# Takes a character vector of sequences
# and writes to disk the sequences 
# in a fastq format. If temporary=TRUE,
# the fastq file is written in a 
# temporary folder. 
#' @importFrom utils write.table
.fastqfy <- function(sequences,
                     temporary=TRUE,
                     file=NULL
){
    lines <- list()
    lines[[1]] <- paste0("@", sequences)
    lines[[2]] <- sequences
    lines[[3]] <- paste0("+", sequences)
    lines[[4]] <- vapply(nchar(sequences), function(x){
                      paste0(rep('~', x), collapse='')
                  }, FUN.VALUE="a")    
    temp <- split(do.call(cbind, lines),
                  f=sequences)
    temp <- matrix(unlist(temp), ncol=1)
    if (temporary){
        file <- tempfile()
    } else {
        if (is.null(file)){
            stop("If temporary=FALSE, 'file' must be provided.")
        }
    }
    write.table(temp, 
                file=file,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE)
    return(file)
}




# Takes a character vector of sequences
# and writes to disk the sequences 
# in a fasta format. If temporary=TRUE,
# the fasta file is written in a 
# temporary folder. 
#' @importFrom utils write.table
.fastafy <- function(sequences,
                     temporary=TRUE,
                     file=NULL
){
    lines <- list()
    lines[[1]] <- paste0(">", sequences)
    lines[[2]] <- sequences
    temp <- split(do.call(cbind, lines),
                  f=sequences)
    temp <- matrix(unlist(temp), ncol=1)
    if (temporary){
        file <- tempfile()
    } else {
        if (is.null(file)){
            stop("If temporary=FALSE, 'file' must be provided.")
        }
    }
    write.table(temp, 
                file=file,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE)
    return(file)
}





.validateBowtieIndex <- function(bowtie_index) {
    check_index_type <- function(prefix, ext) {
        suffixes <- c(
                paste0(".", seq_len(4), ext),
                paste0(".rev.", seq_len(2), ext)
        )
        files <- paste0(prefix, suffixes)
        present <- files[file.exists(files)]
        missing <- files[!file.exists(files)]
        
        list(
                ext = ext,
                files = files,
                present = present,
                missing = missing,
                complete = length(present) == 6,
                partial = length(present) > 0 && length(present) < 6
        )
    }
    
    idx_ebwt  <- check_index_type(bowtie_index, ".ebwt")
    idx_ebwtl <- check_index_type(bowtie_index, ".ebwtl")
    idx_bt2   <- check_index_type(bowtie_index, ".bt2")
    idx_bt2l  <- check_index_type(bowtie_index, ".bt2l")
    
    all_present <- c(
            idx_ebwt$present, idx_ebwtl$present,
            idx_bt2$present, idx_bt2l$present
    )
    
    complete_sets <- c(
            ebwt  = idx_ebwt$complete,
            ebwtl = idx_ebwtl$complete,
            bt2   = idx_bt2$complete,
            bt2l  = idx_bt2l$complete
    )
    
    if (length(all_present) == 0) {
        stop(
                "Bowtie index not found. Please use bowtie-build or bowtie2-build ",
                "to create an index for the reference genome you are using."
        )
    }
    
    if (!any(complete_sets)) {
        missing <- c()
        if (idx_ebwt$partial)  missing <- c(missing, idx_ebwt$missing)
        if (idx_ebwtl$partial) missing <- c(missing, idx_ebwtl$missing)
        if (idx_bt2$partial)   missing <- c(missing, idx_bt2$missing)
        if (idx_bt2l$partial)  missing <- c(missing, idx_bt2l$missing)
        
        stop(
                "Only incomplete Bowtie/Bowtie2 indexes were detected for basename '",
                bowtie_index, "'. The following files are missing:\n",
                paste(missing, collapse = "\n")
        )
    }
    
    has_complete_bowtie1 <- idx_ebwt$complete || idx_ebwtl$complete
    has_complete_bowtie2 <- idx_bt2$complete || idx_bt2l$complete
    
    if (has_complete_bowtie1 && has_complete_bowtie2) {
        message(
                "Both Bowtie and Bowtie2 indexes were detected for basename '",
                bowtie_index,
                "'. Bowtie will use the Bowtie2 index."
        )
    }
    
    bowtie_index
}




#' @importFrom methods is
.checkBSGenomeOrNull <- function(bsgenome){
    if (!is.null(bsgenome)){
        if (!is(bsgenome, "BSgenome")){
            stop("Provided bsgenome argument must be a 'BSgenome' object or NULL. ")
        }
    }
    invisible(NULL)
}



