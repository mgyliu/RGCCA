# Other checking functions
#-----------------------------
check_blockx <- function(x, y, blocks) {
    message = paste0(x, " should be lower than ", length(blocks),
      " (that is the number of blocks)."
    )
    exit_code = 133
    x <- check_integer(x, y, max = length(blocks), exit_code = exit_code,
                       max_message = message)
    return(x)
}

check_boolean <- function(x, y = x, type = "scalar") {

    if (is.null(y))
        y <- x

    if (any(is.na(y)))
        stop_rgcca(paste(x, "should not be NA."))

    if (!is(y, "logical"))
        stop_rgcca(paste(x, "should be TRUE or FALSE."))

    if (type == "scalar" && length(y) != 1)
        stop_rgcca(paste(x, "should be of length 1."))
}

check_colors <- function(colors){

    if (!is.null(colors)){
        colors <- as.vector(colors)
        for (i in colors) {
            if (!is.na(i) && !(i  %in%  colors()) && is.character2(i)
                && regexpr("^#{1}[a-zA-Z0-9]{6,8}$", i) < 1)
                stop_rgcca("colors must be in colors() or a rgb character.")
        }
    }
}

check_compx <- function(x, y, ncomp, blockx) {
    res <- check_integer(x, y, min = 1)
    if (y > ncomp[blockx]) {
        stop_rgcca(
            paste0(
                x,
                " equals to ",
                y,
                " and should be comprise between 1 and ",
                ncomp[blockx],
                " (the number of block component)."
            ),
            exit_code = 128
        )
    }
 return(res)
}

# Check the format of the connection matrix
#
# @inheritParams rgccad
# @inheritParams set_connection
check_connection <- function(C, blocks, software = FALSE) {

    msg <- ifelse(software, "The connection file", "The design matrix C should")

        if (!all(rownames(C) %in% names(blocks)) ||
        !all(colnames(C) %in% names(blocks)))
        stop_rgcca(paste(msg,
                         "have the rownames and the colnames that match with",
                         "the names of the blocks."),
                   exit_code = 108)

    if (!isSymmetric.matrix(unname(C)))
        stop_rgcca(paste(msg, "be symmetric."), exit_code = 103)

    # d <- unique(diag(C))
    # if (length(d) != 1 || d != 0)
    #     stop_rgcca("The diagonal of the connection matrix file should be 0.",
    #         exit_code = 105)

    if (any(is.na(C)))
        stop_rgcca(paste(msg, "not contain NA values."), exit_code = 105)

    x <- C >= 0 & C <= 1
    if (sum(!x) != 0)
        stop_rgcca(paste(msg, "contain numbers between 0 and 1."), exit_code = 106)

    if (all(C == 0))
        stop_rgcca(paste(msg, "not contain only 0."), exit_code = 107)

    invisible(check_size_blocks(blocks, "connection matrix", C))

    if(is.null(rownames(C)) || is.null(colnames(C)))
        rownames(C) <- names(blocks) -> colnames(C)

    return(C)
    # TODO: warning if superblock = TRUE
}
check_file <- function(f) {
    # Check the existence of a path f: A character giving the path of a file

    if (!file.exists(f))
        stop_rgcca(paste(f, "does not exist."), exit_code = 101)

}
# y <- runif(6)
# check_integer("y", y, "vector", T, min = 0)
# y <- matrix(runif(6, 1, 10), 2, 3)
# check_integer("y", y, "matrix", T)
# check_integer("y", y, "matrix")
# check_integer(NA)
# check_integer(.1)
# check_integer(c(1:2))
# check_integer("x", c(0, 0), type = "vector")
check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1,
                          max = Inf, max_message = NULL, exit_code = NULL,
                          min_message = NULL) {

    if (is.null(y))
        y <- x

    if (type %in% c("matrix", "data.frame"))
        y_temp <- y

    y <- tryCatch(
        as.double(as.matrix(y)),
        warning = function(w)
            stop_rgcca(paste(x, "should be numeric."))
    )

    if (any(is.na(y)))
        stop_rgcca(paste(x, "should not be NA."))

    if (!is(y, "numeric"))
        stop_rgcca(paste(x, "should be numeric."))

    if (type == "scalar" && length(y) != 1)
        stop_rgcca(paste(x, "should be of length 1."))

    if (!float) {
      if (any((y %% 1) != 0)) {
        stop_rgcca(paste(x, "should be an integer."))
      }
      y = as.integer(y)
    }

    if (any(y < min))
      if (!is.null(min_message)) {
        stop_rgcca(min_message, exit_code = exit_code)
      } else {
        stop_rgcca(paste0(x, " should be higher than or equal to ", min, "."))
      }

    if (any(y > max))
      if (!is.null(max_message)) {
        stop_rgcca(max_message, exit_code = exit_code)
      } else {
        stop_rgcca(paste0(x, " should be lower than or equal to ", max, "."), exit_code = 150)
      }

    if (type %in% c("matrix", "data.frame"))
        y <- matrix(
            y,
            dim(y_temp)[1],
            dim(y_temp)[2],
            dimnames = dimnames(y_temp)
        )

    if (type == "data.frame")
        y = as.data.frame(y)

    return(y)
}

check_method <- function(method) {
    analysis <- c("rgcca", "sgcca", "pca", "spca", "pls", "spls",
      "cca", "ifa", "ra", "gcca", "maxvar", "maxvar-b",
      "maxvar-a", "mcoa","cpca-1", "cpca-2", "cpca-4",
      "hpca", "maxbet-b", "maxbet", "maxdiff-b", "maxdiff",
      "maxvar-a", "sabscor", "ssqcor", "ssqcor", "ssqcov-1",
      "ssqcov-2", "ssqcov", "sumcor", "sumcov-1", "sumcov-2",
      "sumcov", "sabscov", "sabscov-1", "sabscov-2")

    if (!tolower(method) %in% analysis)
        stop_rgcca(
            paste0("Please select one type among the following
            type: ", paste(analysis, collapse = ", ")),
            exit_code = 112
        )
}

check_nblocks <- function(blocks, method) {
    if (tolower(method) == "pca") {
        if (length(blocks) == 1) return(blocks)
        nb <- 1
        exit_code <- 110
    } else{
        if (length(blocks) == 2) return(blocks)
        nb <- 2
        exit_code <- 111
    }

    stop_rgcca(
        paste0(
            length(blocks),
            " blocks were provided but the number of blocks for ", method,
            " must be ", nb, "."
        ),
        exit_code = exit_code
    )
}
# If less than 2 columns, do not run
# x, a list of matrix
# i_block, the position of the tested matrix in the list
# return an error or NULL
check_ncol <- function(x, i_block) {
    if (NROW(x[[i_block]]) < 2) {
        stop_rgcca(
            "This output is available only for more than one variable."
        )
    }
}

check_ncomp <- function(ncomp, blocks, min = 1) {
    ncomp <- elongate_arg(ncomp, blocks)
    check_size_blocks(blocks, "ncomp", ncomp)
    ncomp <- sapply(
        seq(length(ncomp)),
        function(x){
            msg = paste0("ncomp[", x, "] should be comprise between ", min ," and ",
                         NCOL(blocks[[x]]),
                         " (that is the number of variables of block ", x, ")."
            )
            y <- check_integer("ncomp", ncomp[x], min = min, max_message = msg,
                               max = NCOL(blocks[[x]]), exit_code = 126)
            return(y)
        }
    )
    return(ncomp)
}

# Check if a dataframe contains no qualitative variables
# @inheritParams load_blocks
# @param df A dataframe or a matrix
# @param fo A character giving the name of the tested file
# @param warn_separator A bolean to print warning for bad separator use
check_quantitative <- function(df, fo, header = FALSE, warn_separator = FALSE) {
    qualitative <- is.character2(df, warn_separator = TRUE)

    if (qualitative) {
        msg <- paste(
            fo,
            "contains qualitative data. Please, transform them in a disjunctive table."
        )

        if (!header)
            msg <- paste0(msg, "Possible mistake: header parameter is disabled, check if the file doesn't have one.")

        stop_rgcca(paste(msg, "\n"), exit_code = 100)
    }

}

check_response <- function(response = NULL, df = NULL) {

    if (!is.null(response)) {
        qualitative <- is.character2(response, type = "all")

        # if (length(qualitative) > 1)
        #     stop_rgcca(
        #     "Please, select a response file with either qualitative data only or quantitative data only.",
        #     108
        #     )

        if (!qualitative)
            response <- to_numeric(response)
        if (NCOL(response) > 1) {
            disjunctive <- unique(apply(response, 1, sum))

            if (length(disjunctive) &&
                unique(disjunctive %in% c(0, 1)) && disjunctive) {

                response2 <- factor(apply(response, 1, which.max))

                if (!is.null(colnames(response)))
                    levels(response2) <- colnames(response)

                return(
                    as.matrix(
                        data.frame(
                            as.character(response2),
                            row.names = rownames(response)
                        )))

            } else {
                warning("There is multiple columns in the response block. By default, only the first column will be considered.")
                return(as.matrix(response[, 1]))
            }
        }

        return(response)

    } else
        return(rep(1, NROW(df[[1]])))

}

# Test on the sign of the correlation
check_sign_comp <- function(rgcca_res, w){

    w1 <- rgcca_res$a
    y1  <- lapply(seq_along(w1),
                 function(i) pm(rgcca_res$call$blocks[[i]],
                                        rgcca_res$a[[i]]))
    y <- lapply(seq_along(w1),
                 function(i) pm(rgcca_res$call$blocks[[i]], w[[i]]))


    for (k in seq(length(w))) {
        if(NCOL(w[[k]])>1)
        {
            for (j in seq(NCOL(w[[k]]))) {
                res = ifelse(NROW(rgcca_res$a[[k]]) < NROW(rgcca_res$Y[[k]]),
                             cor(y1[[k]][, j], y[[k]][, j]),
                             cor(w1[[k]][, j], w[[k]][, j]))
                if (!is.na(res) && res  < 0)
                    w[[k]][, j] <- -1 * w[[k]][, j]
            }
        }
        else
        {
          res = ifelse(NROW(rgcca_res$a[[k]]) < NROW(rgcca_res$Y[[k]]),
                       cor(y1[[k]], y[[k]]),
                       cor(w1[[k]], w[[k]])
                       )
          if (!is.na(res) && res  < 0)
          w[[k]] <- -1 * w[[k]]
        }
    }

    return(w)
}

check_size_blocks <- function(blocks, x, y = x) {

    if (identical(x, y))
        x <- ""
    if (any(class(y) %in% c("matrix", "data.frame"))) {
        dim_y <- NCOL(y)
        dim_type <- "number of columns"
    }else{
        dim_y <- length(y)
        dim_type <- "size"
    }

    if (dim_y != length(blocks))
        stop_rgcca(
            paste0(
                x,
                " should have the same ",
                dim_type ,
                " (actually ",
                dim_y,
                ") as the number of blocks (",
                length(blocks),
                ")."
            ),
            exit_code = 130
        )
    else
        return(TRUE)
}


# Print warning if file size over
check_size_file <- function(filename) {
    size <- file.size(filename)
    if (size > 5e+06)
        # warning(paste0('The size of ', filename, ' is over 5 Mo (',
        #  round(size / 1E6, 1), ' Mo). File loading could take some times...'),
        message("File loading in progress ...")
}

check_penalty <- function(penalty, blocks, method = "rgcca", superblock = FALSE) {
  if (superblock) {
    blocks[[length(blocks) + 1]] <- Reduce(cbind,blocks)
    names(blocks)[length(blocks)] = "superblock"
  }
  penalty <- elongate_arg(penalty, blocks)
  name = ifelse(method == "rgcca", "tau", "sparsity")
  check_size_blocks(blocks, name, penalty)
  penalty1 <- penalty

  is_matrix = is(penalty, "matrix")

  # Check value of each penalty
  if (method == "rgcca") penalty <- sapply(penalty, check_tau, USE.NAMES = FALSE)
  if (method == "sgcca") {
    if (is_matrix) divider = NROW(penalty1)
    else divider = 1
    penalty <- sapply(
      seq(length(penalty)),
      function(x) check_spars(penalty[x], blocks[[1 + (x - 1) / divider]]))
  }

  if (is(penalty1, "matrix"))
    penalty <- matrix(penalty, NROW(penalty1), NCOL(penalty1))

  return(penalty)
}

check_spars <- function(sparsity, block) {
  min_sparsity <- 1 / sqrt(NCOL(block))
  min_message <- paste0("Sparsity parameter equals to ", sparsity,
                        ". For SGCCA, it must be greater than ",
                        "1/sqrt(number_column) (i.e., ", min_sparsity, ").")
  sparsity <- check_integer("sparsity", sparsity, float = TRUE,
                            min = min_sparsity, max = 1, min_message = min_message)
  invisible(sparsity)
}

# #' @export
check_superblock <- function(is_supervised = NULL, is_superblock = NULL, verbose = TRUE) {
    if (!is.null(is_supervised)) {
        if (verbose)
            warn_connection("supersized method with a response")
        if (is_superblock) {
            if (!is.null(is_superblock) && verbose)
                warning("In a supervised mode, the superblock corresponds to the response.")
        }
        return(FALSE)
    }else
        return(isTRUE(is_superblock))
}
check_tau <- function(tau) {
    msg <- "tau should be comprise between 0 and 1 or should correspond to the character 'optimal' for automatic setting"
    if (is.na(tau) || tau != "optimal") {
        tryCatch({
            tau <- check_integer("tau", tau, float = TRUE, min = 0, max = 1)
        }, error = function(e) {
            if (class(e)[1] == 150)
                stop_rgcca(paste0(msg, " (currently equals to ", tau, ")."), exit_code = 129)
              else
                stop_rgcca(e$message)
        })
    }
    invisible(tau)
}

check_scheme <- function(scheme) {
  if (
    (mode(scheme) != "function") &&
    (scheme != "horst") &&
    (scheme != "factorial") &&
    (scheme != "centroid")
  ) stop_rgcca(paste0("Choose one of the three following schemes: horst, ",
                      "centroid, factorial or design the g function"))
}
