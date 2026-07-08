#' Join metadata to a Seurat object
#'
#' Perform a left join between Seurat meta.data and an external metadata table.
#' This function supports cell-level, cluster-level, or sample-level annotation
#' by matching a user-defined key column.
#'
#' Matching is performed using a left join:
#' all cells in the Seurat object are retained, and unmatched entries
#' in the metadata are filled with NA.
#'
#' @param seu A Seurat object.
#'
#' @param metadata A data.frame or file path (.csv, .tsv, .txt, .rds)
#' containing metadata to be joined.
#'
#' @param by.x Column name in `seu@meta.data` used as the join key.
#' If NULL, a single common column name between the two tables will be used.
#'
#' @param by.y Column name in `metadata` used as the join key.
#' If NULL, defaults to `by.x`.
#'
#' @param cols Character vector specifying which columns to add from metadata.
#' If NULL, all non-key columns will be added.
#'
#' @param overwrite Logical; whether to overwrite existing columns in
#' `seu@meta.data`. Default is FALSE.
#'
#' @param strict Logical; if TRUE, the function stops when any join key
#' in Seurat has no match in metadata. If FALSE, unmatched entries are filled
#' with NA.
#'
#' @param verbose Logical; whether to print a summary of the join operation.
#'
#' @return A Seurat object with updated `meta.data`. A summary of the join
#' operation is stored in `attr(seu, "join_metadata_summary")`.
#'
#' @details
#' This function performs a base-R left join using `match()` to preserve
#' the original order of cells in the Seurat object.
#'
#' Duplicate keys in the metadata are not allowed and will trigger an error.
#'
#' Common use cases include:
#' - Cell-level annotation (e.g., barcode-based metadata)
#' - Cluster-level annotation (e.g., seurat_clusters mapping)
#' - Sample-level annotation (e.g., clinical metadata)
#'
#' @examples
#' \dontrun{
#' seu <- seurat_join_metadata(
#'     seu,
#'     metadata = "clinical.tsv",
#'     by.x = "orig.ident",
#'     by.y = "sample"
#' )
#'
#' seu <- seurat_join_metadata(
#'     seu,
#'     metadata = cluster_annot,
#'     by.x = "seurat_clusters",
#'     by.y = "cluster",
#'     strict = TRUE
#' )
#' }
#'
#' @export
seurat_join_metadata <- function(
    seu,
    metadata,
    by.x = NULL,
    by.y = NULL,
    cols = NULL,
    overwrite = FALSE,
    strict = FALSE,
    verbose = TRUE
){

    # ------------------------------------------------------------
    # 0. Check Seurat object
    # ------------------------------------------------------------
    if (!inherits(seu, "Seurat")) {
        stop("'seu' must be a Seurat object.")
    }

    # ------------------------------------------------------------
    # 1. Internal helper: read metadata
    # ------------------------------------------------------------
    .read_metadata <- function(x) {

        if (is.data.frame(x)) return(x)

        if (!is.character(x)) {
            stop("metadata must be data.frame or file path")
        }

        if (!file.exists(x)) {
            stop("Metadata file not found: ", x)
        }

        ext <- tolower(tools::file_ext(x))

        df <- switch(
            ext,
            csv = utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE),
            tsv = utils::read.delim(x, stringsAsFactors = FALSE, check.names = FALSE),
            txt = utils::read.delim(x, stringsAsFactors = FALSE, check.names = FALSE),
            rds = readRDS(x),
            stop("Unsupported file format: ", ext)
        )

        if (!is.data.frame(df)) {
            stop("Invalid metadata format")
        }

        return(df)
    }

    # ------------------------------------------------------------
    # 2. Internal helper: validate metadata
    # ------------------------------------------------------------
    .validate_metadata <- function(df) {

        if (anyDuplicated(colnames(df))) {
            stop("Duplicated column names in metadata")
        }

        if (any(colnames(df) == "")) {
            stop("Empty column names in metadata")
        }

        TRUE
    }

    # ------------------------------------------------------------
    # 3. Load metadata
    # ------------------------------------------------------------
    metadata <- .read_metadata(metadata)
    .validate_metadata(metadata)

    meta.x <- seu[[]]

    # ------------------------------------------------------------
    # 4. Determine join keys
    # ------------------------------------------------------------
    if (is.null(by.x)) {

        common <- intersect(colnames(meta.x), colnames(metadata))

        if (length(common) == 0) {
            stop("No common column names found. Please specify by.x and by.y.")
        }

        if (length(common) > 1) {
            stop(
                "Multiple common columns found: ",
                paste(common, collapse = ", "),
                ". Please specify by.x and by.y."
            )
        }

        by.x <- common
    }

    if (is.null(by.y)) {
        by.y <- by.x
    }

    if (!by.x %in% colnames(meta.x)) {
        stop("by.x not found in Seurat meta.data: ", by.x)
    }

    if (!by.y %in% colnames(metadata)) {
        stop("by.y not found in metadata: ", by.y)
    }

    # ------------------------------------------------------------
    # 5. Convert keys to character (avoid factor bugs)
    # ------------------------------------------------------------
    meta.x[[by.x]] <- as.character(meta.x[[by.x]])
    metadata[[by.y]] <- as.character(metadata[[by.y]])

    # ------------------------------------------------------------
    # 6. Duplicate check in metadata key
    # ------------------------------------------------------------
    if (anyDuplicated(metadata[[by.y]])) {

        dup <- unique(metadata[[by.y]][duplicated(metadata[[by.y]])])

        stop(
            "Duplicated keys in metadata (by.y). Example: ",
            paste(head(dup, 10), collapse = ", ")
        )
    }

    # ------------------------------------------------------------
    # 7. Determine columns to add
    # ------------------------------------------------------------
    if (is.null(cols)) {
        cols <- setdiff(colnames(metadata), by.y)
    } else {
        missing.cols <- setdiff(cols, colnames(metadata))
        if (length(missing.cols) > 0) {
            stop("Missing columns in metadata: ",
                 paste(missing.cols, collapse = ", "))
        }
    }

    # overwrite handling
    skipped <- character(0)

    if (!overwrite) {
        skipped <- intersect(cols, colnames(meta.x))
        cols <- setdiff(cols, skipped)
    }

    if (length(cols) == 0) {
        if (verbose) message("No columns added.")
        return(seu)
    }

    # ------------------------------------------------------------
    # 8. Join
    # ------------------------------------------------------------
    idx <- match(meta.x[[by.x]], metadata[[by.y]])

    matched <- sum(!is.na(idx))
    unmatched <- sum(is.na(idx))

    if (strict && unmatched > 0) {
        miss <- unique(meta.x[[by.x]][is.na(idx)])
        stop(
            "Unmatched keys found (strict=TRUE). Example: ",
            paste(head(miss, 10), collapse = ", ")
        )
    }

    # ------------------------------------------------------------
    # 9. Add metadata
    # ------------------------------------------------------------
    for (cc in cols) {
        meta.x[[cc]] <- metadata[[cc]][idx]
    }

    seu[[]] <- meta.x

    # ------------------------------------------------------------
    # 10. Summary
    # ------------------------------------------------------------
#    summary <- list(
#        join = c(by.x = by.x, by.y = by.y),
#        n_cells = nrow(meta.x),
#        n_metadata = nrow(metadata),
#        matched = matched,
#        unmatched = unmatched,
#        added_columns = cols,
#        skipped_columns = skipped
#    )
#
#    attr(seu, "join_metadata_summary") <- summary

    if (verbose) {

        cat("\n")
        cat("====================================\n")
        cat("seurat_join_metadata v1\n")
        cat("====================================\n")

        cat("Join:\n")
        cat("  Seurat:", by.x, "\n")
        cat("  Metadata:", by.y, "\n\n")

        cat("Cells:", nrow(meta.x), "\n")
        cat("Metadata rows:", nrow(metadata), "\n")
        cat("Matched:", matched, "\n")
        cat("Unmatched:", unmatched, "\n\n")

        cat("Added columns:\n")
        if (length(cols) > 0) {
            cat("  ", paste(cols, collapse = ", "), "\n")
        }

        if (length(skipped) > 0) {
            cat("\nSkipped columns:\n")
            cat("  ", paste(skipped, collapse = ", "), "\n")
        }

        cat("====================================\n\n")
    }

    return(seu)
}

