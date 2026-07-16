###############################################################################
# geom-curve2.R
# Custom ggplot2 geom for drawing curved lines with endpoint nodes.
# Extends geom_curve to automatically add circles at start and end points.
# Author: Hou Yun (original), modified for cluster network plots.
###############################################################################

# Helper: load a function from a package if installed
#' @noRd
get_function <- function(pkg, fun) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " package has not been installed", call. = FALSE)
  }
  eval(parse(text = paste0(pkg, "::", fun)))
}

# Helper: modify an aesthetic mapping (aes) with another
#' @noRd
aes_modify <- function(aes1, aes2) {
  aes <- modifyList(as.list(aes1), as.list(aes2))
  class(aes) <- "uneval"
  aes
}

# Helper: create a data frame from a list, recycling length-1 vectors
# (Copied from ggplot2 3.3.6)
#' @noRd
new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    stop("Elements must be named", call. = FALSE)
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      stop("Elements must equal the number of rows or 1", call. = FALSE)
    }
    x[[i]] <- rep(x[[i]], n)
  }
  class(x) <- "data.frame"
  attr(x, "row.names") <- .set_row_names(n)
  x
}

# Helper: extract variable name from a mapping for a given aesthetic
#' @noRd
aes_vars <- function(mapping, vars) {
  if (vars %in% names(mapping)) {
    rlang::as_name(mapping[[vars]])
  } else {
    NULL
  }
}

# Helper: default operator (return b if a is NULL)
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# Helper: ensure list elements have names, using a prefix if missing
#' @noRd
make_list_names <- function(x, pre = "X", sep = "") {
  stopifnot(is.list(x))
  n <- length(x)
  name <- names(x)
  if (!is.null(name) && all(name != "" & !is.na(name)))
    return(x)
  if (is.null(x)) {
    names(x) <- paste0(pre, sep, seq_len(n))
  }
  if (all(name == "" | is.na(name))) {
    names(x) <- paste0(pre, sep, seq_len(n))
  } else {
    idx <- name == "" | is.na(name)
    name[idx] <- paste0(pre, sep, sum(idx))
    names(x) <- make.unique(name)
  }
  return(x)
}

# Helper: check if an object is empty (NULL or zero rows/columns)
#' @noRd
empty <- function(df) {
  if (inherits(df, "data.frame") || inherits(df, "matrix")) {
    is.null(df) || nrow(df) == 0 || ncol(df) == 0
  } else {
    is.null(df) || length(df) == 0
  }
}

# Helper: assign a name to a grid grob
#' @noRd
ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

# Helper: rename columns of a data frame
#' @noRd
rename <- function(data, ...) {
  ll <- list(...)
  if (length(ll) == 0) {
    data
  } else {
    old <- unname(unlist(ll))
    new <- names(ll)
    names(data)[names(data) %in% old] <- new
  }
  data
}

# Helper: get R version string
#' @noRd
r_version <- function() {
  strsplit(R.version.string, " ")[[1]][3]
}

# Helper: split data by group (supports list of regex/functions)
#' @noRd
split_by_group <- function(data, group) {
  n <- nrow(data)
  if (is.list(group)) {
    nm <- names(group) %||% rep("", length(group))
    out <- vector("list", length = length(group))
    
    for (ii in seq_along(group)) {
      g <- group[[ii]]
      if (is.function(g)) {
        if (nm[ii] == "") {
          stop("Invalid group parameters: should be a named list.", call. = FALSE)
        } else {
          out[[ii]] <- g(data)
        }
      } else {
        if (!is.atomic(g)) {
          stop("Invalid group parameters: should be a named list.", call. = FALSE)
        }
        if (length(g) == 1L) {
          if (nm[ii] == "") nm[ii] <- g
          out[[ii]] <- regex_select(regex = g, byrow = TRUE)(data)
        } else {
          if (nm[ii] == "") nm[ii] <- paste_with_na(g, collapse = "-")
          out[[ii]] <- regex_select(regex = g, byrow = TRUE)(data)
        }
      }
    }
    names(out) <- nm
  } else {
    out <- split(data, group)
  }
  out
}

# Helper: convert data frame to matrix (for heatmaps etc.)
#' @noRd
df_to_matrix <- function(x,
                         value,
                         row_id = NULL,
                         col_id = NULL,
                         row_names = NULL,
                         col_names = NULL,
                         missing = NA) {
  row_id <- row_id %||% names(x)[1]
  col_id <- col_id %||% names(x)[2]
  rnm <- row_names %||% unique(x[[row_id]])
  cnm <- col_names %||% unique(x[[col_id]])
  ID <- paste(rep(rnm, length(cnm)), rep(cnm, each = length(rnm)), sep = "--")
  vv <- rep(missing, length(ID))
  ii <- match(paste(x[[row_id]], x[[col_id]], sep = "--"), ID)
  vv[ii] <- x[[value]]
  matrix(vv, nrow = length(rnm), ncol = length(cnm),
         dimnames = list(rnm, cnm))
}

# Helper: get facet variables from a ggplot object
#' @noRd
get_facet_vars <- function(plot) {
  if (inherits(plot$facet, "FacetNull")) {
    return(NULL)
  }
  facet <- plot$facet$params$facets
  names(facet)
}

# ----------------------------------------------------------------------------
# Main geom_curve2 function (user-facing)
# ----------------------------------------------------------------------------
#' Curve Layer with Endpoint Nodes
#'
#' This geom draws curved lines (like geom_curve) but additionally places
#' small circles at both endpoints. Useful for network diagrams.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_curve
#' @section Aesthetics:
#' \code{geom_curve2()} understands the following aesthetics (required in bold):
#'     \itemize{
#'       \item \strong{\code{x}}
#'       \item \strong{\code{y}}
#'       \item \strong{\code{xend}}
#'       \item \strong{\code{yend}}
#'       \item \code{alpha}
#'       \item \code{colour}
#'       \item \code{linetype}
#'       \item \code{size}
#'       \item \code{curvature}
#'   }
#' @importFrom ggplot2 GeomCurve GeomPoint draw_key_path
#' @importFrom grid gTree
#' @rdname geom_curve2
#' @author Hou Yun
geom_curve2 <- function(mapping = NULL,
                        data = NULL,
                        stat = "identity",
                        position = "identity",
                        ...,
                        angle = 90,
                        ncp = 5,
                        arrow = NULL,
                        arrow.fill = NULL,
                        lineend = "butt",
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomCurve2,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      angle = angle,
      ncp = ncp,
      arrow = arrow,
      arrow.fill = arrow.fill,
      lineend = lineend,
      na.rm = na.rm,
      ...
    )
  )
}

# ----------------------------------------------------------------------------
# GeomCurve2 ggproto object (the actual geometry implementation)
# ----------------------------------------------------------------------------
#' @rdname linkET-extensions
#' @format NULL
#' @usage NULL
GeomCurve2 <- ggproto(
  "GeomCurve2", GeomCurve,
  default_aes = aes(colour = "grey35", size = 0.5, linetype = 1,
                    alpha = NA, curvature = 0),
  required_aes = c("x", "y", "xend", "yend"),
  
  # Draw panel: creates curve grobs and adds endpoint points
  draw_panel = function(self, data, panel_params, coord, drop = TRUE,
                        node.shape = 21, node.colour = "blue", node.fill = "red",
                        node.size = 2, angle = 90, ncp = 5, arrow = NULL,
                        arrow.fill = NULL, lineend = "butt", node.color = NULL,
                        na.rm = FALSE) {
    if (empty(data)) {
      return(ggplot2::zeroGrob())
    }
    
    if (!coord$is_linear()) {
      warning("`geom_curve2` is not implemented for non-linear coordinates.",
              call. = FALSE)
    }
    
    trans <- coord$transform(data, panel_params)
    
    arrow.fill <- arrow.fill %||% trans$colour
    
    # Create the curved line grobs
    grobs <- curve2Grob(x1 = trans$x,
                        y1 = trans$y,
                        x2 = trans$xend,
                        y2 = trans$yend,
                        default.units = "native",
                        curvature = trans$curvature,
                        col = scales::alpha(trans$colour, trans$alpha),
                        fill = scales::alpha(arrow.fill, trans$alpha),
                        lwd = trans$size * ggplot2::.pt,
                        lty = trans$linetype,
                        lineend = lineend,
                        angle = angle,
                        ncp = ncp,
                        square = FALSE,
                        squareShape = 1,
                        inflect = FALSE,
                        open = TRUE,
                        arrow = arrow)
    
    # Build data frames for start and end nodes
    aesthetics <- setdiff(names(data), c("x", "y", "xend", "yend", "colour",
                                         "fill", "size", "linetype"))
    if (!is.null(node.color)) {
      node.colour <- node.color
    }
    # Determine start/end node aesthetics (allow separate values)
    start.colour <- node.colour[1]
    end.colour <- if (length(node.colour) > 1) node.colour[2] else node.colour[1]
    start.fill <- node.fill[1]
    end.fill <- if (length(node.fill) > 1) node.fill[2] else node.fill[1]
    start.shape <- node.shape[1]
    end.shape <- if (length(node.shape) > 1) node.shape[2] else node.shape[1]
    start.size <- node.size[1]
    end.size <- if (length(node.size) > 1) node.size[2] else node.size[1]
    
    start.data <- new_data_frame(
      list(x = data$x,
           y = data$y,
           colour = start.colour,
           fill = start.fill,
           shape = start.shape,
           size = start.size,
           stroke = 0.5))
    end.data <- new_data_frame(
      list(x = data$xend,
           y = data$yend,
           colour = end.colour,
           fill = end.fill,
           shape = end.shape,
           size = end.size,
           stroke = 0.5))
    
    # Optionally drop duplicate nodes
    if (isTRUE(drop)) {
      start.data <- cbind(start.data, data[aesthetics])[!duplicated(start.data), , drop = FALSE]
      end.data <- cbind(end.data, data[aesthetics])[!duplicated(end.data), , drop = FALSE]
    } else {
      start.data <- cbind(start.data, data[aesthetics])
      end.data <- cbind(end.data, data[aesthetics])
    }
    
    # Combine curves and points into a single grob tree
    ggname(
      "geom_curve2",
      grid::gTree(
        children = grid::gList(
          grobs,
          GeomPoint$draw_panel(start.data, panel_params, coord),
          GeomPoint$draw_panel(end.data, panel_params, coord)
        )
      )
    )
  },
  draw_key = draw_key_path
)

# ----------------------------------------------------------------------------
# Internal function: create a single curve grob (vectorised)
# ----------------------------------------------------------------------------
#' @importFrom grid grobTree
#' @noRd
curve2Grob <- function(x1, y1, x2, y2,
                       curvature = 1,
                       col = "black",
                       fill = "grey50",
                       lwd = 0.5,
                       lty = 1,
                       lineend = "butt",
                       ...) {
  n <- max(length(x1), length(y1), length(x2), length(y2))
  
  # Recycle all arguments to length n
  x1 <- rep_len(x1, n)
  y1 <- rep_len(y1, n)
  x2 <- rep_len(x2, n)
  y2 <- rep_len(y2, n)
  curvature <- rep_len(curvature, n)
  col <- rep_len(col, n)
  fill <- rep_len(fill, n)
  lwd <- rep_len(lwd, n)
  lty <- rep_len(lty, n)
  
  # Create a grob for each curve and combine
  grobs <- lapply(seq_len(n), function(.n) {
    grid::curveGrob(x1 = x1[.n],
                    y1 = y1[.n],
                    x2 = x2[.n],
                    y2 = y2[.n],
                    curvature = curvature[.n],
                    gp = gpar(col = col[.n],
                              fill = fill[.n],
                              lwd = lwd[.n],
                              lty = lty[.n],
                              lineend = lineend),
                    ...)
  })
  do.call("grobTree", grobs)
}
