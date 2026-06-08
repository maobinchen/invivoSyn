#' Match a column name in a data frame, case-insensitive, with aliases
#'
#' @param df data frame
#' @param target canonical target column name (e.g. "Treatment")
#' @param aliases character vector of acceptable aliases
#'
#' @return the actual column name found in `df`, or `NULL` if none match
#' @keywords internal
match_col <- function(df, target, aliases = character()) {
  cols_lc <- tolower(names(df))
  hit <- which(cols_lc == tolower(target))
  if (length(hit) > 0L) {
    return(names(df)[hit[1]])
  }
  hit <- which(cols_lc %in% tolower(aliases))
  if (length(hit) > 0L) {
    return(names(df)[hit[1]])
  }
  return(NULL)
}

#' Resolve role arguments into a canonical roles list
#'
#' @param vehicle scalar character (treatment name for vehicle)
#' @param singles character vector of single-agent treatment names
#' @param combo scalar character (treatment name for the combination arm)
#' @param roles optional pre-built list (overrides individual args if supplied)
#' @param drug1,drug2 deprecated 2-drug shorthand
#'
#' @return a named list with `vehicle`, `singles`, `combo`, or `NULL` if
#'   nothing was supplied (caller should fall back to positional defaults)
#' @keywords internal
resolve_role_args <- function(vehicle = NULL, singles = NULL, combo = NULL,
                              roles = NULL, drug1 = NULL, drug2 = NULL) {
  if (!is.null(roles)) {
    stopifnot(is.list(roles), all(c("vehicle", "singles", "combo") %in% names(roles)))
    return(roles[c("vehicle", "singles", "combo")])
  }
  if (is.null(singles) && (!is.null(drug1) || !is.null(drug2))) {
    warning("`drug1`/`drug2` are deprecated; pass `singles = c(drug1, drug2)` instead.",
            call. = FALSE)
    singles <- c(drug1, drug2)
  }
  if (is.null(vehicle) && is.null(singles) && is.null(combo)) {
    return(NULL)
  }
  if (is.null(vehicle) || is.null(singles) || is.null(combo)) {
    stop("To specify roles, all of `vehicle`, `singles`, and `combo` must be supplied.",
         call. = FALSE)
  }
  return(list(vehicle = vehicle, singles = singles, combo = combo))
}

#' Build a fully-resolved roles list (with Group labels) from treatment-name roles
#'
#' @param role_args output of [resolve_role_args()], or `NULL`
#' @param trt_levels character vector of treatment names in factor-level order
#'   (i.e. the order that becomes "Group 1", "Group 2", ...)
#'
#' @return a list with `vehicle`, `singles`, `combo`, `vehicle_grp`,
#'   `single_grps`, `combo_grp`
#' @keywords internal
build_roles <- function(role_args, trt_levels) {
  if (is.null(role_args)) {
    if (length(trt_levels) < 4L) {
      stop("Need at least 4 treatments (vehicle + 2 singles + combo). Found: ",
           length(trt_levels), ".", call. = FALSE)
    }
    role_args <- list(
      vehicle = trt_levels[1L],
      singles = trt_levels[2:(length(trt_levels) - 1L)],
      combo   = trt_levels[length(trt_levels)]
    )
  }
  missing <- setdiff(c(role_args$vehicle, role_args$singles, role_args$combo), trt_levels)
  if (length(missing) > 0L) {
    stop("Treatment(s) not found in data: ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(trt_levels, collapse = ", "), call. = FALSE)
  }
  if (length(role_args$singles) < 2L) {
    stop("`singles` must contain at least 2 single-agent treatments (got ",
         length(role_args$singles), ").", call. = FALSE)
  }
  ordered_trts <- c(role_args$vehicle, role_args$singles, role_args$combo)
  grp_labels <- paste0("Group ", seq_along(ordered_trts))
  return(list(
    vehicle      = role_args$vehicle,
    singles      = role_args$singles,
    combo        = role_args$combo,
    vehicle_grp  = grp_labels[1L],
    single_grps  = grp_labels[2:(length(grp_labels) - 1L)],
    combo_grp    = grp_labels[length(grp_labels)],
    ordered_trts = ordered_trts
  ))
}

#' Sanity-check a tumor-volume tibble
#'
#' Emits warnings (not errors) for non-fatal issues such as small group sizes
#' or missing baselines. Errors only on truly broken structure.
#'
#' @param tv_long long-format tumor-volume tibble (output of [expand_tv()])
#' @param roles fully-resolved roles list (output of [build_roles()]), or `NULL`
#'   to skip role-specific checks
#'
#' @return `invisible(tv_long)`
#' @keywords internal
validate_tv <- function(tv_long, roles = NULL) {
  dup <- tv_long %>%
    dplyr::select(Mouse, Treatment) %>%
    dplyr::distinct() %>%
    dplyr::count(Mouse) %>%
    dplyr::filter(n > 1L)
  if (nrow(dup) > 0L) {
    stop("Mouse ID(s) assigned to multiple treatments: ",
         paste(dup$Mouse, collapse = ", "), call. = FALSE)
  }
  if (!is.null(roles)) {
    participating <- c(roles$vehicle_grp, roles$single_grps, roles$combo_grp)
    sizes <- tv_long %>%
      dplyr::filter(Group %in% participating) %>%
      dplyr::distinct(Group, Mouse) %>%
      dplyr::count(Group)
    thin <- sizes %>% dplyr::filter(n < 3L)
    if (nrow(thin) > 0L) {
      warning("Group(s) with <3 mice — bootstrap CIs may be unreliable: ",
              paste0(thin$Group, " (n=", thin$n, ")", collapse = ", "),
              call. = FALSE)
    }
  }
  no_baseline <- tv_long %>%
    dplyr::group_by(Mouse) %>%
    dplyr::summarise(has0 = any(Day == 0L), .groups = "drop") %>%
    dplyr::filter(!has0) %>%
    dplyr::pull(Mouse)
  if (length(no_baseline) > 0L) {
    warning("Mouse ID(s) without a Day-0 baseline (RTV/DeltaTV will be NA): ",
            paste(no_baseline, collapse = ", "), call. = FALSE)
  }
  return(invisible(tv_long))
}

#' Read tumor-volume data from a wide-format CSV
#'
#' Reads a CSV in which each row is one mouse and each numeric column is
#' tumor volume on a given day. Treatments are mapped to canonical group
#' labels ("Group 1", ...) either positionally (legacy behavior) or by
#' supplying explicit role arguments.
#'
#' @param file path to the input CSV (wide format: one column per day).
#' @param treatment_col name of the treatment column. Case-insensitive
#'   matches against `"Treatment"`, `"Group"`, `"Trt"` are accepted.
#' @param mouse_col name of the mouse-identifier column. Case-insensitive
#'   matches against `"Mouse"`, `"Subject"`, `"ID"`, `"Animal"` are accepted.
#' @param vehicle treatment name to be mapped to `"Group 1"`.
#' @param singles character vector of single-agent treatment names, in the
#'   order they should appear as `"Group 2"`, `"Group 3"`, ... Must have
#'   length >= 2.
#' @param combo treatment name for the N-drug combination arm, mapped to
#'   the last group label.
#' @param roles alternatively, a list `list(vehicle = ..., singles = ...,
#'   combo = ...)` supplying all three at once.
#' @param drug1,drug2 deprecated 2-drug shorthand; use `singles` instead.
#'
#' @return a long-format tibble with columns `Group`, `Treatment`, `Mouse`,
#'   `Day`, `TV`, `TV0`, `RTV`, `DeltaTV`, `logTV`. The fully-resolved roles
#'   list is attached as `attr(.,"roles")`.
#'
#' @export
#' @importFrom magrittr %>%
#' @import tidyr
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Legacy positional layout (vehicle first, combo last):
#' tv <- read_tv(system.file("extdata", "test.csv", package = "invivoSyn"))
#'
#' # Explicit role mapping (recommended; row order in the CSV no longer matters):
#' tv <- read_tv(
#'   system.file("extdata", "test.csv", package = "invivoSyn"),
#'   vehicle = "Vehicle",
#'   singles = c("Irinotecan", "Rabusertib"),
#'   combo   = "Irinotecan+Rabusertib"
#' )
#'
#' # Three-drug combination:
#' # tv <- read_tv("triple.csv",
#' #               vehicle = "PBS",
#' #               singles = c("A", "B", "C"),
#' #               combo   = "A+B+C")
#' }
read_tv <- function(file,
                    treatment_col = "Treatment",
                    mouse_col     = "Mouse",
                    vehicle = NULL,
                    singles = NULL,
                    combo   = NULL,
                    roles   = NULL,
                    drug1   = NULL,
                    drug2   = NULL) {
  dat <- utils::read.csv(file, check.names = FALSE)

  trt_actual   <- match_col(dat, treatment_col, c("treatment", "group", "trt"))
  mouse_actual <- match_col(dat, mouse_col, c("mouse", "subject", "id", "animal"))
  if (is.null(trt_actual)) {
    stop("Could not find treatment column (looked for: ", treatment_col,
         "). Available columns: ", paste(names(dat), collapse = ", "),
         call. = FALSE)
  }
  if (is.null(mouse_actual)) {
    stop("Could not find mouse column (looked for: ", mouse_col,
         "). Available columns: ", paste(names(dat), collapse = ", "),
         call. = FALSE)
  }
  names(dat)[names(dat) == trt_actual]   <- "Treatment"
  names(dat)[names(dat) == mouse_actual] <- "Mouse"

  role_args <- resolve_role_args(vehicle = vehicle, singles = singles,
                                 combo = combo, roles = roles,
                                 drug1 = drug1, drug2 = drug2)

  if (!is.null(role_args)) {
    trt_levels <- c(role_args$vehicle, role_args$singles, role_args$combo)
    roles_full <- build_roles(role_args, trt_levels)
    missing_in_data <- setdiff(trt_levels, unique(dat$Treatment))
    if (length(missing_in_data) > 0L) {
      stop("Treatment(s) not present in data: ",
           paste(missing_in_data, collapse = ", "),
           "\nAvailable: ", paste(unique(dat$Treatment), collapse = ", "),
           call. = FALSE)
    }
    extra <- setdiff(unique(dat$Treatment), trt_levels)
    if (length(extra) > 0L) {
      message("Ignoring ", length(extra), " non-participating treatment(s): ",
              paste(extra, collapse = ", "))
    }
    dat <- dat %>% dplyr::filter(Treatment %in% trt_levels)
  } else {
    trt_levels <- unique(dat$Treatment)
    if (length(trt_levels) == 4L) {
      message("Using positional 4-group layout (vehicle=", trt_levels[1L],
              ", singles=", paste(trt_levels[2:3], collapse = "/"),
              ", combo=", trt_levels[4L],
              "). Pass `vehicle=`, `singles=`, `combo=` to make this explicit.")
    }
    roles_full <- build_roles(NULL, trt_levels)
  }

  tv <- dat %>%
    dplyr::mutate(Treatment = factor(Treatment, levels = roles_full$ordered_trts)) %>%
    dplyr::mutate(Group = factor(paste0("Group ", as.numeric(Treatment)),
                                 levels = paste0("Group ", seq_along(roles_full$ordered_trts)))) %>%
    dplyr::select(Group, Treatment, Mouse, dplyr::everything())

  tv_long <- tv %>%
    tidyr::pivot_longer(-c(1:3), names_to = "Day", values_to = "TV") %>%
    dplyr::filter(!is.na(TV))
  tv_long <- expand_tv(tv_long)

  validate_tv(tv_long, roles_full)
  attr(tv_long, "roles") <- roles_full
  return(tv_long)
}

#' Expand a long-format tumor-volume tibble with derived columns
#'
#' Adds `TV0` (per-mouse baseline tumor volume), `logTV`, `RTV` (relative
#' tumor volume = TV / TV0), and `DeltaTV` (TV - TV0). Re-bases `Day` so
#' that each mouse's first measurement is Day 0.
#'
#' @param tv_long long-format tibble with at least `Group`, `Mouse`, `Day`, `TV`.
#'
#' @return the input tibble with `TV0`, `logTV`, `RTV`, `DeltaTV` added.
#' @export
expand_tv <- function(tv_long) {
  tv_long <- tv_long %>% dplyr::mutate(Day = as.numeric(Day))
  tv_long <- tv_long %>%
    dplyr::left_join(
      tv_long %>% dplyr::group_by(Mouse) %>% dplyr::summarise(Min_day = min(Day), .groups = "drop"),
      by = "Mouse"
    ) %>%
    dplyr::mutate(Day = Day - min(Day))
  tv_long <- tv_long %>% dplyr::select(-Min_day) %>% dplyr::filter(!is.na(TV))
  tv_long <- tv_long %>%
    dplyr::left_join(
      tv_long %>% dplyr::filter(Day == 0) %>% dplyr::select(Mouse, TV) %>% dplyr::rename(TV0 = TV),
      by = "Mouse"
    )
  tv_long <- tv_long %>%
    dplyr::mutate(logTV   = log(TV + 1),
                  RTV     = ifelse(TV0 == 0, NA_real_, TV / TV0),
                  DeltaTV = TV - TV0)
  return(tv_long)
}

#' Retrieve the resolved roles list attached to a tumor-volume tibble
#'
#' If `tv` was produced by [read_tv()] with explicit `vehicle`/`singles`/`combo`,
#' the roles list is read from `attr(tv, "roles")`. Otherwise the function
#' falls back to a positional interpretation of `levels(tv$Group)`.
#'
#' @param tv tibble returned by [read_tv()]
#'
#' @return a list with `vehicle`, `singles`, `combo`, `vehicle_grp`,
#'   `single_grps`, `combo_grp`, `ordered_trts`.
#' @export
get_roles <- function(tv) {
  r <- attr(tv, "roles")
  if (!is.null(r)) {
    return(r)
  }
  grps <- levels(tv$Group)
  if (length(grps) < 4L) {
    stop("Need at least 4 groups (vehicle + >=2 singles + combo). Found: ",
         length(grps), ". Re-run read_tv() with explicit ",
         "`vehicle`/`singles`/`combo`.", call. = FALSE)
  }
  trts <- tv %>%
    dplyr::select(Group, Treatment) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Group) %>%
    dplyr::pull(Treatment) %>%
    as.character()
  return(list(
    vehicle      = trts[1L],
    singles      = trts[2:(length(trts) - 1L)],
    combo        = trts[length(trts)],
    vehicle_grp  = grps[1L],
    single_grps  = grps[2:(length(grps) - 1L)],
    combo_grp    = grps[length(grps)],
    ordered_trts = trts
  ))
}

#' Re-tag the roles attribute on an existing tumor-volume tibble
#'
#' Useful when a single master CSV contains several combos and you want
#' to run synergy on a different subset without re-reading the file.
#' Re-orders the `Treatment` and `Group` factors to match the new roles
#' and drops non-participating rows.
#'
#' @param tv tibble returned by [read_tv()]
#' @param vehicle,singles,combo as in [read_tv()]
#'
#' @return a tibble with updated `Group`/`Treatment` factors and a new
#'   `"roles"` attribute.
#' @export
set_roles <- function(tv, vehicle, singles, combo) {
  role_args <- resolve_role_args(vehicle = vehicle, singles = singles, combo = combo)
  trt_levels <- c(role_args$vehicle, role_args$singles, role_args$combo)
  roles_full <- build_roles(role_args, trt_levels)
  trt_in_data <- unique(as.character(tv$Treatment))
  missing <- setdiff(trt_levels, trt_in_data)
  if (length(missing) > 0L) {
    stop("Treatment(s) not present in tv: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  out <- tv %>%
    dplyr::filter(as.character(Treatment) %in% trt_levels) %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment), levels = roles_full$ordered_trts)) %>%
    dplyr::mutate(Group = factor(paste0("Group ", as.numeric(Treatment)),
                                 levels = paste0("Group ", seq_along(roles_full$ordered_trts))))
  attr(out, "roles") <- roles_full
  return(out)
}
