# Code for Atkinson et al. CJFAS Maunscript
# Manuscript Title: Hydrology-mediated shifts in home range size and spatial overlap for two mesopredators in the coastal Florida Everglades

# ==== 1) Packages ====
pkgs <- c(
  "dplyr","purrr","readr","lubridate","sf","amt","ggplot2","scales",
  "ggh4x","viridisLite","mgcv","tidyr","lme4","stringr","tibble"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ==== 2) Paths & data ====

#Load data
bass_rds <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis/AllBass_POR_0612025.rds"
snook_rds <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis/AllSnook_POR_0612025.rds"
SharkRiver <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis/SharkRiver_UTMZone17N.shp"
data_dir <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis"
water_csv <- file.path(data_dir, "MO-215_Jan1990-Sep2025.csv")

df_bass  <- readRDS(bass_rds);  if (!"Species" %in% names(df_bass))  df_bass$Species  <- "Bass"
df_snook <- readRDS(snook_rds); if (!"Species" %in% names(df_snook)) df_snook$Species <- "Snook"
df_water <- read_csv(water_csv, show_col_types = FALSE)


# Output paths
outdir_root  <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Outputs/For Manuscript_Round1"
outdir_plots <- file.path(outdir_root, "plots")
outdir_spatial <- file.path(outdir_root, "spatial")
outdir_csv <- file.path(outdir_root, "csv")
out_gpkg <- file.path(outdir_spatial, "home_ranges_monthly.gpkg")
out_csv_all <- file.path(outdir_csv, "home_ranges_monthly_summary.csv")
out_csv_bass <- file.path(outdir_csv, "home_ranges_monthly_bass.csv")
out_csv_snook <- file.path(outdir_csv, "home_ranges_monthly_snook.csv")
out_csv_indiv <- file.path(outdir_csv, "home_ranges_by_individual.csv")
out_csv_indiv_b <- file.path(outdir_csv, "home_ranges_by_individual_bass.csv")
out_csv_indiv_s <- file.path(outdir_csv, "home_ranges_by_individual_snook.csv")


# ==== 3) Spatial setup ====
crs_target_epsg <- 26917  # NAD83 / UTM 17N (meters)
min_coas <- 30
kde_levels <- c(0.50, 0.95)  
kde_bw <- "href"            
river_buffer_m <- 100        

# ==== 4) HOURLY COAs ====
prep_and_coa <- function(df) {
  needed <- c("Transmitter","Datetime_UTC","Station","Latitude","Longitude")
  miss <- setdiff(needed, names(df)); if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  if ("Distance" %in% names(df)) df$Distance <- suppressWarnings(as.numeric(df$Distance))
  df <- df %>%
    mutate(
      Datetime_UTC = suppressWarnings(lubridate::as_datetime(Datetime_UTC, tz = "UTC")),
      Latitude = suppressWarnings(as.numeric(Latitude)),
      Longitude = suppressWarnings(as.numeric(Longitude))
    ) %>%
    filter(!is.na(Datetime_UTC), !is.na(Latitude), !is.na(Longitude))
  
  pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326) %>% st_transform(crs_target_epsg)
  xy <- st_coordinates(pts); pts$x_m <- xy[,1]; pts$y_m <- xy[,2]
  has_distance <- "Distance" %in% names(df)
  
  pts %>%
    mutate(COA_time = lubridate::floor_date(Datetime_UTC, unit = "hour")) %>%
    st_drop_geometry() %>%
    group_by(Species, Transmitter, COA_time) %>%
    summarise(
      n_detections = n(),
      x_m = mean(x_m, na.rm = TRUE),
      y_m = mean(y_m, na.rm = TRUE),
      Station_count = n_distinct(Station),
      Distance_mean_km = if (has_distance) mean(Distance, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("x_m","y_m"), crs = crs_target_epsg) %>%
    mutate(Month = lubridate::floor_date(COA_time, unit = "month"))
}
coas_bass  <- prep_and_coa(df_bass)
coas_snook <- prep_and_coa(df_snook)
coas_all   <- bind_rows(coas_bass, coas_snook)

# ==== 5) MCPs & KDEs ====
grp_list <- coas_all %>% group_by(Species, Transmitter, Month) %>% group_split()

hr_results <- purrr::map(grp_list, function(g){
  if (nrow(g) < min_coas) return(NULL)
  coords <- st_coordinates(g)
  n_unique <- nrow(unique(round(coords, 0)))
  dat <- g %>%
    mutate(x = coords[,1], y = coords[,2], t = COA_time) %>%
    st_drop_geometry()
  
  trk <- amt::make_track(dat, x, y, t, crs = paste0("EPSG:", crs_target_epsg))
  
  
  mcp_poly <- if (n_unique >= 3) {
    tmp <- try(amt::hr_isopleths(amt::hr_mcp(trk, levels = 0.95)), silent = TRUE)
    if (inherits(tmp,"try-error")) NULL else tmp
  } else NULL
  
 
  kde_try <- NULL
  if (n_unique >= 2) {
    kde_try <- try(amt::hr_isopleths(amt::hr_kde(trk, h = kde_bw, levels = kde_levels)), silent = TRUE)
    if (inherits(kde_try,"try-error"))
      kde_try <- try(amt::hr_isopleths(amt::hr_kde(trk, h = 150, levels = kde_levels)), silent = TRUE)
    if (inherits(kde_try,"try-error")) kde_try <- NULL
  }
  
  list(
    key = list(Species=unique(g$Species), Transmitter=unique(g$Transmitter),
               Month=unique(g$Month), n_COAs=nrow(g)),
    mcp = mcp_poly, kde = kde_try
  )
}) %>% purrr::compact()

# Bind
bind_with_keys <- function(lst, which = c("mcp","kde")) {
  which <- match.arg(which)
  parts <- purrr::map(lst, function(x){
    obj <- x[[which]]; if (is.null(obj)) return(NULL)
    obj$Species <- x$key$Species; obj$Transmitter <- x$key$Transmitter
    obj$Month <- x$key$Month; obj$n_COAs <- x$key$n_COAs; obj
  }) %>% purrr::compact()
  if (!length(parts)) return(NULL); do.call(rbind, parts)
}
mcp95_raw <- bind_with_keys(hr_results, "mcp")
kde_raw   <- bind_with_keys(hr_results, "kde")

normalize_levels <- function(sfobj){
  if (is.null(sfobj) || nrow(sfobj)==0) return(sfobj)
  lvl <- suppressWarnings(as.numeric(as.character(sfobj$level)))
  if (max(lvl, na.rm = TRUE) > 1) lvl <- lvl / 100
  sfobj$level <- lvl; sfobj
}
mcp95_raw <- normalize_levels(mcp95_raw)
kde_raw   <- normalize_levels(kde_raw)
kde50_raw <- if (!is.null(kde_raw) && nrow(kde_raw)) dplyr::filter(kde_raw, dplyr::near(level, 0.5, tol=1e-8)) else NULL
kde95_raw <- if (!is.null(kde_raw) && nrow(kde_raw)) dplyr::filter(kde_raw, dplyr::near(level, 0.95, tol=1e-8)) else NULL

# ==== 6) Clip to river ====
river <- st_read(SharkRiver, quiet = TRUE) %>% st_transform(crs_target_epsg)
river <- if ("st_make_valid" %in% getNamespaceExports("sf")) suppressWarnings(st_make_valid(river)) else suppressWarnings(st_buffer(river, 0))
types <- unique(as.character(st_geometry_type(river, by_geometry = TRUE)))
if (!all(types %in% c("POLYGON","MULTIPOLYGON"))) {
  message("SharkRiver geometry is ", paste(types, collapse=", "), " — buffering by ", river_buffer_m, " m to create a clipping mask.")
  river <- river %>% st_buffer(river_buffer_m) %>% st_union() %>% st_cast("MULTIPOLYGON")
}


clip_end <- function(x){
  if (is.null(x) || nrow(x)==0) return(NULL)
  out <- try(suppressWarnings(st_intersection(x, river)), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out) || nrow(out)==0) return(NULL)
  out
}
mcp95_all <- clip_end(mcp95_raw)
kde50_all <- clip_end(kde50_raw)
kde95_all <- clip_end(kde95_raw)

# ==== 7) HR area tables ====
area_tbl <- function(sf_poly, method_label){
  if (is.null(sf_poly) || nrow(sf_poly)==0) return(tibble())
  sf_poly %>%
    mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
    st_drop_geometry() %>%
    transmute(Species, Transmitter, Month, n_COAs,
              level = suppressWarnings(as.numeric(level)),
              method = method_label, area_km2)
}
summary_tbl <- bind_rows(
  area_tbl(mcp95_all, "MCP95"),
  area_tbl(kde50_all, "KDE50"),
  area_tbl(kde95_all, "KDE95")
) %>% arrange(Species, Transmitter, Month, method)


write_csv(summary_tbl, out_csv_all)
write_csv(filter(summary_tbl, Species=="Bass"), out_csv_bass)
write_csv(filter(summary_tbl, Species=="Snook"), out_csv_snook)

# By individual
by_indiv <- summary_tbl %>%
  group_by(Species, Transmitter, method) %>%
  summarise(
    n_months = n_distinct(Month),
    n_COAs_total = sum(n_COAs, na.rm = TRUE),
    mean_area_km2 = mean(area_km2, na.rm = TRUE),
    median_area_km2 = median(area_km2, na.rm = TRUE),
    sd_area_km2 = sd(area_km2, na.rm = TRUE),
    min_area_km2 = min(area_km2, na.rm = TRUE),
    max_area_km2 = max(area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Species, Transmitter, method)
write_csv(by_indiv, out_csv_indiv)
write_csv(filter(by_indiv, Species=="Bass"),  out_csv_indiv_b)
write_csv(filter(by_indiv, Species=="Snook"), out_csv_indiv_s)

# ==== 8) Output HR spatial outputs ====
if (file.exists(out_gpkg)) file.remove(out_gpkg)
write_layer <- function(obj, layer_name){
  if (!is.null(obj) && nrow(obj)>0) st_write(obj, out_gpkg, layer=layer_name, delete_layer=TRUE, quiet=TRUE)
}
write_layer(mcp95_all, "MCP95")
write_layer(kde50_all, "KDE50")
write_layer(kde95_all, "KDE95")
write_layer(mcp95_raw, "MCP95_raw")
write_layer(kde50_raw, "KDE50_raw")
write_layer(kde95_raw, "KDE95_raw")

cat("\nWrote CSVs to: ", normalizePath(outdir_csv, winslash="/"),
    "\nWrote GPKG : ", normalizePath(out_gpkg, winslash="/"), "\n", sep="")

# ==== 9) Interspecific overlap ====
dissolve_species_month <- function(sfobj){
  if (is.null(sfobj) || nrow(sfobj)==0) return(NULL)
  sfobj |> mutate(Month = as.Date(Month)) |>
    group_by(Species, Month) |>
    summarise(geometry = st_union(geometry), .groups="drop")
}
mcp95_sm <- dissolve_species_month(mcp95_all)
kde50_sm <- dissolve_species_month(kde50_all)
kde95_sm <- dissolve_species_month(kde95_all)

compute_overlap_series <- function(sm_obj, method_name){
  if (is.null(sm_obj) || nrow(sm_obj) == 0) return(tibble())
  B <- filter(sm_obj, Species == "Bass") |> select(Month, geometry)
  S <- filter(sm_obj, Species == "Snook") |> select(Month, geometry)
  if (!nrow(B) || !nrow(S)) return(tibble())
  
  B_tbl <- tibble(Month = B$Month, geom_B = st_geometry(B))
  S_tbl <- tibble(Month = S$Month, geom_S = st_geometry(S))
  df <- inner_join(B_tbl, S_tbl, by = "Month")
  
  out <- df %>%
    rowwise() %>%
    mutate(
      A_B = as.numeric(st_area(st_union(geom_B))) / 1e6,
      A_S = as.numeric(st_area(st_union(geom_S))) / 1e6,
      A_I = {
        gi <- try(suppressWarnings(st_intersection(geom_B, geom_S)), silent = TRUE)
        if (inherits(gi,"try-error") || length(gi) == 0) 0 else {
          gi <- gi[!st_is_empty(gi)]
          if (length(gi) == 0) 0 else as.numeric(st_area(st_union(gi))) / 1e6
        }
      },
      A_U = pmax(A_B + A_S - A_I, 0)
    ) %>% ungroup() %>%
    mutate(
      jaccard = ifelse(A_U > 0, A_I / A_U, 0),
      pct_overlap_bass  = ifelse(A_B > 0, A_I / A_B, 0),
      pct_overlap_snook = ifelse(A_S > 0, A_I / A_S, 0),
      overlap_coef = ifelse(pmin(A_B, A_S) > 0, A_I / pmin(A_B, A_S), 0),
      method = method_name
    ) %>% select(Month, method, A_B, A_S, A_I, A_U, jaccard,
                 pct_overlap_bass, pct_overlap_snook, overlap_coef)
  out
}
overlap_mcp95 <- compute_overlap_series(mcp95_sm, "MCP95")
overlap_kde50 <- compute_overlap_series(kde50_sm, "KDE50")
overlap_kde95 <- compute_overlap_series(kde95_sm, "KDE95")


write_csv(overlap_kde50, file.path(outdir_csv, "overlap_kde50_species_month.csv"))
write_csv(overlap_kde95, file.path(outdir_csv, "overlap_kde95_species_month.csv"))
write_csv(overlap_mcp95, file.path(outdir_csv, "overlap_mcp95_species_month.csv"))

# plot
plot_overlap_ts <- function(df, method, yvar, ylab, file_stub){
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df <- arrange(df, Month)
  p <- ggplot(df, aes(Month, .data[[yvar]])) +
    geom_line(linewidth = 0.7, colour = "black") +
    geom_point(size = 1.8, colour = "black") +
    scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = paste(method, ylab), x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_overlap_ts
  
  fn <- file.path(outdir_plots, paste0(file_stub, "_", tolower(method), ".tiff"))
  
  ggsave(fn, p, device = "tiff", width = 9, height = 7, units = "in",
         dpi = 600, compression = "lzw")
}
plot_overlap_ts(overlap_kde50, "KDE50", "jaccard", "Jaccard overlap", "overlap_jaccard")
plot_overlap_ts(overlap_kde95, "KDE95", "jaccard", "Jaccard overlap", "overlap_jaccard")
plot_overlap_ts(overlap_mcp95, "MCP95", "jaccard", "Jaccard overlap", "overlap_jaccard")



# ==== 10) Intraspecific overlap ====
.area_km2 <- function(g){ if (length(g)==0 || all(st_is_empty(g))) return(0); as.numeric(st_area(g))/1e6 }
per_indiv_month <- function(sf_in){
  if (is.null(sf_in) || nrow(sf_in) == 0) return(NULL)
  sf_in |> mutate(Month = as.Date(Month)) |>
    group_by(Species, Transmitter, Month) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
}
compute_intra_pairs <- function(sf_in, method_label, species_name){
  pim <- per_indiv_month(sf_in) |> filter(Species == species_name)
  if (is.null(pim) || nrow(pim) == 0) return(tibble())
  groups <- pim |> group_by(Month) |> group_split()
  purrr::map_dfr(groups, function(g){
    if (nrow(g) < 2) return(tibble())
    cmb <- utils::combn(seq_len(nrow(g)), 2, simplify = FALSE)
    purrr::map_dfr(cmb, function(ix){
      a <- g[ix[1], ]; b <- g[ix[2], ]
      A_A <- .area_km2(a$geometry); A_B <- .area_km2(b$geometry)
      Iab <- try(suppressWarnings(st_intersection(a$geometry, b$geometry)), silent = TRUE)
      A_I <- if (inherits(Iab,"try-error") || length(Iab)==0) 0 else {
        Iab <- Iab[!st_is_empty(Iab)]; if (length(Iab)==0) 0 else .area_km2(st_union(Iab))
      }
      A_U <- max(A_A + A_B - A_I, 0)
      tibble(
        Species = species_name, Month = a$Month,
        method = method_label, A_id = a$Transmitter, B_id = b$Transmitter,
        A_km2 = A_A, B_km2 = A_B, I_km2 = A_I, U_km2 = A_U,
        jaccard = ifelse(A_U > 0, A_I / A_U, 0),
        pct_A = ifelse(A_A > 0, A_I / A_A, 0),
        pct_B = ifelse(A_B > 0, A_I / A_B, 0),
        overlap_coef = ifelse(pmin(A_A, A_B) > 0, A_I / pmin(A_A, A_B), 0)
      )
    })
  })
}
intra_monthly_overlap_polys <- function(kde_sf, species, min_overlap_m2 = 1){
  if (is.null(kde_sf) || !inherits(kde_sf, "sf") || nrow(kde_sf)==0) return(NULL)
  x <- kde_sf %>% suppressWarnings(st_make_valid()) %>%
    suppressWarnings(st_collection_extract("POLYGON"))
  x <- x[!st_is_empty(x), , drop = FALSE]; if (!nrow(x)) return(NULL)
  need <- c("Species","Transmitter","Month")
  miss <- setdiff(need, names(x)); if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  x <- filter(x, .data$Species == species)
  if (!nrow(x)) return(NULL)
  if (!inherits(x$Month, "Date")) x$Month <- as.Date(x$Month)
  x_u <- x %>% group_by(Month, Transmitter) %>% summarise(geometry = st_union(st_geometry(.)), .groups="drop")
  x_u <- x_u[!st_is_empty(x_u), , drop = FALSE]; if (!nrow(x_u)) return(NULL)
  x_u$area_m2 <- as.numeric(st_area(x_u))
  crs_use <- st_crs(x_u)
  out_list <- list(); by_m <- split(x_u, x_u$Month)
  for (mn in names(by_m)) {
    dfm <- by_m[[mn]]; if (nrow(dfm) < 2) next
    idx <- utils::combn(nrow(dfm), 2)
    for (k in seq_len(ncol(idx))) {
      ia <- idx[1,k]; ib <- idx[2,k]
      A <- dfm[ia, , drop=FALSE]; B <- dfm[ib, , drop=FALSE]
      inter <- try(suppressWarnings(st_intersection(A, B)), silent=TRUE)
      if (inherits(inter,"try-error") || is.null(inter) || nrow(inter)==0) next
      inter <- suppressWarnings(st_make_valid(inter))
      inter <- suppressWarnings(st_collection_extract(inter, "POLYGON"))
      inter <- inter[!st_is_empty(inter), , drop = FALSE]; if (nrow(inter)==0) next
      g_union <- suppressWarnings(st_union(st_geometry(inter)))
      if (length(g_union)==0 || any(st_is_empty(g_union))) next
      overlap_m2 <- as.numeric(st_area(g_union)); if (!is.finite(overlap_m2) || overlap_m2 < min_overlap_m2) next
      ab_union <- try(suppressWarnings(st_union(st_geometry(A), st_geometry(B))), silent=TRUE)
      union_m2 <- if (inherits(ab_union,"try-error")) NA_real_ else as.numeric(st_area(ab_union))
      geom_sfc <- if (inherits(g_union,"sfc")) g_union else st_sfc(g_union, crs = crs_use)
      row_df <- tibble(
        Species = species, Month = A$Month[[1]],
        Transmitter_A = A$Transmitter[[1]], Transmitter_B = B$Transmitter[[1]],
        area_A_m2 = A$area_m2[[1]], area_B_m2 = B$area_m2[[1]],
        overlap_m2 = overlap_m2,
        prop_A = overlap_m2 / A$area_m2[[1]],
        prop_B = overlap_m2 / B$area_m2[[1]],
        jaccard = if (is.finite(union_m2) && union_m2 > 0) overlap_m2 / union_m2 else NA_real_
      )
      out_list[[length(out_list)+1]] <- st_sf(row_df, geometry = geom_sfc, crs = crs_use)
    }
  }
  if (!length(out_list)) return(NULL); do.call(rbind, out_list)
}
plot_intra_maps <- function(overlap_sf, species_sm, species_name, method_label, outline_col){
  if (is.null(overlap_sf) || nrow(overlap_sf)==0) return(invisible(NULL))
  tmp <- overlap_sf |> mutate(A_km2 = as.numeric(st_area(geometry))/1e6) |>
    arrange(desc(A_km2)) |> slice_head(n=6)
  if (!nrow(tmp)) return(invisible(NULL))
  outline <- species_sm |> filter(Species == species_name, Month %in% tmp$Month) |>
    mutate(Month_lab = format(Month, "%Y-%m"))
  fill <- tmp |> mutate(Month_lab = format(Month, "%Y-%m"))
  p <- ggplot() +
    { if (exists("river") && !is.null(river) && nrow(river)) geom_sf(data=river, fill=NA, color="grey70", linewidth=0.3) } +
    geom_sf(data = fill, fill = "#6a51a3", color = NA, alpha = 0.35) +
    geom_sf(data = outline, fill = NA, color = outline_col, linewidth = 0.9) +
    facet_wrap(~ Month_lab, ncol = 3) + coord_sf(expand = FALSE) +
    labs(title = paste0(species_name, " intra-species overlap (", method_label, "): top months"),
         subtitle = "Shaded = any pairwise overlap; outline = species union",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"))
  fn <- file.path(outdir_plots, paste0("intra_", tolower(species_name), "_maps_", tolower(method_label), "_top6.tiff"))
  ggsave(fn, p, device = "tiff", width = 10, height = 7, units = "in", dpi = 600, compression = "lzw")
  message("Saved: ", fn)
}
plot_intra_heatmap <- function(overlap_sf, method_label, species_name, cell_m = 50){
  if (is.null(overlap_sf) || nrow(overlap_sf) == 0) return(invisible(NULL))
  grid_sf <- overlap_heatmap_grid(overlap_sf, river, cell_m = cell_m)
  if (is.null(grid_sf) || nrow(grid_sf) == 0) return(invisible(NULL))
  grid_sf <- mutate(grid_sf, fill_val = ifelse(pct_months > 0, pct_months, NA_real_))
  p <- ggplot() +
    { if (exists("river") && !is.null(river) && nrow(river)) geom_sf(data=river, fill=NA, color="grey70", linewidth=0.3) } +
    geom_sf(data = grid_sf, aes(fill = fill_val), color = NA) +
    scale_fill_viridis_c(name = "% months overlapped", labels = percent_format(accuracy = 1),
                         limits = c(0, 1), oob = scales::squish, na.value = NA) +
    coord_sf(expand = FALSE) +
    labs(title = paste0(species_name, " intra-species overlap frequency (", method_label, ")"),
         subtitle = paste0("Grid = ", cell_m, " m; colored where ≥ 2 individuals overlapped in a month"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank(), legend.position = "right")
  fn <- file.path(outdir_plots, paste0("intra_", tolower(species_name), "_heatmap_", tolower(method_label), "_", cell_m, "m.tiff"))
  ggsave(fn, p, device = "tiff", width = 9, height = 7, units = "in", dpi = 600, compression = "lzw")
  message("Saved: ", fn)
  st_write(grid_sf, out_gpkg, layer = paste0("intra_", tolower(species_name), "_heat_", tolower(method_label), "_", cell_m, "m"),
           delete_layer = TRUE, quiet = TRUE)
}

# Intraspecific: bass KDE50/95 
bass_pairs_50 <- compute_intra_pairs(kde50_all, "KDE50", "Bass")
write_csv(bass_pairs_50, file.path(outdir_csv, "intra_bass_pairwise_kde50.csv"))
bass_month_50 <- bass_pairs_50 |> group_by(Species, Month, method) |>
  summarise(n_pairs = n(), median_jaccard = median(jaccard, na.rm = TRUE),
            mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop") |>
  arrange(Month)
write_csv(bass_month_50, file.path(outdir_csv, "intra_bass_monthly_kde50.csv"))
bass_overlap_polys_50 <- intra_monthly_overlap_polys(kde50_all, "Bass")
plot_intra_maps(bass_overlap_polys_50, kde50_sm, "Bass", "KDE50", outline_col = "#377eb8")
plot_intra_heatmap(bass_overlap_polys_50, "KDE50", "Bass", cell_m = 50)

bass_pairs_95 <- compute_intra_pairs(kde95_all, "KDE95", "Bass")
write_csv(bass_pairs_95, file.path(outdir_csv, "intra_bass_pairwise_kde95.csv"))
bass_month_95 <- bass_pairs_95 |> group_by(Species, Month, method) |>
  summarise(n_pairs = n(), median_jaccard = median(jaccard, na.rm = TRUE),
            mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop") |>
  arrange(Month)
write_csv(bass_month_95, file.path(outdir_csv, "intra_bass_monthly_kde95.csv"))
bass_overlap_polys_95 <- intra_monthly_overlap_polys(kde95_all, "Bass")
plot_intra_maps(bass_overlap_polys_95, kde95_sm, "Bass", "KDE95", outline_col = "#377eb8")
plot_intra_heatmap(bass_overlap_polys_95, "KDE95", "Bass", cell_m = 50)

# Intraspecific: snook KDE50/95 
snook_pairs_50 <- compute_intra_pairs(kde50_all, "KDE50", "Snook")
write_csv(snook_pairs_50, file.path(outdir_csv, "intra_snook_pairwise_kde50.csv"))
snook_month_50 <- snook_pairs_50 |> group_by(Species, Month, method) |>
  summarise(n_pairs = n(), median_jaccard = median(jaccard, na.rm = TRUE),
            mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop") |>
  arrange(Month)
write_csv(snook_month_50, file.path(outdir_csv, "intra_snook_monthly_kde50.csv"))
snook_overlap_polys_50 <- intra_monthly_overlap_polys(kde50_all, "Snook")
plot_intra_maps(snook_overlap_polys_50, kde50_sm, "Snook", "KDE50", outline_col = "#e6550d")
plot_intra_heatmap(snook_overlap_polys_50, "KDE50", "Snook", cell_m = 50)

snook_pairs_95 <- compute_intra_pairs(kde95_all, "KDE95", "Snook")
write_csv(snook_pairs_95, file.path(outdir_csv, "intra_snook_pairwise_kde95.csv"))
snook_month_95 <- snook_pairs_95 |> group_by(Species, Month, method) |>
  summarise(n_pairs = n(), median_jaccard = median(jaccard, na.rm = TRUE),
            mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop") |>
  arrange(Month)
write_csv(snook_month_95, file.path(outdir_csv, "intra_snook_monthly_kde95.csv"))
snook_overlap_polys_95 <- intra_monthly_overlap_polys(kde95_all, "Snook")
plot_intra_maps(snook_overlap_polys_95, kde95_sm, "Snook", "KDE95", outline_col = "#e6550d")
plot_intra_heatmap(snook_overlap_polys_95, "KDE95", "Snook", cell_m = 50)

# ==== 11) Water level ====

wl_monthly <- df_water %>%
  dplyr::mutate(
    WaterDate = lubridate::mdy(.data[["Date"]]),  
    WaterVal  = suppressWarnings(as.numeric(.data[["WaterLevel(cm)"]]))
  ) %>%
  dplyr::filter(!is.na(WaterDate), is.finite(WaterVal)) %>%
  dplyr::mutate(Month = lubridate::floor_date(WaterDate, "month")) %>%
  dplyr::group_by(Month) %>%
  dplyr::summarise(WaterLevel = mean(WaterVal, na.rm = TRUE), .groups = "drop")


# ==== 12) HR vs. water level (GAMs)  ====

std_home_range <- function(df){
  nm <- names(df); nml <- tolower(trimws(nm)); names(df) <- trimws(nm)
  if (!"month" %in% nml) stop("df_home_range is missing 'Month'.")
  mcol <- nm[which(nml == "month")[1]]; df$Month <- as.Date(df[[mcol]])
  if ("area_km2" %in% nml) {
    acol <- nm[which(nml == "area_km2")[1]]; df$area_km2 <- suppressWarnings(as.numeric(df[[acol]]))
  } else if ("area_m2" %in% nml) {
    acol <- nm[which(nml == "area_m2")[1]]; df$area_km2 <- suppressWarnings(as.numeric(df[[acol]])) / 1e6
  } else stop("Need 'area_km2' or 'area_m2'.")
  if (!"method" %in% nml) stop("df_home_range missing 'method'.")
  mth <- nm[which(nml == "method")[1]]; df[[mth]] <- toupper(as.character(df[[mth]])); names(df)[names(df)==mth] <- "method"
  if ("transmitter" %in% nml) names(df)[which(nml=="transmitter")[1]] <- "Transmitter"
  if (!"species" %in% nml) stop("df_home_range missing 'Species'.")
  df
}
df_home_range <- std_home_range(summary_tbl)

gam_dat <- df_home_range %>%
  filter(method %in% c("KDE50","KDE95")) %>%
  mutate(Month = as.Date(Month)) %>%
  left_join(wl_monthly, by = "Month") %>%
  filter(!is.na(WaterLevel), !is.na(area_km2)) %>%
  mutate(
    Species = factor(Species),
    method = factor(method, levels = c("KDE50","KDE95")),
    area_log1p = log1p(area_km2)
  )

fit_gam <- function(df, species_name, method_name, k_smooth = 3){
  d <- df %>% filter(Species == species_name, method == method_name)
  if (!nrow(d)) stop("No rows for ", species_name, " / ", method_name, " after joins/filters.")
  d$WaterLevel <- as.numeric(d$WaterLevel)
  nuniq <- length(unique(d$WaterLevel[is.finite(d$WaterLevel)]))
  k_use <- max(3, min(k_smooth, nuniq - 1))
  use_re <- "Transmitter" %in% names(d) && length(unique(d$Transmitter)) > 1
  if (use_re) d$Transmitter <- factor(d$Transmitter)
  form <- if (use_re) area_log1p ~ s(WaterLevel, k=k_use) + s(Transmitter, bs="re") else area_log1p ~ s(WaterLevel, k=k_use)
  m <- mgcv::gam(form, data = d, method = "REML")
  newx <- tibble(WaterLevel = seq(min(d$WaterLevel), max(d$WaterLevel), length.out = 200))
  if (use_re) newx$Transmitter <- levels(d$Transmitter)[1]
  sm_labels <- if (length(m$smooth)) vapply(m$smooth, function(s) s$label, character(1)) else character(0)
  excl_idx <- if (use_re) which(grepl("Transmitter", sm_labels, fixed = TRUE)) else integer(0)
  pr <- if (length(excl_idx)) predict(m, newdata = newx, se.fit = TRUE, exclude = excl_idx) else predict(m, newdata = newx, se.fit = TRUE)
  newx$fit <- pr$fit; newx$se <- pr$se.fit
  newx$fit_km2 <- pmax(0, exp(newx$fit) - 1)
  newx$lwr_km2 <- pmax(0, exp(newx$fit - 1.96*newx$se) - 1)
  newx$upr_km2 <- pmax(0, exp(newx$fit + 1.96*newx$se) - 1)
  list(model = m, data = d, pred = newx)
}
gam_stats_label <- function(m, water_term_pattern = "WaterLevel"){
  sm <- summary(m)
  r2_adj <- suppressWarnings(sm$r.sq)
  dev_expl <- suppressWarnings(sm$dev.expl)
  mf <- model.frame(m); y <- model.response(mf); yhat <- fitted(m)
  r2_plain <- suppressWarnings(cor(y, yhat, use = "complete.obs")^2)
  pval <- NA_real_
  if (!is.null(sm$s.table)) {
    ridx <- grep(water_term_pattern, rownames(sm$s.table)); if (length(ridx)) pval <- sm$s.table[ridx[1], ncol(sm$s.table)]
  } else if (!is.null(sm$p.table)) {
    ridx <- grep(water_term_pattern, rownames(sm$p.table)); if (length(ridx)) {
      pcol <- grep("Pr", colnames(sm$p.table))[1]; if (length(pcol)) pval <- sm$p.table[ridx[1], pcol]
    }
  }
  paste0("R²(adj) = ", ifelse(is.finite(r2_adj), sprintf("%.2f", r2_adj), "NA"),
         " R² = ", ifelse(is.finite(r2_plain), sprintf("%.2f", r2_plain), "NA"),
         " Dev.expl = ", ifelse(is.finite(dev_expl), sprintf("%.0f%%", 100*dev_expl), "NA"),
         " p-value = ", ifelse(is.finite(pval), format.pval(pval, digits = 2, eps = 1e-3), "NA"))
}

fit_bass_50  <- fit_gam(gam_dat, "Bass",  "KDE50")
fit_snook_50 <- fit_gam(gam_dat, "Snook", "KDE50")
plot_gam_species(fit_bass_50,  "Bass",  "KDE50")
plot_gam_species(fit_snook_50, "Snook", "KDE50")
plot_gam_both(fit_bass_50, fit_snook_50, "KDE50")
fit_bass_95  <- fit_gam(gam_dat, "Bass",  "KDE95")
fit_snook_95 <- fit_gam(gam_dat, "Snook", "KDE95")
plot_gam_species(fit_bass_95,  "Bass",  "KDE95")
plot_gam_species(fit_snook_95, "Snook", "KDE95")
plot_gam_both(fit_bass_95, fit_snook_95, "KDE95")

# GAM summaries
dir.create(outdir_csv, recursive = TRUE, showWarnings = FALSE)

readr::write_lines(capture.output(summary(fit_bass_50$model)),
                   file.path(outdir_csv, "GAM_summary_Bass_KDE50.csv"))
readr::write_lines(capture.output(summary(fit_snook_50$model)),
                   file.path(outdir_csv, "GAM_summary_Snook_KDE50.csv"))
readr::write_lines(capture.output(summary(fit_bass_95$model)),
                   file.path(outdir_csv, "GAM_summary_Bass_KDE95.csv"))
readr::write_lines(capture.output(summary(fit_snook_95$model)),
                   file.path(outdir_csv, "GAM_summary_Snook_KDE95.csv"))
# ==== 13) Overlap vs. water level (GAM) ====

prep_overlap <- function(df){
  df %>% mutate(Month = as.Date(Month), method = toupper(as.character(method)),
                jaccard = as.numeric(jaccard),
                pct_overlap_bass  = as.numeric(pct_overlap_bass),
                pct_overlap_snook = as.numeric(pct_overlap_snook)) %>%
    filter(!is.na(Month))
}
ov50 <- prep_overlap(overlap_kde50) %>% left_join(wl_monthly, by = "Month")
ov95 <- prep_overlap(overlap_kde95) %>% left_join(wl_monthly, by = "Month")

fit_gam_prop <- function(dat, ycol, k_smooth = 6){
  d <- dat %>% select(WaterLevel, !!sym(ycol)) %>% rename(prop = !!sym(ycol)) %>%
    filter(is.finite(WaterLevel), is.finite(prop))
  eps <- 1e-6; d$prop <- pmin(pmax(d$prop, eps), 1-eps)
  nuniq <- length(unique(d$WaterLevel[is.finite(d$WaterLevel)]))
  k_use <- max(3, min(k_smooth, nuniq - 1))
  m <- mgcv::gam(prop ~ s(WaterLevel, k = k_use), family = mgcv::betar(link="logit"),
                 data = d, method = "REML")
  newx <- tibble(WaterLevel = seq(min(d$WaterLevel), max(d$WaterLevel), length.out = 200))
  pr <- predict(m, newdata = newx, se.fit = TRUE, type = "link")
  invlogit <- function(z) 1/(1+exp(-z))
  newx$fit <- invlogit(pr$fit)
  newx$lwr <- invlogit(pr$fit - 1.96*pr$se.fit)
  newx$upr <- invlogit(pr$fit + 1.96*pr$se.fit)
  list(model = m, data = d, pred = newx)
}
plot_overlap_species <- function(fitobj, species_name, method_name, col_line="black"){
  d <- fitobj$data; newx <- fitobj$pred
  sm <- summary(fitobj$model)
  lab <- paste0("R²(adj) = ", sprintf("%.2f", sm$r.sq),
                " Dev.expl = ", sprintf("%.0f%%", 100*sm$dev.expl))
  p <- ggplot(d, aes(WaterLevel, prop)) +
    geom_point(alpha = 0.35, size = 1.6) +
    geom_ribbon(data = newx, aes(x = WaterLevel, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, alpha = 0.20, fill = col_line) +
    geom_line(data = newx, aes(x = WaterLevel, y = fit),
              inherit.aes = FALSE, linewidth = 1.0, color = col_line) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.02, vjust = 1.5, size = 3.6) +
    labs(title = paste0(method_name, " directional overlap vs water: ", species_name),
         x = "Monthly mean water level (MO-215)", y = "Percent of species area overlapped") +
    theme_bw(base_size = 11)
  fn <- file.path(outdir_plots, paste0("GAM_OverlapDir_vs_Water_", method_name, "_", species_name, ".tiff"))
  ggsave(fn, p, device = "tiff", width = 9, height = 7, units = "in", dpi = 600, compression = "lzw")
  message("Saved: ", fn); p
}
plot_overlap_both <- function(fit_bass, fit_snook, method_name, col_bass="#377eb8", col_snook="black"){
  sm_b <- summary(fit_bass$model); sm_s <- summary(fit_snook$model)
  lab_b <- paste0("R²(adj) = ", sprintf("%.2f", sm_b$r.sq), " Dev.expl = ", sprintf("%.0f%%", 100*sm_b$dev.expl))
  lab_s <- paste0("R²(adj) = ", sprintf("%.2f", sm_s$r.sq), " Dev.expl = ", sprintf("%.0f%%", 100*sm_s$dev.expl))
  p <- ggplot() +
    geom_ribbon(data = fit_bass$pred,  aes(x=WaterLevel, ymin=lwr, ymax=upr), fill = col_bass,  alpha = 0.15, show.legend = FALSE) +
    geom_ribbon(data = fit_snook$pred, aes(x=WaterLevel, ymin=lwr, ymax=upr), fill = col_snook, alpha = 0.15, show.legend = FALSE) +
    geom_line(data = transform(fit_bass$pred,  Species="Bass"),  aes(WaterLevel, fit, colour = Species), linewidth = 1.0) +
    geom_line(data = transform(fit_snook$pred, Species="Snook"), aes(WaterLevel, fit, colour = Species), linewidth = 1.0) +
    scale_color_manual(name="Species", values=c(Bass=col_bass, Snook=col_snook)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    annotate("text", x=Inf, y=Inf, label=lab_b, hjust=1.02, vjust=2.8, size=3.6, color=col_bass) +
    annotate("text", x=Inf, y=Inf, label=lab_s, hjust=1.02, vjust=1.5, size=3.6, color=col_snook) +
    labs(title = paste0(method_name, " directional overlap vs water: Bass & Snook"),
         x = "Monthly mean water level (MO-215)", y = "Percent of species area overlapped") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  fn <- file.path(outdir_plots, paste0("GAM_OverlapDir_vs_Water_", method_name, "_both.tiff"))
  ggsave(fn, p, device = "tiff", width = 9, height = 7, units = "in", dpi = 600, compression = "lzw")
  message("Saved: ", fn); p
}
plot_overlap_jaccard <- function(fitobj, method_name){
  d <- fitobj$data; newx <- fitobj$pred; sm <- summary(fitobj$model)
  lab <- paste0("R²(adj) = ", sprintf("%.2f", sm$r.sq), " Dev.expl = ", sprintf("%.0f%%", 100*sm$dev.expl))
  p <- ggplot(d, aes(WaterLevel, prop)) +
    geom_point(alpha = 0.35, size = 1.6) +
    geom_ribbon(data = newx, aes(x=WaterLevel, ymin=lwr, ymax=upr), inherit.aes=FALSE, alpha=0.2, fill="grey70") +
    geom_line(data = newx, aes(x=WaterLevel, y=fit), inherit.aes=FALSE, linewidth=1.0) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    annotate("text", x=Inf, y=Inf, label=lab, hjust=1.02, vjust=1.5, size=3.6) +
    labs(title = paste0(method_name, " Jaccard overlap vs water"),
         x = "Monthly mean water level (MO-215)", y = "Jaccard overlap") +
    theme_bw(base_size = 11)
  fn <- file.path(outdir_plots, paste0("GAM_OverlapJaccard_vs_Water_", method_name, ".tiff"))
  ggsave(fn, p, device = "tiff", width = 9, height = 7, units = "in", dpi = 600, compression = "lzw")
  message("Saved: ", fn); p
}

# Build long frames and fit
dir50 <- bind_rows(
  ov50 %>% transmute(Month, method, Species = "Bass",  prop = pct_overlap_bass,  WaterLevel),
  ov50 %>% transmute(Month, method, Species = "Snook", prop = pct_overlap_snook, WaterLevel)
) %>% filter(is.finite(prop), is.finite(WaterLevel))
dir95 <- bind_rows(
  ov95 %>% transmute(Month, method, Species = "Bass",  prop = pct_overlap_bass,  WaterLevel),
  ov95 %>% transmute(Month, method, Species = "Snook", prop = pct_overlap_snook, WaterLevel)
) %>% filter(is.finite(prop), is.finite(WaterLevel))
jac50 <- ov50 %>% transmute(Month, method, prop = jaccard, WaterLevel) %>% filter(is.finite(prop), is.finite(WaterLevel))
jac95 <- ov95 %>% transmute(Month, method, prop = jaccard, WaterLevel) %>% filter(is.finite(prop), is.finite(WaterLevel))

fit_bass_50_dir <- fit_gam_prop(filter(dir50, Species=="Bass",  method=="KDE50"), "prop")
fit_snook_50_dir<- fit_gam_prop(filter(dir50, Species=="Snook", method=="KDE50"), "prop")
fit_jac_50 <- fit_gam_prop(filter(jac50, method=="KDE50"), "prop")
plot_overlap_species(fit_bass_50_dir,  "Bass",  "KDE50")
plot_overlap_species(fit_snook_50_dir, "Snook", "KDE50")
plot_overlap_both(fit_bass_50_dir, fit_snook_50_dir, "KDE50")
plot_overlap_jaccard(fit_jac_50, "KDE50")

fit_bass_95_dir <- fit_gam_prop(filter(dir95, Species=="Bass",  method=="KDE95"), "prop")
fit_snook_95_dir<- fit_gam_prop(filter(dir95, Species=="Snook", method=="KDE95"), "prop")
fit_jac_95 <- fit_gam_prop(filter(jac95, method=="KDE95"), "prop")
plot_overlap_species(fit_bass_95_dir,  "Bass",  "KDE95")
plot_overlap_species(fit_snook_95_dir, "Snook", "KDE95")
plot_overlap_both(fit_bass_95_dir, fit_snook_95_dir, "KDE95")
plot_overlap_jaccard(fit_jac_95, "KDE95")

save_gam_summary <- function(fit, stub){
  write_lines(capture.output(print(summary(fit$model))), file.path(outdir_csv, paste0("GAM_overlap_summary_", stub, ".txt")))
}
save_gam_summary(fit_bass_50_dir,  "BassDir_KDE50")
save_gam_summary(fit_snook_50_dir, "SnookDir_KDE50")
save_gam_summary(fit_jac_50,       "Jaccard_KDE50")
save_gam_summary(fit_bass_95_dir,  "BassDir_KDE95")
save_gam_summary(fit_snook_95_dir, "SnookDir_KDE95")
save_gam_summary(fit_jac_95,       "Jaccard_KDE95")






# ==== 14) Figures ====
# Following Code is for figure generation. This chunk was wrote so that it can 
# be ran independetly off the outputs generated above.

# Packages
pkgs <- c("dplyr","readr","lubridate","ggplot2","mgcv","scales","tidyr","patchwork")
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) install.packages(miss, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Paths
in_dir_csv   <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis"

outdir_plots <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Outputs/For Manuscript_Round1/plots"


# Water Level data and water year classification
water_csv    <- "C:/Users/atkin/Florida International University/Mangrove Coastal Fisheries - 2. Cameron Atkinson/manuscripts/Snook Home Range/R/Data for Analysis/MO-215_Jan1990-Sep2025.csv"

water_raw <- readr::read_csv(water_csv, show_col_types = FALSE) %>%
  dplyr::mutate(
    Date = lubridate::mdy(.data[["Date"]]),
    WL   = suppressWarnings(as.numeric(.data[["WaterLevel(cm)"]]))
  ) %>%
  dplyr::filter(!is.na(Date), is.finite(WL))

water_year <- function(d, start_month = 5L) {
  d <- as.Date(d)
  lubridate::year(d) + as.integer(lubridate::month(d) >= start_month)
}


# Home-range summaries (monthly)
hr_csv <- file.path(in_dir_csv, "home_ranges_monthly_summary.csv")
stopifnot(file.exists(hr_csv))
df_hr <- readr::read_csv(hr_csv, show_col_types = FALSE) %>%
  dplyr::mutate(
    Month    = as.Date(Month),
    method   = toupper(as.character(method)),
    Species  = as.factor(Species),
    area_km2 = suppressWarnings(as.numeric(area_km2))
  )

# Inter-species overlap by month (KDE50 / KDE95)
olap50_csv <- file.path(in_dir_csv, "overlap_kde50_species_month.csv")
olap95_csv <- file.path(in_dir_csv, "overlap_kde95_species_month.csv")
inter50 <- readr::read_csv(olap50_csv, show_col_types = FALSE) %>% mutate(Month = as.Date(Month))
inter95 <- readr::read_csv(olap95_csv, show_col_types = FALSE) %>% mutate(Month = as.Date(Month))

# Intra-species overlap by momnth 
bass_m50_csv  <- file.path(in_dir_csv, "intra_bass_monthly_kde50.csv")
snook_m50_csv <- file.path(in_dir_csv, "intra_snook_monthly_kde50.csv")
bass_m50  <- readr::read_csv(bass_m50_csv,  show_col_types = FALSE) %>% mutate(Month = as.Date(Month))
snook_m50 <- readr::read_csv(snook_m50_csv, show_col_types = FALSE) %>% mutate(Month = as.Date(Month))


# Monthly means for joins with HR/overlap
wl_monthly <- water_raw %>%
  dplyr::mutate(Month = lubridate::floor_date(Date, "month")) %>%
  dplyr::group_by(Month) %>%
  dplyr::summarise(WaterLevel = mean(WL, na.rm = TRUE), .groups = "drop")

# Daily water level
wl_daily0 <- water_raw %>%
  dplyr::transmute(
    Date,
    WL,
    WY = water_year(Date, start_month = 5L),
    Season = dplyr::if_else(
      Date >= lubridate::make_date(lubridate::year(Date), 5, 15) &
        Date <= lubridate::make_date(lubridate::year(Date),10,15),
      "Wet","Dry"
    ))
# ==== 15) Figure S1: Hydrograph ====
# Two panel: (a) Water level time series  (b) Marsh depth boxplots 
# Time frame (Corresponds to movment POR)
wy_min <- 2012L; wy_max <- 2025L

date_min <- as.Date(sprintf("%d-05-01", wy_min - 1))  
date_max <- as.Date(sprintf("%d-04-30", wy_max))     

breaks_wy  <- as.Date(sprintf("%d-05-01", (wy_min-1):(wy_max-1))) 
labels_wy  <- wy_min:wy_max

#GAM trendline for water level over time
wl_daily_a <- wl_daily0 %>%
  dplyr::filter(Date >= date_min, Date <= date_max) %>%
  dplyr::mutate(t_year = as.numeric(Date - date_min) / 365.25)


fit_wl_a <- mgcv::gam(WL ~ s(t_year, k = 3), data = wl_daily_a, method = "REML")
newx_a <- tibble::tibble(t_year = seq(min(wl_daily_a$t_year), max(wl_daily_a$t_year), length.out = 500)) |>
  dplyr::mutate(Date = date_min + t_year * 365.25)
pr_a <- predict(fit_wl_a, newdata = newx_a, se.fit = TRUE)
newx_a$fit <- pr_a$fit; newx_a$lwr <- pr_a$fit - 1.96*pr_a$se.fit; newx_a$upr <- pr_a$fit + 1.96*pr_a$se.fit


#  Panel a
left_pad <- lubridate::days(150)

p_figS1a <- ggplot(wl_daily_a, aes(Date, WL)) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(data = newx_a, aes(x = Date, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = 0.20) +
  geom_line(data = newx_a, aes(x = Date, y = fit), inherit.aes = FALSE, linewidth = 1.0) +
  scale_x_date(breaks = breaks_wy, labels = labels_wy, expand = expansion(mult = c(0,0))) +
  coord_cartesian(xlim = c(date_min - left_pad, date_max)) +
  labs(x = "Water Year", y = "Water Level (cm)") +
  theme_bw(base_size = 14) +
  theme(
    axis.title   = element_text(size = 22),
    axis.text    = element_text(size = 18),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text   = element_text(size = 13, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_figS1a

ggsave(file.path(outdir_plots, "MO215_WY2012-2025_Hydrograph.tiff"),
       p_figS1a, width = 25, height = 15, units = "in",
       dpi = 600, compression = "lzw")

# Panel b
wy_levels <- as.character(seq(wy_min, wy_max))

wl_fig2b <- wl_daily0 |>
  dplyr::filter(WY >= wy_min, WY <= wy_max) |>
  dplyr::mutate(
    z_wl  = as.numeric(scale(WL)),                    
    WY_f  = factor(as.character(WY), levels = wy_levels),
    Season = factor(Season, levels = c("Dry","Wet"))
  )

p_figS1b <- ggplot(wl_fig2b, aes(x = WY_f, y = z_wl)) +
  geom_boxplot(outlier.size = 0.9, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~ Season, ncol = 2) +
  scale_x_discrete(
    drop = FALSE,
    breaks = levels(wl_fig2b$WY_f),
    labels = levels(wl_fig2b$WY_f),
    expand = expansion(add = 0.6)
  ) +
  labs(x = "Water Year", y = "Water Level (z-scored)") +
  theme_bw(base_size = 13) +
  theme(
    axis.title   = element_text(size = 22),
    axis.text    = element_text(size = 18),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey95"),
    strip.text   = element_text(size = 13, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_figS1b

ggsave(file.path(outdir_plots, "M0215_WY2012-2025_ZScored.tiff"),
       p_figS1b, width = 25, height = 15, units = "in",
       dpi = 600, compression = "lzw")

# Combine & Export
p_figS1_ab <- p_figS1a / p_figS1b + patchwork::plot_annotation(tag_levels = "a")

ggsave(file.path(outdir_plots, "M0215_WY2012-2025_TwoPanel.tiff"),
       p_figS1_ab, width = 20, height = 15, units = "in",
       dpi = 600, compression = "lzw")

# ==== 16) Figure 2: HR by water year (KDE50) ====

#Set up 
method_to_plot <- "KDE50"
pal_species <- c(Bass = "#2ca25f", Snook = "#3182bd")
species_labels <- c(Bass = "Florida Bass", Snook = "Common Snook")

hr2 <- df_hr %>%
  dplyr::filter(method == method_to_plot, is.finite(area_km2)) %>%
  dplyr::mutate(WY = water_year(Month))

# Medians
med_tbl <- hr2 %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(med_area = median(area_km2, na.rm = TRUE), .groups = "drop")

# Plot
p_fig2 <- ggplot(hr2, aes(x = factor(WY), y = area_km2)) +
  geom_boxplot(aes(fill = Species), outlier.shape = NA, width = 0.7, color = "black") +
  geom_point(position = position_jitter(width = 0.25, height = 0),
             size = 0.8, alpha = 0.35, color = "black") +
  geom_hline(data = med_tbl, aes(yintercept = med_area),
             color = "black", linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE) +
  facet_wrap(~ Species, ncol = 1, labeller = as_labeller(species_labels)) +
  scale_fill_manual(values = pal_species, guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(x = "Water year", y = expression("Home range (km"^2*")")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2

ggsave(file.path(outdir_plots, "Figure2_HR_byWaterYear_KDE50.tiff"),
       p_fig2, width = 9, height = 7.5, units = "in", dpi = 600, compression = "lzw")


# ==== 17) Figure 3: GAM: HR vs monthly water (KDE50) ====

# Setup
method_to_plot <- "KDE50"
pal_species_named <- c("Florida Bass" = "#2ca25f", "Common Snook" = "#3182bd")

# GAM
hr_gam <- df_hr %>%
  dplyr::filter(method == method_to_plot, is.finite(area_km2)) %>%
  dplyr::left_join(wl_monthly, by = "Month") %>%
  dplyr::filter(is.finite(WaterLevel)) %>%
  dplyr::mutate(area_log1p = log1p(area_km2))

fit_one <- function(dat, sp, k = 6) {
  d <- dplyr::filter(dat, Species == sp); stopifnot(nrow(d) > 3)
  k_use <- max(3, min(k, dplyr::n_distinct(d$WaterLevel) - 1))
  m <- mgcv::gam(area_log1p ~ s(WaterLevel, k = k_use), data = d, method = "REML")
  grid <- tibble::tibble(WaterLevel = seq(min(d$WaterLevel), max(d$WaterLevel), length.out = 300))
  pr   <- predict(m, newdata = grid, se.fit = TRUE)
  tibble::tibble(Species = sp, WaterLevel = grid$WaterLevel,
                 fit_km2 = pmax(0, exp(pr$fit) - 1),
                 lwr_km2 = pmax(0, exp(pr$fit - 1.96*pr$se.fit) - 1),
                 upr_km2 = pmax(0, exp(pr$fit + 1.96*pr$se.fit) - 1))
}

pred_both <- dplyr::bind_rows(
  fit_one(hr_gam, sp = "Bass"),
  fit_one(hr_gam, sp = "Snook")
) %>%
  dplyr::mutate(Species = factor(Species, levels = c("Bass","Snook"),
                                 labels = c("Florida Bass","Common Snook")))

# Plot
p_fig3 <- ggplot() +
  geom_ribbon(data = pred_both,
              aes(WaterLevel, ymin = lwr_km2, ymax = upr_km2, fill = Species),
              alpha = 0.20) +
  geom_line(data = pred_both, aes(WaterLevel, fit_km2, colour = Species), linewidth = 1.0) +
  scale_color_manual(values = pal_species_named, name = NULL) +
  scale_fill_manual(values  = pal_species_named, name = NULL) +
  labs(x = "Water level (cm)", y = expression("Home range (km"^2*")")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position      = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.direction     = "vertical",
    legend.background    = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.key           = element_rect(fill = NA, color = NA),
    legend.text          = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_fig3

ggsave(file.path(outdir_plots, paste0("Figure3_HR_vs_Water_GAM_only_", method_to_plot, ".tiff")),
       p_fig3, width = 9, height = 5.5, units = "in", dpi = 600, compression = "lzw")

# ==== 18) Figure 4: Overlap by water Year (KDE50) ====

# Intra- & inter-specific comparisons
type_levels <- c(
  "Florida Bass – Florida Bass",
  "Common Snook – Common Snook",
  "Florida Bass – Common Snook"
)
pal_type <- c(
  "Florida Bass – Florida Bass"   = "#2ca25f",  
  "Common Snook – Common Snook"   = "#3182bd",  
  "Florida Bass – Common Snook"   = "#2F928E" 
)

# Ovelap info for plotting
olap5_wy <- dplyr::bind_rows(
  bass_m50  %>% dplyr::transmute(Month, type = "Florida Bass – Florida Bass", prop = mean_jaccard),
  snook_m50 %>% dplyr::transmute(Month, type = "Common Snook – Common Snook", prop = mean_jaccard),
  inter50 %>% dplyr::transmute(Month, type = "Florida Bass – Common Snook",  prop = jaccard)
) %>%
  dplyr::filter(is.finite(prop)) %>%
  dplyr::mutate(Month = as.Date(Month),
                WY    = water_year(Month),
                type  = factor(type, levels = type_levels))

# Median
med_tbl <- olap5_wy %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(mu = median(prop, na.rm = TRUE), .groups = "drop")

# Plot
p_fig4 <- ggplot(olap5_wy, aes(x = factor(WY), y = prop)) +
  geom_boxplot(aes(fill = type), outlier.shape = NA, width = 0.7, color = "black") +
  geom_point(position = position_jitter(width = 0.25, height = 0),
             size = 0.8, alpha = 0.35, color = "black") +
  geom_hline(data = med_tbl, aes(yintercept = mu),
             color = "black", linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE) +
  facet_wrap(~ type, ncol = 1, scales = "fixed") +
  scale_fill_manual(values = pal_type, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.75),
                     breaks = seq(0, 0.75, by = 0.25),
                     expand = c(0.02, 0)) +
  labs(x = "Water year", y = "Percent Overlap") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"),
        strip.text.x = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_fig4

ggsave(file.path(outdir_plots, "Figure4_Overlap_byWaterYear_KDE50_points.tiff"),
       p_fig5, width = 8.8, height = 9.2, units = "in", dpi = 600, compression = "lzw")

# ==== 19) Figure 5: GAM: Overlap vs. water level ====

# Interaction/species combos and colors
pal_type <- c(
  "Florida Bass – Florida Bass" = "#2ca25f",
  "Common Snook – Common Snook" = "#3182bd",
  "Florida Bass – Common Snook" = "#2F928E"   
)
type_levels <- names(pal_type)

# GAM
olap <- olap5_wy %>%
  dplyr::mutate(Month = as.Date(Month),
                type  = factor(type, levels = type_levels)) %>%
  dplyr::left_join(wl_monthly, by = "Month") %>%
  dplyr::filter(is.finite(prop), is.finite(WaterLevel)) %>%
  dplyr::mutate(prop = pmin(pmax(prop, 1e-6), 1 - 1e-6))


invlogit <- function(z) 1/(1 + exp(-z))
pred6 <- olap %>%
  dplyr::group_split(type) %>%
  lapply(function(d){
    m <- mgcv::gam(prop ~ s(WaterLevel, k = 6),
                   family = mgcv::betar(link = "logit"),
                   data = d, method = "REML")
    grid <- tibble::tibble(WaterLevel = seq(min(d$WaterLevel), max(d$WaterLevel), length.out = 400))
    pr   <- predict(m, newdata = grid, se.fit = TRUE, type = "link")
    tibble::tibble(
      type = d$type[1],
      WaterLevel = grid$WaterLevel,
      fit = invlogit(pr$fit),
      lwr = invlogit(pr$fit - 1.96 * pr$se.fit),
      upr = invlogit(pr$fit + 1.96 * pr$se.fit)
    )
  }) %>% dplyr::bind_rows() %>%
  dplyr::mutate(type = factor(type, levels = type_levels))

# Plot
p_fig5 <- ggplot() +
  geom_ribbon(data = pred6,
              aes(WaterLevel, ymin = lwr, ymax = upr, fill = type),
              alpha = 0.25, show.legend = FALSE) +
  geom_line(data = pred6,
            aes(WaterLevel, fit, colour = type),
            linewidth = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = pal_type) +
  scale_colour_manual(values = pal_type) +
  facet_wrap(~ type, ncol = 1, scales = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.75),
                     breaks = seq(0, 0.75, by = 0.25),
                     expand = c(0.02, 0)) +
  labs(x = "Water Level (cm)", y = "Percent Overlap") +
  theme_bw(base_size = 11) +
  theme(strip.text.x = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_fig5

ggsave(file.path(outdir_plots, "Figure5_Overlap_vs_Water_GAMonly_KDE50_colored.tiff"),
       p_fig5, width = 8.8, height = 9.2, units = "in", dpi = 600, compression = "lzw")


