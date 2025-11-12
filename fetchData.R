# ---- paleo_streamflow_slr_auto.R ----
# Paleostreamflow reconstruction using NADA scPDSI grids as predictors
# Following methodology of Tootle et al. (2023)
# Author: (you)

# ---- Libraries ----
library(terra)
library(dataRetrieval)
library(dplyr)
library(tidyr)
library(lubridate)
library(MASS)
library(car)
library(boot)
library(qmap)

# ---- 1. Helper functions ----
get_gage_coordinates <- function(gage_id) {
  site <- readNWISsite(gage_id)
  lat <- as.numeric(site$dec_lat_va[1])
  lon <- as.numeric(site$dec_long_va[1])
  list(lat = lat, lon = lon)
}

get_gage_date_range <- function(gage_id) {
  # Query site info to get available date range
  site_info <- readNWISsite(gage_id)
  # Get the earliest and latest dates available
  # Note: readNWISsite may not always have date info, so we'll try to get it from data
  # For now, request a very wide range and let readNWISdv return what's available
  list(start = "1900-01-01", end = as.character(Sys.Date()))
}

add_lagged_predictors <- function(df, cell_cols, lags = c(1, 2)) {
  # Add lagged values for specified columns
  # df: dataframe with 'year' column and predictor columns
  # cell_cols: vector of column names to create lags for
  # lags: vector of lag periods (e.g., c(1, 2) for 1-year and 2-year lags)
  
  df <- df |> arrange(year)
  
  for (col in cell_cols) {
    if (col %in% names(df)) {
      for (lag in lags) {
        lag_col_name <- paste0(col, "_lag", lag)
        df[[lag_col_name]] <- dplyr::lag(df[[col]], n = lag)
      }
    }
  }
  
  return(df)
}

# ---- 2. User settings ----
gage_id <- "02361000"   # Choctawhatchee River near Newton, AL
radius <- 5             # degrees around gage for scPDSI extraction
# Note: Using ALL available data for calibration (no date restrictions)
sig_level <- 0.01        # correlation significance threshold

# Output files
calibration_csv <- paste0("gage", gage_id, "_NADA_calibration.csv")
reconstruction_csv <- paste0("gage", gage_id, "_NADA_reconstruction.csv")
all_cells_csv <- paste0("gage", gage_id, "_NADA_all_cells.csv")

# ---- 3. Get gage coordinates ----
coords <- get_gage_coordinates(gage_id)
lat <- coords$lat
lon <- coords$lon
cat("Gage:", gage_id, "at (", lat, ",", lon, ")\n")

# ---- 4. Download NADA scPDSI ----
nada_url <- "https://www.ncei.noaa.gov/pub/data/paleo/drought/NAmericanDroughtAtlas.v2/NADAv2-2008.nc"
nada_file <- basename(nada_url)

if (!file.exists(nada_file)) {
  message("Downloading NADA scPDSI data (~200 MB)...")
  download.file(nada_url, nada_file, mode = "wb")
}

# ---- 5. Load and crop NADA ----
r <- rast(nada_file)
bbox <- ext(lon - radius, lon + radius, lat - radius, lat + radius)
r_sub <- crop(r, bbox)
layer_names <- names(r_sub)

# Diagnostic: Check layer name format
cat("First few layer names:", head(layer_names, 5), "\n")
cat("Total layers:", length(layer_names), "\n")

# Extract years from time dimension or construct from NADA time range
# NADA v2 covers years 0-2000 CE (2001 years total, but may have fewer layers)
# Check if raster has time information
if (has.time(r_sub)) {
  time_vals <- time(r_sub)
  # Convert time to years (assuming time is in years or can be converted)
  if (inherits(time_vals, "Date") || inherits(time_vals, "POSIXt")) {
    years <- as.integer(format(time_vals, "%Y"))
  } else {
    years <- as.integer(time_vals)
  }
  cat("Extracted years from time dimension\n")
} else {
  n_layers <- nlyr(r_sub)
  start_year <- 0
  end_year <- start_year + n_layers - 1
  years <- start_year:(start_year + n_layers - 1)
  cat("Constructed years from NADA time range (0 CE to", end_year, "CE)\n")
}

cat("Extracted years range:", min(years, na.rm = TRUE), "to", max(years, na.rm = TRUE), "\n")
cat("Sample extracted years:", head(years[!is.na(years)], 10), "\n")

# ---- 6. Convert raster to tidy data frame ----
df_cells <- as.data.frame(r_sub, xy = TRUE)
df_long <- df_cells |>
  pivot_longer(cols = -c(x, y), names_to = "layer", values_to = "scpdsi") |>
  mutate(year = rep(years, times = nrow(df_cells))) |>
  drop_na(scpdsi)

cat("Years in df_long after processing:", 
    min(df_long$year, na.rm = TRUE), "to", max(df_long$year, na.rm = TRUE), "\n")
cat("Unique years in df_long:", length(unique(df_long$year)), "\n")

# ---- 7. Get USGS streamflow data ----
# Request ALL available data for the gage (not just calibration period)
# Using a very wide date range to get all available data
cat("\n=== Fetching ALL available discharge data for gage", gage_id, "===\n")
flow_raw <- readNWISdv(siteNumbers = gage_id,
                       parameterCd = "00060",
                       startDate = "1900-01-01",  # Very early date to get all available data
                       endDate = as.character(Sys.Date())) |>  # Current date to get most recent data
  renameNWISColumns()

# Diagnostic: Check Date column format
cat("\n=== Discharge Data Diagnostics ===\n")
cat("Requested date range: 1900-01-01 to", as.character(Sys.Date()), "(to get all available data)\n")
cat("Note: Using ALL available overlapping years for calibration (no date restrictions)\n")
cat("Column names:", paste(names(flow_raw), collapse = ", "), "\n")
cat("Date column class:", class(flow_raw$Date), "\n")
cat("First few Date values:", head(flow_raw$Date, 5), "\n")

# Ensure Date is in proper Date format
if (!inherits(flow_raw$Date, "Date")) {
  flow_raw$Date <- as.Date(flow_raw$Date)
  cat("Converted Date column to Date format\n")
}

cat("Raw discharge data date range:", 
    as.character(min(flow_raw$Date, na.rm = TRUE)), "to", 
    as.character(max(flow_raw$Date, na.rm = TRUE)), "\n")
cat("Total daily records:", nrow(flow_raw), "\n")
cat("Records with flow data:", sum(!is.na(flow_raw$Flow)), "\n")
cat("Records with missing flow:", sum(is.na(flow_raw$Flow)), "\n")

# Process to annual (March-October mean)
flow <- flow_raw |>
  mutate(year = year(Date),
         month = month(Date)) |>
  filter(month >= 3 & month <= 10) |>
  group_by(year) |>
  summarise(flow_cfs = mean(Flow, na.rm = TRUE),
            n_months = n(),
            date_min = min(Date),
            date_max = max(Date))

cat("\nAnnual flow summary (March-October mean):\n")
cat("Flow data years:", min(flow$year), "to", max(flow$year), "\n")
cat("Number of flow years:", nrow(flow), "\n")
first_year_idx <- which(flow$year == min(flow$year))[1]
last_year_idx <- which(flow$year == max(flow$year))[1]
cat("First year date range:", 
    as.character(flow$date_min[first_year_idx]), "to",
    as.character(flow$date_max[first_year_idx]), "\n")
cat("Last year date range:", 
    as.character(flow$date_min[last_year_idx]), "to",
    as.character(flow$date_max[last_year_idx]), "\n")
cat("Years with complete data (6+ months):", 
    sum(flow$n_months >= 6), "\n")
cat("Years with incomplete data (<6 months):", 
    sum(flow$n_months < 6), "\n")
cat("=====================================\n\n")

# ---- 8. Correlation screening (flow vs scPDSI cells) ----
flow_years <- flow$year
df_overlap <- df_long |> filter(year %in% flow_years)

cat("Overlapping years found:", length(unique(df_overlap$year)), "\n")
if (nrow(df_overlap) > 0) {
  cat("Sample overlapping years:", head(unique(df_overlap$year), 10), "\n")
}

# Identify cells by x,y coordinates
cell_coords <- df_overlap |> 
  distinct(x, y) |>
  mutate(cell_id = paste0("cell_", row_number()))

# Check if we have any cells
if (nrow(cell_coords) == 0) {
  stop("No overlapping years found between flow data and NADA data.")
}

df_overlap <- df_overlap |>
  left_join(cell_coords, by = c("x", "y"))

corrs <- data.frame(cell_id = cell_coords$cell_id, x = cell_coords$x, 
                    y = cell_coords$y, r = NA, p = NA)

for (i in seq_len(nrow(cell_coords))) {
  cdat <- df_overlap |> filter(cell_id == cell_coords$cell_id[i])
  # Select only year and scpdsi columns for merging
  cdat_subset <- cdat[, c("year", "scpdsi")]
  merged <- inner_join(flow, cdat_subset, by = "year")
  if (nrow(merged) > 10) {
    ct <- cor.test(merged$flow_cfs, merged$scpdsi)
    corrs$r[i] <- ct$estimate
    corrs$p[i] <- ct$p.value
  }
}

# Keep significant positive cells
sig_cells <- corrs |> filter(p <= sig_level & r > 0, !is.na(p))
message(nrow(sig_cells), " grid cells passed correlation screening.")

if (nrow(sig_cells) == 0)
  stop("No significant grid cells found. Try increasing radius or p threshold.")

# ---- 8b. Create spreadsheet with ALL cells and observed flow ----
# This includes all PDSI cells for ALL years in NADA (starting from earliest year)
cat("\n=== Creating spreadsheet with all PDSI cells (all years) ===\n")

# Get all PDSI data for all cells in the cropped region (all years, not just overlapping)
all_cells_all_years <- df_long |>
  left_join(cell_coords, by = c("x", "y")) |>
  filter(!is.na(cell_id)) |>  # Only cells in the cropped region
  dplyr::select(year, cell_id, scpdsi) |>
  pivot_wider(id_cols = year, names_from = cell_id, values_from = scpdsi) |>
  arrange(year)

# Get the earliest year in PDSI data
earliest_pdsi_year <- min(all_cells_all_years$year, na.rm = TRUE)
cat("Earliest PDSI year:", earliest_pdsi_year, "\n")
cat("Latest PDSI year:", max(all_cells_all_years$year, na.rm = TRUE), "\n")

# Join with observed flow (flow will be NA for years without flow data)
all_cells_with_flow <- all_cells_all_years |>
  left_join(flow |> dplyr::select(year, flow_cfs), by = "year") |>
  # Reorder columns: year and flow_cfs first, then all PDSI cells
  dplyr::select(year, flow_cfs, everything()) |>
  arrange(year)

# Get list of PDSI cell column names (all columns except year and flow_cfs)
pdsi_cols <- setdiff(names(all_cells_with_flow), c("year", "flow_cfs"))

# Add lagged predictors (t-1 and t-2 years) for all PDSI cells
cat("Adding lagged predictors (1-year and 2-year lags) for", length(pdsi_cols), "PDSI cells...\n")
all_cells_with_flow <- add_lagged_predictors(all_cells_with_flow, pdsi_cols, lags = c(1, 2))

cat("Created dataframe with", nrow(all_cells_with_flow), "years (from", 
    earliest_pdsi_year, "CE) and", length(pdsi_cols), "PDSI cells\n")
cat("Total columns:", ncol(all_cells_with_flow), "(including", length(pdsi_cols) * 2, "lagged predictors)\n")
cat("Years with observed flow:", sum(!is.na(all_cells_with_flow$flow_cfs)), "\n")
cat("Years without observed flow:", sum(is.na(all_cells_with_flow$flow_cfs)), "\n")

# Create cell metadata with coordinates and correlation info for reference
cell_metadata <- cell_coords |>
  left_join(corrs, by = c("cell_id", "x", "y")) |>
  arrange(cell_id) |>
  mutate(cell_index = row_number(),
         is_significant = ifelse(cell_id %in% sig_cells$cell_id, "Yes", "No")) |>
  dplyr::select(cell_index, cell_id, x, y, r, p, is_significant)

cat("Total grid cells:", nrow(cell_metadata), "\n")

# ---- 9. Build predictor matrix ----
sel_cell_ids <- sig_cells$cell_id
predictor_df <- df_long |>
  left_join(cell_coords, by = c("x", "y")) |>
  filter(cell_id %in% sel_cell_ids) |>
  dplyr::select(year, cell_id, scpdsi) |>
  pivot_wider(id_cols = year, names_from = cell_id, values_from = scpdsi)

df_cal <- inner_join(flow, predictor_df, by = "year") |> drop_na()

# ---- 10. Stepwise regression ----
formula_base <- as.formula(paste("flow_cfs ~", paste(sel_cell_ids, collapse = " + ")))
full_model <- lm(formula_base, data = df_cal)
step_model <- stepAIC(full_model, direction = "both", trace = FALSE)

summary(step_model)

# VIF
n_predictors <- length(coef(step_model)) - 1  # Exclude intercept
if (n_predictors >= 2) {
  cat("\nVariance Inflation Factors:\n")
  print(vif(step_model))
} else {
  cat("\nVIF not calculated: model has only", n_predictors, "predictor(s)\n")
}

dwtest <- car::durbinWatsonTest(step_model)
cat("\nDurbin-Watson test:\n")
print(dwtest)

# ---- 11. Leave-One-Out Cross-Validation ----
# Check if model has valid predictions
pred_cal <- predict(step_model, newdata = df_cal)
if (any(is.na(pred_cal))) {
  cat("Warning: Model produces NA predictions for some calibration data\n")
  cat("Skipping LOOCV\n")
} else {
  cv_results <- cv.glm(df_cal, step_model, K = nrow(df_cal))
  if (is.na(cv_results$delta[1])) {
    cat("LOOCV MSE: Could not be calculated (possibly due to model issues)\n")
  } else {
    cat("LOOCV MSE:", cv_results$delta[1], "\n")
  }
}

# ---- 12. Apply model to full scPDSI record ----
predictor_full <- df_long |>
  left_join(cell_coords, by = c("x", "y")) |>
  filter(cell_id %in% sel_cell_ids) |>
  dplyr::select(year, cell_id, scpdsi) |>
  pivot_wider(id_cols = year, names_from = cell_id, values_from = scpdsi) |>
  arrange(year)

recon <- predict(step_model, newdata = predictor_full)
reconstruction_df <- data.frame(year = predictor_full$year, recon_flow = recon)

# Check for NA predictions
n_na_recon <- sum(is.na(recon))
if (n_na_recon > 0) {
  cat("Warning:", n_na_recon, "years have NA reconstructed values (likely missing PDSI data)\n")
}

# ---- 13. Quantile Mapping Bias Correction ----
# Only use years with valid reconstructed values for calibration period
cal_years_mask <- reconstruction_df$year %in% df_cal$year
recon_cal <- recon[cal_years_mask]

# Filter out NAs for quantile mapping
valid_mask <- !is.na(recon_cal) & !is.na(df_cal$flow_cfs)
if (sum(valid_mask) < 10) {
  cat("Warning: Too few valid values for quantile mapping. Skipping bias correction.\n")
  reconstruction_df$recon_bc <- reconstruction_df$recon_flow
} else {
  tryCatch({
    fit <- fitQmapRQUANT(obs = df_cal$flow_cfs[valid_mask],
                         mod = recon_cal[valid_mask])
    # Apply bias correction only to non-NA values
    recon_bc <- reconstruction_df$recon_flow
    valid_recon_mask <- !is.na(reconstruction_df$recon_flow)
    recon_bc[valid_recon_mask] <- doQmapRQUANT(reconstruction_df$recon_flow[valid_recon_mask], fit)
    reconstruction_df$recon_bc <- recon_bc
    cat("Quantile mapping bias correction applied successfully\n")
  }, error = function(e) {
    cat("Warning: Quantile mapping failed:", e$message, "\n")
    cat("Using uncorrected reconstructed values\n")
    reconstruction_df$recon_bc <<- reconstruction_df$recon_flow
  })
}

# ---- 14. Save outputs ----
write.csv(df_cal, calibration_csv, row.names = FALSE, na = "")
write.csv(reconstruction_df, reconstruction_csv, row.names = FALSE, na = "")
write.csv(all_cells_with_flow, all_cells_csv, row.names = FALSE, na = "")

# Also save cell metadata (cell_id to coordinates mapping)
cell_metadata_csv <- paste0("gage", gage_id, "_NADA_cell_metadata.csv")
write.csv(cell_metadata, cell_metadata_csv, row.names = FALSE, na = "")

cat("Saved calibration:", calibration_csv, "\n")
cat("Saved reconstruction:", reconstruction_csv, "\n")
cat("Saved all cells data:", all_cells_csv, "\n")
cat("Saved cell metadata:", cell_metadata_csv, "\n")

# ---- 15. Plot diagnostic ----
plot(reconstruction_df$year, reconstruction_df$recon_bc, type = "l",
     main = paste("Paleostreamflow Reconstruction (USGS", gage_id, ")"),
     xlab = "Year", ylab = "Estimated Flow (cfs)", col = "blue")
points(df_cal$year, df_cal$flow_cfs, col = "red", pch = 16)
legend("topright",
       legend = c("Reconstructed (bias-corrected)", "Observed"),
       col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 16))
