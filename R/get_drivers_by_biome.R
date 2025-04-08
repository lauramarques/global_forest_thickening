get_drivers_by_biome <- function(v_biomes_forests){
  
  # read aridity index data
  rasta <- rast("~/data/aridityindex_zomer_2022/Global-AI_ET0_v3_annual/ai_v3_yr.tif")
  
  # values are provided as integers, multiplied by 1e-4
  values(rasta) <- values(rasta) * 1e-4
  
  # aggregate to 0.5 deg
  rasta_agg <- aggregate(
    rasta,
    fact = res(rasta)[1]^-1*0.5,
    fun = "mean"
  )
  
  # Crop and mask the raster to the polygon
  r_masked <- mask(rasta_agg, v_biomes_forests)
  
  # Compute cell areas in hectares
  area_raster <- cellSize(r_masked, unit="ha")
  
  # Extract cell areas at these points
  areas <- extract(area_raster, v_biomes_forests)
  
  # Convert the raster to a data frame with coordinates
  r_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE) |>
    rename(ai = awi_pm_sr_yr,
           lon = x,
           lat = y) 
  
  points <- vect(r_df, geom=c("lon", "lat"), crs=crs(r_masked))
  
  # Extract cell areas at these points
  areas <- extract(area_raster, points)
  
  r_df <- r_df |>
    mutate(area_ha = areas[,2]) |>
    relocate(area_ha, .after = lat)
  
  # reset row numbers
  rownames(r_df) <- NULL
  
  r_df <- r_df |>
    rownames_to_column("sitename")
  
  siteinfo <- r_df |>
    select(sitename, lon, lat) |>
    as_tibble()
  
  # PBR
  settings_gsde <- list(varnam = c("PBR"), layer = 1:3)
  
  df_gsde <- ingest(
    siteinfo,
    source    = "gsde",
    settings  = settings_gsde,
    dir       = "~/data/soil/shangguan"
  )
  
  df_gsde <- df_gsde |>
    unnest(cols = data) 
  
  # ORGC
  settings_wise <- get_settings_wise(varnam = c("ORGC"), layer = 1:3)
  
  df_wise <- ingest(
    siteinfo,
    source    = "wise",
    settings  = settings_wise,
    dir       = "~/data/soil/wise"
  )
  
  df_wise <- df_wise |>
    unnest(cols = data)
  
  # Ndep
  siteinfo <- siteinfo |>
    mutate(year_start = 1970,
           year_end = 2009)
  
  df_ndep <- ingest(
    siteinfo,
    source    = "ndep",
    timescale = "y",
    dir       = "~/data/ndep_lamarque/",
    verbose   = FALSE
  )
  
  df_ndep <- df_ndep |>
    unnest(cols = data) |>
    group_by(sitename) |>
    # Calculate the mean for the last 40 years available (1970 to 2009)
    summarise(noy = mean(noy, na.rm=T), nhx = mean(nhx, na.rm=T))
  
  # Join all data
  data_agg <- r_df |> 
    left_join(df_gsde) |>
    left_join(df_wise) |>
    left_join(df_ndep) |>
    mutate(ndep = noy+nhx) |>
    rename(plotID = sitename)
  
  return(data_agg)
}
