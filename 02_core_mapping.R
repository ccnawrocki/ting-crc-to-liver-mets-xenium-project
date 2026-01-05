rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

# Outs directory
outs_dir <- "02_core_mapping_outs"
dir.create(path = "02_core_mapping_outs")

# Need the core to patient map (CORE.NOTES columns is the patient ID)
coremap_defacto <- openxlsx::read.xlsx(xlsxFile = "../../ting-crc-to-liver-mets-data/maps.xlsx", sheet = 7, cols = 9:15)


# ## prim1 -----------------------------------------------------------------------
# 
# prim1 <- data.table::fread("../../ting-crc-to-liver-mets-data/prim1/stow/cells.csv.gz") |> 
#   as.data.frame()
# prim1$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = prim1$segmentation_method)
# rownames(prim1) <- prim1$cell_id
# prim1$run <- 2
# prim1$slide <- "prim1"
# prim1$sample_type <- "primary"
# 
# prim1_cores <- sf::read_sf("../../ting-crc-to-liver-mets-data/prim1_cores.geojson") |> dplyr::select(name, geometry)
# sf::st_crs(prim1_cores) <- NA
# 
# library(ggplot2)
# ggplot() + 
#   geom_sf(data = prim1_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") + 
#   scattermore::geom_scattermore(data = prim1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) + 
#   theme_void()
# ggsave(filename = file.path(outs_dir, "prim1_cores.pdf"), width = 6, height = 8)
# 
# cell_pts <- sf::st_as_sf(x = prim1 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
# cellstatus <- sf::st_within(x = cell_pts, y = prim1_cores, sparse = T)
# cellstatus <- as.vector(cellstatus)
# cellstatus[(lengths(cellstatus) == 0)] <- NA
# coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(prim1_cores), to = prim1_cores$name)
# prim1$core <- coreids
# 
# ggplot() + 
#   geom_sf(data = prim1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
#   scattermore::geom_scattermore(data = prim1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
#   geom_text(data = prim1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500), 
#             mapping = aes(x = X, y = Y, label = core)) +
#   theme_void() + 
#   ggprism::scale_color_prism() +
#   Seurat::NoLegend()
# ggsave(filename = file.path(outs_dir, "prim1_cores_labeled.pdf"), width = 6, height = 8)
# 
# prim1 <- prim1[!is.na(prim1$core),]
# prim1$patient <- plyr::mapvalues(x = paste(prim1$slide, prim1$core), from = paste(coremap_defacto$`TMA.#`, coremap_defacto$POS), to = coremap_defacto$CORE.NOTES)
# 
# ggplot() + 
#   geom_sf(data = prim1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
#   scattermore::geom_scattermore(data = prim1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
#   geom_text(data = prim1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)), 
#             mapping = aes(x = X, y = Y, label = patient)) +
#   theme_void() + 
#   ggprism::scale_color_prism() +
#   Seurat::NoLegend()
# ggsave(filename = file.path(outs_dir, "prim1_cores_by_patient.pdf"), width = 6, height = 8)


## prim2 -----------------------------------------------------------------------

prim2 <- data.table::fread("../../ting-crc-to-liver-mets-data/prim2/stow/cells.csv.gz") |> 
  as.data.frame()
prim2$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = prim2$segmentation_method)
rownames(prim2) <- prim2$cell_id
prim2$run <- 2
prim2$slide <- "prim2"
prim2$sample_type <- "primary"

prim2_cores <- sf::read_sf("../../ting-crc-to-liver-mets-data/prim2_cores.geojson") |> dplyr::select(name, geometry)
sf::st_crs(prim2_cores) <- NA

library(ggplot2)
ggplot() + 
  geom_sf(data = prim2_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") + 
  scattermore::geom_scattermore(data = prim2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) + 
  theme_void()
ggsave(filename = file.path(outs_dir, "prim2_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = prim2 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = prim2_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(prim2_cores), to = prim2_cores$name)
prim2$core <- coreids

ggplot() + 
  geom_sf(data = prim2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
  scattermore::geom_scattermore(data = prim2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = prim2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500), 
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() + 
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "prim2_cores_labeled.pdf"), width = 6, height = 8)

prim2 <- prim2[!is.na(prim2$core),]
prim2$patient <- plyr::mapvalues(x = paste(prim2$slide, prim2$core), from = paste(coremap_defacto$`TMA.#`, coremap_defacto$POS), to = coremap_defacto$CORE.NOTES)

ggplot() + 
  geom_sf(data = prim2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
  scattermore::geom_scattermore(data = prim2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = prim2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)), 
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() + 
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "prim2_cores_by_patient.pdf"), width = 6, height = 8)


# ## mets1 -----------------------------------------------------------------------
# 
# mets1 <- data.table::fread("../../ting-crc-to-liver-mets-data/mets1/stow/cells.csv.gz") |> 
#   as.data.frame()
# mets1$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = mets1$segmentation_method)
# rownames(mets1) <- mets1$cell_id
# mets1$run <- 2
# mets1$slide <- "mets1"
# mets1$sample_type <- "mets"
# 
# mets1_cores <- sf::read_sf("../../ting-crc-to-liver-mets-data/mets1_cores.geojson") |> dplyr::select(name, geometry)
# sf::st_crs(mets1_cores) <- NA
# 
# library(ggplot2)
# ggplot() + 
#   geom_sf(data = mets1_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") + 
#   scattermore::geom_scattermore(data = mets1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) + 
#   theme_void()
# ggsave(filename = file.path(outs_dir, "mets1_cores.pdf"), width = 6, height = 8)
# 
# cell_pts <- sf::st_as_sf(x = mets1 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
# cellstatus <- sf::st_within(x = cell_pts, y = mets1_cores, sparse = T)
# cellstatus <- as.vector(cellstatus)
# cellstatus[(lengths(cellstatus) == 0)] <- NA
# coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(mets1_cores), to = mets1_cores$name)
# mets1$core <- coreids
# 
# ggplot() + 
#   geom_sf(data = mets1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
#   scattermore::geom_scattermore(data = mets1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
#   geom_text(data = mets1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500), 
#             mapping = aes(x = X, y = Y, label = core)) +
#   theme_void() + 
#   ggprism::scale_color_prism() +
#   Seurat::NoLegend()
# ggsave(filename = file.path(outs_dir, "mets1_cores_labeled.pdf"), width = 6, height = 8)
# 
# mets1 <- mets1[!is.na(mets1$core),]
# mets1$patient <- plyr::mapvalues(x = paste(mets1$slide, mets1$core), from = paste(coremap_defacto$`TMA.#`, coremap_defacto$POS), to = coremap_defacto$CORE.NOTES)
# 
# ggplot() + 
#   geom_sf(data = mets1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
#   scattermore::geom_scattermore(data = mets1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
#   geom_text(data = mets1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)), 
#             mapping = aes(x = X, y = Y, label = patient)) +
#   theme_void() + 
#   ggprism::scale_color_prism() +
#   Seurat::NoLegend()
# ggsave(filename = file.path(outs_dir, "mets1_cores_by_patient.pdf"), width = 6, height = 8)


## mets2 -----------------------------------------------------------------------

mets2 <- data.table::fread("../../ting-crc-to-liver-mets-data/mets2/stow/cells.csv.gz") |> 
  as.data.frame()
mets2$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = mets2$segmentation_method)
rownames(mets2) <- mets2$cell_id
mets2$run <- 2
mets2$slide <- "mets2"
mets2$sample_type <- "mets"

mets2_cores <- sf::read_sf("../../ting-crc-to-liver-mets-data/mets2_cores.geojson") |> dplyr::select(name, geometry)
sf::st_crs(mets2_cores) <- NA

library(ggplot2)
ggplot() + 
  geom_sf(data = mets2_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") + 
  scattermore::geom_scattermore(data = mets2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) + 
  theme_void()
ggsave(filename = file.path(outs_dir, "mets2_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = mets2 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = mets2_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(mets2_cores), to = mets2_cores$name)
mets2$core <- coreids

ggplot() + 
  geom_sf(data = mets2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
  scattermore::geom_scattermore(data = mets2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = mets2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500), 
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() + 
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "mets2_cores_labeled.pdf"), width = 6, height = 8)

mets2 <- mets2[!is.na(mets2$core),]
mets2$patient <- plyr::mapvalues(x = paste(mets2$slide, mets2$core), from = paste(coremap_defacto$`TMA.#`, coremap_defacto$POS), to = coremap_defacto$CORE.NOTES)

ggplot() + 
  geom_sf(data = mets2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") + 
  scattermore::geom_scattermore(data = mets2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = mets2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)), 
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() + 
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "mets2_cores_by_patient.pdf"), width = 6, height = 8)


## Saving ----------------------------------------------------------------------
qs2::qs_save(
  object = list(
  # "prim1" = prim1, 
  "prim2" = prim2, 
  # "mets1" = mets1, 
  "mets2" = mets2
), 
file = file.path(outs_dir, "core-mapped_metadata.qs2")
)

