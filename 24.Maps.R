# Load libraries
library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)  # for annotation_scale
library(scales)     # for alpha()

#------------------------
# 1. Data Import and Preparation
#------------------------

load("24.Maps.RData")
# Read sample data
merged_data <- read.csv("Nigeria_Guinea_Wide_ST.csv")

# Fix village naming
merged_data$Village[merged_data$Village == "Aba gboro (Ile Ife)"] <- "Abagboro"

# Clean and summarize
merged_data <- merged_data %>%
  mutate(LASV_positive = as.character(LASV_positive)) %>%
  filter(!is.na(Village) & !is.na(lattitude) & !is.na(longitude) & !is.na(State))

village_summary <- merged_data %>%
  group_by(Village, Country, lattitude, longitude, State) %>%
  summarise(
    LASV_positive_count = sum(LASV_positive == "pos"),
    total_samples = n(),
    LASV_positive_ratio = LASV_positive_count / total_samples,
    .groups = "drop"
  )

# sf objects
village_summary_sf <- st_as_sf(village_summary, coords = c("longitude", "lattitude"), crs = 4326)
village_summary_sf$ID <- 1:nrow(village_summary_sf)

# World map
all_africa <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")

# Niger & Benue Rivers
rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
niger_benue_rivers <- rivers %>% filter(grepl("Niger|Benue", name, ignore.case = TRUE))

# River labels
rivers_labels <- data.frame(
  name = c("Niger River", "Benue River"),
  longitude = c(6.3, 7.3),
  latitude = c(7.7, 7.7)
)

# City for Nigeria
cities_data <- data.frame(
  City = c("Benin City"),
  longitude = c(5.6026),
  latitude = c(6.33609)
)
cities_sf <- st_as_sf(cities_data, coords = c("longitude", "latitude"), crs = 4326)

# Faranah for Guinea
faranah_coords <- data.frame(
  City = "Faranah",
  latitude = 10.042539092871408,
  longitude = -10.740036417471334
)
Faranah_sf <- st_as_sf(faranah_coords, coords = c("longitude", "latitude"), crs = 4326)

# For pie chart placeholders
village_pie_data <- village_summary_sf  # adjust if you later add actual pie charts!

#------------------------
# 2. Guinea Plot
#------------------------

guinea_plot <- ggplot(data = all_africa) +
  geom_sf(color = "black", linewidth = 1) +
  
  # Add pie charts (currently points)
  geom_sf(data = village_pie_data, shape = 4, linewidth = .1, stroke = .5, color = "black") +
  geom_sf(data = Faranah_sf, shape = 20, linewidth = 2, stroke = 1.2, color = "grey30") +
  geom_sf(data = niger_benue_rivers, color = "navy", size = 3) +
  
  # Country labels
  geom_text(data = data.frame(label = "Guinea", x = -11.5, y = 10.5),
            aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", color = "black") +
  geom_text(data = data.frame(label = "Sierra Leone", x = -11.5, y = 9.5),
            aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", color = "black") +
  
  # Village labels
  geom_label(data = village_summary_sf, 
             aes(x = st_coordinates(geometry)[,1], 
                 y = st_coordinates(geometry)[,2], 
                 label = Village), 
             size = 3, 
             color = "black", 
             fill = alpha("lightgrey", 0.85), 
             label.size = 0.2,
             fontface = "bold",
             hjust = -.1) +
  
  # Faranah label
  geom_label_repel(data = Faranah_sf %>% filter(City == "Faranah"),
                   aes(x = st_coordinates(geometry)[,1], 
                       y = st_coordinates(geometry)[,2], 
                       label = City), 
                   size = 3, 
                   color = alpha("grey30"),
                   fill = alpha("lightgrey", 0), 
                   label.size = 0,  
                   fontface = "italic",
                   segment.color = "black", 
                   segment.size = 0.5,
                   nudge_y = 0.05,
                   nudge_x = -0.2) +
  
  annotation_scale(location = "tr", width_hint = 0.3, text_col = "black", text_cex = 0.8) +
  
  geom_text(data = rivers_labels %>% filter(name == "Niger River"), 
            aes(x = -10.5, y = 10.6, label = name), 
            size = 3, 
            color = "navy", 
            fontface = "italic") +
  
  coord_sf(xlim = c(-12, -10), ylim = c(9, 11), expand = FALSE) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

#------------------------
# 3. Nigeria Plot
#------------------------

nigeria_plot <- ggplot(data = all_africa) +
  geom_sf(color = "black", linewidth = 1) +
  
  # Pie chart placeholders
  geom_sf(data = village_pie_data, shape = 4, linewidth = .1, stroke = .5, color = "black") +
  geom_sf(data = cities_sf, shape = 20, linewidth = 2, stroke = 1.2, color = "grey30") +
  geom_sf(data = niger_benue_rivers, color = "navy", size = 3) +
  
  # Country label
  geom_text(data = data.frame(label = "Nigeria", x = 5.5, y = 8),
            aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", color = "black") +
  
  scale_size(range = c(1, 5), name = "Total Samples") +
  
  # Village labels
  geom_label(data = village_summary_sf, 
             aes(x = st_coordinates(geometry)[,1], 
                 y = st_coordinates(geometry)[,2], 
                 label = Village), 
             size = 3, 
             color = "black", 
             fill = alpha("lightgrey", 0.85), 
             label.size = 0.2,
             fontface = "bold",
             hjust = -.1,
             vjust = 1.2) +
  
  # River label
  geom_text(data = rivers_labels %>% filter(name == "Niger River"), 
            aes(x = longitude + 0.1, y = latitude, label = name), 
            size = 3, 
            color = "navy", 
            fontface = "italic") +
  
  # City label
  geom_label_repel(data = cities_sf,
                   aes(x = st_coordinates(geometry)[,1], 
                       y = st_coordinates(geometry)[,2], 
                       label = City), 
                   size = 3, 
                   color = alpha("grey30"),
                   fill = alpha("lightgrey", 0), 
                   label.size = 0,  
                   fontface = "italic",
                   segment.color = "black", 
                   segment.size = 0.5,
                   nudge_y = 0.05,
                   nudge_x = -0.2) +
  
  annotation_scale(location = "tr", width_hint = 0.3, text_col = "black", text_cex = 0.8) +
  
  coord_sf(xlim = c(4, 7), ylim = c(6, 9), expand = FALSE) +
  theme_bw() +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


#------------------------
# 4. Plot Output
#------------------------

guinea_plot
nigeria_plot

