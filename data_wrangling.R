library(tidyverse)
library(zoo)

df <- read.table("selected_indices.txt", sep=",")

################################
# Indices computed as colnames #
################################
INDICES=c("EVNsp", "MFC", "ROItotal", "ROIcover", "ACI", "BI", "roughness")
FILENAMES="filename"
headers <- c(INDICES, FILENAMES)
colnames(df) <- headers 

#################################
# Extract location and datetime #
#################################
df <- df %>%
  mutate(
    location = sub(".*proj_sound-of-norway/([^/]+)/.*", "\\1", filename),
    datetime = sub(".*proj_sound-of-norway/.*/(\\d{4}-\\d{2}-\\d{2}T\\d{2}_\\d{2}_\\d{2}).*\\.mp3", "\\1", filename),
    datetime = gsub("_", ":", datetime)  # Replace underscores with colons in the time
  )

df <- df %>%  mutate(file = sub(".*?proj_sound-of-norway/(.*)", "proj_sound-of-norway/\\1", filename))  # Extract the necessary part of the filename)

# Convert datetime to POSIXct
df$datetime <- as.POSIXct(df$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

########################
# COMBINE THE DATASET #
#######################
sub_indices <- df %>% 
  #select(file, location, datetime, NDSI, Ht, Hf, ACI, MFC, AEI, BI, ADI) %>% 
  drop_na() %>% 
  filter(location != "bugg_RPiID-1000000036c12906" & 
           location != "bugg_RPiID-10000000dda66772") 

# In correlation with bird calls
annotations = read_csv("birdnet_lite_detections-proj_sound-of-norway-yr2complete.csv")
annotations_filtered <- annotations %>% 
  select(audio_link, tags, confidence, detected_time, longitude, latitude) %>% 
  mutate(file = sub(".*?proj_sound-of-norway/(.*)", "proj_sound-of-norway/\\1", audio_link)) %>% 
  select(!audio_link)

# Get site and location
site_recorder_correspondance <- annotations %>% 
  select(location = recorder, site) %>% 
  distinct()

# Join the sites
df_annotated <- full_join(sub_indices, annotations_filtered, by = "file") %>% 
  drop_na(BI)

df_annotated <- full_join(df_annotated, site_recorder_correspondance, by = "location") %>% 
  distinct()

df_annotated <- df_annotated %>% 
  mutate(months = months(datetime)) %>% 
  mutate(months = factor(months, levels=month.name)) %>% 
  mutate(week = strftime(datetime, format="%V")) %>% 
  mutate(doy = lubridate::yday(datetime)) %>% 
  mutate(tags = replace_na(tags, "None")) %>% 
  mutate(confidence = replace_na(confidence, 1))

#############################################
# FILTER OUT SPECIES BIRDNET IS NOT GOOD ON #
#############################################

sure_species <- c("Barn Swallow", "Willow Tit", "White Wagtail", "Mew Gull", "Green Sandpiper", "Great Tit",
                  "Garden Warbler", "Willow Warbler", "European Goldfinch", "European Siskin", "Eurasian Nuthatch",
                  "European Magpie", "Eurasian Green Woodpecker", "European Greenfinch", "Yellowhammer", "Brambling",
                  "Eurasian Blue Tit", "Dunnock", "Common Chiffchaff", "Common Wood-Pigeon", "Common Swift", "Eurasian Jay",
                  "Long-eared Owl", "Hooded Crow", "Herring Gull", "Eurasian Wren", "Eurasian Robin", "Greylag Goose", 
                  "Eurasian Jackdaw", "Fieldfare", "Spotted Flycatcher", "Common Raven", "Great Spotted Woodpecker",
                  "Black Woodpecker", "Eurasian Bullfinch", "Grey Heron", "Eurasian Curlew", "Common Crane", "Barnacle Goose",
                  "Eurasian Tree Sparrow", "Eurasian Treecreeper", "Eurasian Kestrel", "Redwing", "Hawfinch", "None")

df_annotated <- df_annotated %>%
  filter(tags %in% sure_species)

####################################
# INTEGRATE NDVI VALUES TO DATASET #
####################################

ndvi <- read_csv("son_modis_ndvi_2015_2023.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  filter(format(date, "%Y") == "2022") %>% 
  mutate(week = strftime(date, format="%V")) %>% 
  mutate(ndvi = ndvi * 0.0001)

####################################################
# AVERAGE NOT GOOD, IT SHOULD BE A 2 WEEKS AVERAGE #
####################################################
ndvi_mean <- ndvi %>% 
  group_by(recorder, week) %>% 
  summarise(mean_ndvi = mean(ndvi)) %>% 
  rename(location = recorder)

ndvi_mean <- ndvi %>%
  arrange(date) %>%  # Ensure the data is sorted by date
  group_by(recorder) %>%  # Group by id and category to calculate the rolling mean within each group
  mutate(NDVI_2week_avg = rollapply(ndvi, width = 2, FUN = mean, align = 'right', fill = NA, by = 1)) %>% 
  rename(location = recorder) %>% 
  ungroup()

ggplot(ndvi_mean, aes(x = week, y=NDVI_2week_avg)) + 
  geom_point() + facet_wrap(~recorder)

# Add NDVI to the dataset
dd <- full_join(df_annotated, ndvi_mean, by = c("location", "week")) %>% 
  drop_na(NDVI_2week_avg)
