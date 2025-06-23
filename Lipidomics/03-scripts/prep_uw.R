library(tidyverse)
library(here)

source(here("03-scripts", "seelipids_helpers.R"))
source(here("03-scripts", "parse_uw.R"))

# load the parsed data
file_tidy = here("02-tidydata", "uw_tidy.tsv")
# skip if parsing was just done and uwdata_long is already in the global environment
uwdata_long = read_tsv(file_tidy)

# this file annotates unique extract IDs (`eid`) with metadata (`sp`, `treatment`, etc)
file_meta = here("01-rawdata", "uw_metadata.csv")
# load metadata
metadata = file_meta %>% read_csv()

# join metadata to lipid data
uwdata = uwdata_blanked %>%
  ungroup() %>% 
  left_join(metadata, by = "eid") %>% 
  # then make compound ID a factor so it plots in a uniform order
  mutate(
    # put headgroups in the order specified in the color palette,
    # which can be changed in seelipids_helpers.R
    class = factor(class, levels = names(chroma_cl)),
    # then order things by total chain length and unsaturation
    id = factor(
      id,
      levels = cur_data() %>%
        distinct(id, class, carbon, dbonds) %>%
        arrange(class, carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>% 
  # sometimes it makes plotting easier by casting metadata
  # vars to factors with fixed order
  mutate(trt = trt %>% factor(., levels = c("C", "T"))) %>% 
  # sort the rows for easy viewing
  arrange(eid, id)

# store it
uwdata %>% write_tsv(here("02-tidydata", "uw_wmeta.tsv"))

# extract just the phospholipids
uwdata_pl = uwdata %>% 
  filter(str_detect(id, "P")) %>% 
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar))