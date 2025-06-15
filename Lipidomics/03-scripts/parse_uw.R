library(tidyverse)
library(here)
library(readxl)

source(here("03-scripts", "seelipids_helpers.R"))

# this directory should contain *Lipid Species* datafiles in XLSX format
dir_data   = here("01-rawdata", "uw_core")
pat_data   = "csv" # filename pattern; UNIX glob
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T) %>% 
  # avoid MS autosave files
  .[which(!str_detect(., '~'))]

# how many datafiles did we get?
str_glue("Found {length(files_data)} datafiles:") %>% message()
files_data

# once parsed, data will be saved here for faster downstream loading
file_tidy  = here("02-tidydata", "uw_tidy.tsv")

# this block reads and binds multiple Lipid Maps Excel files
# if you have just one file, that's fine, length(files_data) will =1.
uwdata_long = files_data %>% 
  # file reader helper function melts data to tidy format
  read_uw() %>% 
  # now parse and normalize
  # mark non-detection as zero
  replace_na(list(rab = 0)) %>% 
  # see seelipids_helpers.R to understand/check outputs
  parse_uw_id()

## check parsing
#checkparse = uwdata_long %>% 
#  group_by(across(!matches(c("eid", "rab", "frac_molar")))) %>% 
#  summarize()

# a little report
str_glue("Found {length(unique(uwdata_long$id))} compounds in {length(unique(uwdata_long$eid))} samples") %>% message()

# subtract the process blank
uwdata_blanked = left_join(
  uwdata_long %>% filter(eid != "prepared_blank"),
  uwdata_long %>% filter(eid == "prepared_blank") %>% 
    select(id, rab),
  by = "id",
  suffix = c('', "_blank")
) %>% 
  mutate(
    eid = as.numeric(eid),
    rab = rab - rab_blank,
    rab = pmax(rab, 0)
  ) %>% 
  select(-rab_blank) %>% 
  # normalize to total lipid abundance
  group_by(eid) %>%
  mutate(frac_molar = rab/sum(rab)) %>%
  # remove zeroes to avoid plotting empty classes
  filter(frac_molar > 0) %>% 
  group_by(id)

# save parsed data to a new TSV
uwdata_long %>% 
  write_tsv(file_tidy)
