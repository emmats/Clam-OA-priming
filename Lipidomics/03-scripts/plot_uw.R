library(tidyverse)
library(here)
library(gridExtra)

source(here("03-scripts", "seelipids_helpers.R"))
source(here("03-scripts", "prep_uw.R"))

#### PRODUCTION PLOTS

# average class composition barplots: total lipids and PLs
data_5ab = bind_rows(
  uwdata %>% 
    mutate(subs = "Total lipids"),
  uwdata_pl %>% 
    mutate(subs = "Phospholipids")
) %>% 
  # sum by class in each individual
  group_by(eid, class, trt, subs) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  # average across individuals
  group_by(class, trt, subs) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # and renormalize
  group_by(trt, subs) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
  # order factors for plotting
  mutate(
    subs  = subs %>% factor(., levels = c("Total lipids", "Phospholipids")),
    class = class %>% factor(., levels = names(chroma_cl))
  )

# panel padding in mm
pad = 2

panel_5ab = data_5ab %>% 
  gg_headgp(
    darkmode = FALSE,
    label_frac = 1,
    thres_draw = 0.01,
    # aesthetic mappings are passed straight thru
    x = trt,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~subs) +
  theme_tiny() +
  theme(
    legend.position = c(-1, -1),
    axis.line.x = element_blank(),
    panel.spacing = unit(2*pad, "mm")
  ) +
  labs(
    x = "Treatment",
    y = "Mole fraction",
    fill = "Lipid classes\n>0.01 mole frac."
  ) +
  scale_x_discrete(labels = c("Amb.", "OA")) +
  scale_y_continuous(breaks = c(0,1), labels = c(0, 1))
  #scale_y_continuous(labels = function(x){sprintf("%.1f", 2*x)})

# CHAIN AVERAGES: by class and for TL/PL
# define the class subsets to evaluate
class_foci = list(
  "TG"            = "TG",
  "Total lipids"  = uwdata$class %>% unique() %>% as.character(),
  "Phospholipids" = uwdata_pl$class %>% unique() %>% as.character(),
  # two-legged analogs only
  "PCs"           = c("PC", "O-PC", "P-PC"),
  "PEs"           = c("PE", "O-PE", "P-PE"),
  "PI"            = "PI"
) %>% 
  tibble(subs = names(.), class = .) %>% 
  unnest(class) %>% 
  mutate(subs = subs %>% factor(levels = unique(subs)))

chroma_foci = c(
  "TG"            = chroma_cl[["TG"]],
  "Total lipids"  = "black",
  "Phospholipids" = "white",
  "PCs"           = chroma_cl[["PC"]],
  "PEs"           = chroma_cl[["PE"]],
  "PI"            = chroma_cl[["PI"]]
)

olines_foci = c(
  "TG"            = "white",
  "Total lipids"  = "white",
  "Phospholipids" = "black",
  "PCs"           = "white",
  "PEs"           = "white",
  "PI"            = "white"
)

# CHAIN LENGTH PLOTS
  
data_cli = class_foci %>% 
  left_join(uwdata, by = "class", relationship = "many-to-many") %>% 
  # sum compounds by saturation level within each individual
  group_by(eid, trt, subs, carbon) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and trt
  group_by(trt, subs, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/max(frac_molar)) %>% 
  mutate(trt = c("C" = "Amb.", "T" = "OA")[trt]) %>% 
  filter(!is.na(trt))

panel_cli_tls = data_cli %>% 
  filter(subs %in% c(
    "TG"           ,
    "Total lipids" ,
    "Phospholipids"
  )) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = 0.5,
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = subs,     color = subs,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(trt),
    rows = vars(subs),
    scales = "free_y"
  ) +
  scale_fill_manual(values = chroma_foci) + scale_color_manual(values = olines_foci) +
  theme_tiny() +
  scale_y_continuous(n.breaks = 2) +
  theme(
    panel.spacing = unit(pad, "mm"),
    legend.position = "none",
    #axis.text.x = element_text(angle=90, vjust=0.5),
    axis.text.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  ) +
  guides(
    #x = guide_axis(cap = "none"),
    y = guide_axis(cap = "none")
  ) +
  coord_flip() +
  labs(
    y = "Frequency",
    x = "Total chain length"
  )

panel_cli_pls = data_cli %>% 
  filter(subs %in% c(
    "PCs"           ,
    "PEs" ,
    "PI"
  )) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = 0.5,
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = subs,     color = subs,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(trt),
    rows = vars(subs),
    scales = "free_y"
  ) +
  scale_fill_manual(values = chroma_foci) + scale_color_manual(values = olines_foci) +
  theme_tiny() +
  scale_y_continuous(n.breaks = 2) +
  theme(
    panel.spacing = unit(pad, "mm"),
    legend.position = "none",
    #axis.text.x = element_text(angle=90, vjust=0.5),
    axis.text.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  ) +
  guides(
    #x = guide_axis(cap = "none"),
    y = guide_axis(cap = "none")
  ) +
  coord_flip() +
  labs(
    y = "Frequency",
    x = "Total chain length"
  )

# UNSAT PLOTS

data_dbi = class_foci %>% 
  left_join(uwdata, by = "class", relationship = "many-to-many") %>% 
  # sum compounds by saturation level within each individual
  group_by(eid, trt, subs, dbonds) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and trt
  group_by(trt, subs, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/max(frac_molar)) %>% 
  mutate(trt = c("C" = "Amb.", "T" = "OA")[trt]) %>% 
  filter(!is.na(trt))

panel_dbi_tls = data_dbi %>% 
  filter(subs %in% c(
    "TG"           ,
    "Total lipids" ,
    "Phospholipids"
  )) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = 0.5,
    #alpha_odd = 0.5,
    x = dbonds,
    y = frac_molar,
    fill = subs,     color = subs,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(trt),
    rows = vars(subs),
    scales = "free_y"
  ) +
  scale_fill_manual(values = chroma_foci) + scale_color_manual(values = olines_foci) +
  theme_tiny() +
  scale_y_continuous(n.breaks = 2) +
  theme(
    panel.spacing = unit(pad, "mm"),
    legend.position = "none",
    #axis.text.x = element_text(angle=90, vjust=0.5),
    axis.text.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  ) +
  guides(
    #x = guide_axis(cap = "none"),
    y = guide_axis(cap = "none")
  ) +
  coord_flip() +
  labs(
    y = "Frequency",
    x = "Total unsaturations"
  )

panel_dbi_pls = data_dbi %>% 
  filter(subs %in% c(
    "PCs"           ,
    "PEs" ,
    "PI"
  )) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = 0.5,
    #alpha_odd = 0.5,
    x = dbonds,
    y = frac_molar,
    fill = subs,     color = subs,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(trt),
    rows = vars(subs),
    scales = "free_y"
  ) +
  scale_fill_manual(values = chroma_foci) + scale_color_manual(values = olines_foci) +
  theme_tiny() +
  scale_y_continuous(n.breaks = 2) +
  theme(
    panel.spacing = unit(pad, "mm"),
    legend.position = "none",
    #axis.text.x = element_text(angle=90, vjust=0.5),
    axis.text.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  ) +
  guides(
    #x = guide_axis(cap = "none"),
    y = guide_axis(cap = "none")
  ) +
  coord_flip() +
  labs(
    y = "Frequency",
    x = "Total unsaturations"
  )

### FIGURE ARRANGEMENT
arrangeGrob(
  grobs = list(
    panel_5ab,
    panel_cli_tls,
    panel_dbi_tls,
    panel_cli_pls,
    panel_dbi_pls
  ),
  widths = c(1, 0, 1, 0.1, 1.5, 0.1, 1, 0, 1),
  heights = c(0.9, 0.1),
  layout_matrix = rbind(
    c(2, NA, 3, NA, 1, NA, 4, NA, 5),
    c(2, NA, 3, NA, NA, NA, 4, NA, 5)
  )
) %>% 
  ggsave(here("04-pdf", "Fig5_20250613a.pdf"), ., width = 120, height = 60, units = "mm")

#### SCRATCH/SUPP PLOTS

## PLOT1: all lipid species, plotted individually
# not very informative :^P
uwdata %>% 
  mutate(eid = as.factor(eid)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # this can help declutter a complex plot by removing(!) low-abundance species
    # but if not left at 0, it should be set very low to avoid misrepresenting the data!
    #thres_draw = 0.002, 
    # how abundant a compound must be to get a text label
    label_frac = 0.01,
    # aesthetic mappings are passed straight thru
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~trt, scales = "free_x") +
  labs(
    title = "Headgroup composition by treat and individual",
    x = "Extract ID",
    y = "Mole fraction total lipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_alllipids.pdf"), width = 10, height = 8)
ggsave(here("05-png", "uw_alllipids.png"), width = 10, height = 8)

## PLOT2: PLs, plotted individually
uwdata_pl %>% 
  mutate(eid = as.factor(eid)) %>% 
  gg_headgp(
    darkmode = FALSE,
    # this can help declutter a complex plot by removing(!) low-abundance species
    # but if not left at 0, it should be set very low to avoid misrepresenting the data!
    #thres_draw = 0.002, 
    # how abundant a compound must be to get a text label
    label_frac = 0.01,
    # aesthetic mappings are passed straight thru
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~trt, scales = "free_x") +
  labs(
    title = "Headgroup composition by treatment and individual",
    x = "Extract ID",
    y = "Mole fraction phospholipids",
    fill = "Headgroup"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_PLs.pdf"), width = 10, height = 8)
ggsave(here("05-png", "uw_PLs.png"), width = 10, height = 8)

## PLOT3: PL acyl carbon distributions
## (change filter to look at other classes)
uwdata_pl %>% 
  # sum compounds by saturation level within each individual
  group_by(eid, class, dbonds) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and location
  group_by(trt, class, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 1.0,
    x = dbonds,
    y = frac_molar,
    fill = class,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(class),
    rows = vars(trt)
  ) +
  labs(
    title = "PL double bond distributions",
    x = "Total acyl double bonds",
    y = "Mole fraction phospholipids"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_PLacyldbonds.pdf"), width = 12, height = 4)
ggsave(here("05-png", "uw_PLacyldbonds.png"), width = 12, height = 4)

## PLOT4: PL acyl chain lengths
## (change filter to look at other classes)
uwdata_pl %>% 
  # sum compounds by saturation level within each individual
  group_by(eid, class, carbon) %>% 
  mutate(frac_molar = sum(frac_molar)) %>% 
  # then average the mole fraction within each species and location
  group_by(trt, class, carbon) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # renormalization is needed due to lipotype's precision issues
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  gg_acylch(
    darkmode = FALSE,
    meanline = TRUE,
    alpha_odd = 0.5,
    x = carbon,
    y = frac_molar,
    fill = class,
    # often I like to make the columns in the grid different
    # phospholipid classes, but here they are "species"
    cols = vars(class),
    rows = vars(trt)
  ) +
  labs(
    title = "PL total chain length distributions",
    x = "Total acyl carbons",
    y = "Mole fraction phospholipids"
  )
# save vector and raster images
ggsave(here("04-pdf", "uw_PLacylcarbon.pdf"), width = 12, height = 4)
ggsave(here("05-png", "uw_PLacylcarbon.png"), width = 12, height = 4)


