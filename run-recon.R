# fit RR model and plot results

library(tidyverse)
library(ggsidekick) # Sean A's helper functions for ggplot;  devtools::install_github("seananderson/ggsidekick")

# source functions
source("initRR.R")

# generate data file if inputs updated
# source("processData.R")
# processData()

# fit RR model
source("fitRR.R")
rpt <- fitRR()

# process and save time-series of spawner abundance
cdn_harvest <- read.csv("data/YkCk_Harvest_CA_Data.csv") %>%
  filter(Type == "CA_Mainstem",
         Year > 1984) %>%
  select(Year, Estimate)

border_passage <- rpt[["runSize_t"]]
cdn_er <- cdn_harvest$Estimate/(cdn_harvest$Estimate+border_passage)

CU_border_passage <- exp(rpt$lnRunSize_st)*1e-3
CU_spwn <- CU_border_passage * (1-cdn_er)

colnames(CU_spwn) <- seq(1985,2023)
rownames(CU_spwn) <- c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","YukonR.Teslinheadwaters","MiddleYukonR.andtribs.","UpperYukonR.")

CU_spwn_df <- as.data.frame(CU_spwn)
CU_spwn_df$CU <- rownames(CU_spwn)


CU_spawn_long <- CU_spwn_df %>% 
  pivot_longer(!CU, names_to = 'Year', values_to = 'spawn')

Rse <- filter(rpt$sdrpt,par=="runSize_st")
CU_spawn_long$lwr <- Rse$lCI*1e-3
CU_spawn_long$upr <- Rse$uCI*1e-3
CU_spawn_long$se <- Rse$se*1e-3


g <- ggplot(CU_spawn_long, aes(x = as.factor(Year), y = spawn, fill = CU)) + 
  geom_bar(stat = "identity", width=1) +
  scale_fill_manual(values = c("#440154FF", "#277F8EFF","#46337EFF", "#365C8DFF", "#1FA187FF", "#4AC16DFF","#9FDA3AFF", "#FDE725FF")) +
  xlab("Year") +
  ylab("Spawner abundance (000s)") +
  scale_x_discrete(breaks = c("1985","1990","1995","2000","2005","2010","2015","2020")) +
  scale_y_continuous(position = "left") +
  theme_sleek() +
  facet_wrap(~CU, ncol=3, scales = "free_y") +
  theme(strip.text = element_text(size=6),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
g

CU_spawn_long$CU <- factor(CU_spawn_long$CU, levels = c("NorthernYukonR.andtribs.","Whiteandtribs.","Pelly","Stewart","Nordenskiold","MiddleYukonR.andtribs.","YukonR.Teslinheadwaters","UpperYukonR."))

p <- ggplot(CU_spawn_long, aes(x = as.numeric(Year), y = spawn)) +
  geom_line(color="grey") +
  geom_point() +
  facet_wrap(~CU, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x="Year", y="Spawner abundance (thousands)") +
  theme_sleek() +
  theme(strip.background=element_blank(), strip.text = element_text(size = 10))
print(p)

ggsave("output/figures/CU-recon-2024.jpeg", width = 10, height = 5, units = "in")

p <- ggplot(CU_spawn_long, aes(x = as.numeric(Year), y = spawn)) +
  geom_ribbon(aes(x = as.numeric(Year), ymin = spawn_low, ymax = spawn_high), alpha=0.2) +
  geom_line(aes(x = as.numeric(Year), y = spawn), size = 0.7, alpha=0.8) + 
  geom_point(alpha=0.8) +
  facet_wrap(~CU, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x="Year", y="Spawner abundance (thousands)") +
  theme_sleek() +
  theme(strip.background=element_blank(), strip.text = element_text(size = 10))
print(p)

#### COSEWIC ----

porcupine <- read_delim('data/YkCk_PorcupineSonar_Data.csv', comment = '#') %>%
  mutate(CUID = "1209",
         Species = "Chinook",
         Region = "Yukon") %>%
  rename(Spawner.Abundance = Estimate) %>%
  select(CUID, Species, Year, Region, Spawner.Abundance)

  
yukon_COSEWIC <- CU_spawn_long %>%
  mutate(CUID = case_when(CU == "NorthernYukonR.andtribs." ~ "1207",
                          CU == "Whiteandtribs." ~ "1206",
                          CU == "Pelly" ~ "1203",
                          CU == "Stewart" ~ "1205",
                          CU == "Nordenskiold" ~ "1202",
                          CU == "YukonR.Teslinheadwaters" ~ "1212",
                          CU == "MiddleYukonR.andtribs." ~ "1204",
                          CU == "UpperYukonR." ~ "1211"),
         Year = as.numeric(Year)) %>%
  group_by(CUID, Year)%>%
  summarise(Spawner.Abundance = mean(as.numeric(spawn))*1000) %>%
  mutate(Species = "Chinook",
         Region = "Yukon") %>%
  select(CUID, Species, Year, Region, Spawner.Abundance)

yukon_COSEWIC <- rbind(yukon_COSEWIC,porcupine) %>%
  rename(LGL.counts = Spawner.Abundance)

write.csv(yukon_COSEWIC, "output/yukon_chinook_rr_escape.4Sep2024.csv")


#### CSAS ----
mssr_spwn <- CU_spawn_long %>%
  mutate(year = as.numeric(Year),
         cv = se/spawn,
         stock = CU,
         mean = as.numeric(spawn)*1000,
         obs = 1) %>%
    select(stock, year, mean, se, cv, obs)

write.csv(mssr_spwn, "output/esc-data.csv",row.names = F)
