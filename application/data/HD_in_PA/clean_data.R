######################################
# heart disease related deaths in PA #
######################################

# death count from EDDIE system
# based on single race throughout 2013-2022

# population from CDC WONDER system
# 2018-2022 single race data from underlying cause of death
# https://wonder.cdc.gov/Deaths-by-Underlying-Cause.html
# 2013-2017 single-race population estimates 2010-2017
# https://wonder.cdc.gov/Single-Race-v2017.HTML

# EDDIE used single race
# CDC WONDER used single race since 2018
# before they have a separate table for single-race pop estimates (above)

# set up -----------------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(tigris)
library(spdep)
library(stringr)
library(ggplot2)
select <- dplyr::select




# read in ----------------------------------------------------------------------

###############
# death count #
###############

# death data was obtained from EDDIE system

death <- read.csv(
  here("application/data/death-count-hd-2013-2022.csv"),
  skip = 3, header = T, sep = ","
)

head(death)
tail(death)



#############
# pop count #
#############

# four files for one year
# 2018-2022 : All / Hispanic
# 2013-2017: All / Hispanic


# race
# 2018-22

years <- 2018:2022
temp <- NULL # all variables
pop22 <- NULL

for (year in years) {
  path <- here(paste0("application/data/pop-", year, ".txt"))
  temp <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  temp <- temp[complete.cases(temp$County.Code), ]
  temp$Year <- year
  temp[, c("Deaths", "Population", "Crude.Rate")] <- lapply(temp[, c("Deaths", "Population", "Crude.Rate")], as.numeric)
  pop22 <- rbind(pop22, temp)
}

head(pop22)
tail(pop22)

# race
# 2013-17

years <- 2013:2017
temp <- NULL # all variables
pop17 <- NULL

for (year in years) {
  path <- here(paste0("application/data/pop-", year, ".txt"))
  temp <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, na.strings = "")
  temp <- temp[complete.cases(temp$County.Code), ]
  temp <- temp[complete.cases(temp$Gender.Code), ]
  temp <- temp[complete.cases(temp$Race.Code), ]
  temp <- temp[complete.cases(temp$Ethnicity.Code), ]
  temp$Year <- year
  pop17 <- rbind(pop17, temp)
}

head(pop17)

pop17 <- pop17 %>%
  group_by(County.Code, Race, Age.Group.Code, Year, Gender) %>%
  summarize(Population = sum(Population)) # sum over ethnicity


# hispanic
# 2018-2022

years <- 2018:2022
temp <- NULL # all variables
hispop22 <- NULL

for (year in years) {
  path <- here(paste0("application/data/pop-hispanic-", year, ".txt"))
  temp <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  temp <- temp[complete.cases(temp$County.Code), ]
  temp$Year <- year
  temp[, c("Deaths", "Population", "Crude.Rate")] <- lapply(temp[, c("Deaths", "Population", "Crude.Rate")], as.numeric)
  hispop22 <- rbind(hispop22, temp)
}

head(hispop22)
tail(hispop22)



# hispanic
# 2017-2022

# hispanic pop 2013-2017 (this can be of any race)

years <- 2013:2017
temp <- NULL # all variables
hispop17 <- NULL

for (year in years) {
  path <- here(paste0("application/data/pop-hispanic-", year, ".txt"))
  temp <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  temp <- temp[complete.cases(temp$County.Code), ]
  temp$Year <- year
  temp[, c("Deaths", "Population", "Crude.Rate")] <- lapply(temp[, c("Deaths", "Population", "Crude.Rate")], as.numeric)
  hispop17 <- rbind(hispop17, temp)
}

head(hispop17)
tail(hispop17)

# compare colnames
colnames(hispop22)
colnames(hispop17)








# confirmation -----------------------------------------------------------------

if (FALSE) {
  ########
  # 2022 #
  ########
  
  # Black
  pop22 %>%
    filter(County == "Berks County, PA" &
             Single.Race.6 == "Black or African American" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2022) %>%
    select(Gender, Deaths, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Black" &
             Age == "35 to 39" &
             Year == 2022) %>%
    select(Sex, Obs_Count, Population)
  
  # White
  pop22 %>%
    filter(County == "Berks County, PA" &
             Single.Race.6 == "White" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2022) %>%
    select(Gender, Deaths, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "White" &
             Age == "35 to 39" &
             Year == 2022) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop22 %>%
    filter(County == "Philadelphia County, PA" &
             Single.Race.6 %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2022) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2022) %>%
    select(Sex, Obs_Count, Population)
  
  # Hispanic
  hispop22 %>%
    filter(County == "Berks County, PA" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2022) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Hispanic" &
             Age == "35 to 39" &
             Year == 2022) %>%
    select(Sex, Obs_Count, Population)
  
  
  
  
  
  ########
  # 2018 #
  ########
  
  # Black
  pop22 %>%
    filter(County == "Berks County, PA" &
             Single.Race.6 == "Black or African American" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2018) %>%
    select(Gender, Deaths, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Black" &
             Age == "35 to 39" &
             Year == 2018) %>%
    select(Sex, Obs_Count, Population)
  
  # White
  pop22 %>%
    filter(County == "Berks County, PA" &
             Single.Race.6 == "White" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2018) %>%
    select(Gender, Deaths, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "White" &
             Age == "35 to 39" &
             Year == 2018) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop22 %>%
    filter(County == "Philadelphia County, PA" &
             Single.Race.6 %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2018) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2018) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Hispanic
  hispop22 %>%
    filter(County == "Berks County, PA" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2018) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Hispanic" &
             Age == "35 to 39" &
             Year == 2018) %>%
    select(Sex, Obs_Count, Population)
  
  
  
  
  
  ########
  # 2017 #
  ########
  
  # Black
  pop17 %>%
    filter(County.Code == "42011" &
             Race == "Black or African American" &
             Age.Group.Code == "35-39" &
             Year == 2017) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Black" &
             Age == "35 to 39" &
             Year == 2017) %>%
    select(Sex, Obs_Count, Population)
  
  
  # White
  pop17 %>%
    filter(County.Code == "42011" &
             Race == "White" &
             Age.Group.Code == "35-39" &
             Year == 2017) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "White" &
             Age == "35 to 39" &
             Year == 2017) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop17 %>%
    filter(County.Code == "42101" &
             Race %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Age.Group.Code == "35-39" &
             Year == 2017) %>%
    select(Gender, Population) %>%
    group_by(Gender) %>%
    summarize(Population = sum(Population))
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2017) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Hispanic
  hispop17 %>%
    filter(County.Code == "42011" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2017) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Hispanic" &
             Age == "35 to 39" &
             Year == 2017) %>%
    select(Sex, Obs_Count, Population)
  
  
  
  
  
  ########
  # 2016 #
  ########
  
  # Black
  pop17 %>%
    filter(County.Code == "42101" &
             Race == "Black or African American" &
             Age.Group.Code == "35-39" &
             Year == 2016) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Black" &
             Age == "35 to 39" &
             Year == 2016) %>%
    select(Sex, Obs_Count, Population)
  
  
  # White
  pop17 %>%
    filter(County.Code == "42011" &
             Race == "White" &
             Age.Group.Code == "35-39" &
             Year == 2016) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "White" &
             Age == "35 to 39" &
             Year == 2016) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop17 %>%
    filter(County.Code == "42101" &
             Race %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Age.Group.Code == "35-39" &
             Year == 2016) %>%
    select(Gender, Population) %>%
    group_by(Gender) %>%
    summarize(Population = sum(Population))
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2016) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop17 %>%
    filter(County.Code == "42011" &
             Race %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Age.Group.Code == "35-39" &
             Year == 2016) %>%
    select(Gender, Population) %>%
    group_by(Gender) %>%
    summarize(Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2016) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Hispanic
  hispop17 %>%
    filter(County.Code == "42011" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2016) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Hispanic" &
             Age == "35 to 39" &
             Year == 2016) %>%
    select(Sex, Obs_Count, Population)
  
  
  
  ########
  # 2015 #
  ########
  
  # Black
  pop17 %>%
    filter(County.Code == "42101" &
             Race == "Black or African American" &
             Age.Group.Code == "35-39" &
             Year == 2015) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Black" &
             Age == "35 to 39" &
             Year == 2015) %>%
    select(Sex, Obs_Count, Population)
  
  
  # White
  pop17 %>%
    filter(County.Code == "42011" &
             Race == "White" &
             Age.Group.Code == "35-39" &
             Year == 2015) %>%
    select(Gender, Population)
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "White" &
             Age == "35 to 39" &
             Year == 2015) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop17 %>%
    filter(County.Code == "42011" &
             Race %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Age.Group.Code == "35-39" &
             Year == 2015) %>%
    select(Gender, Population) %>%
    group_by(Gender) %>%
    summarize(Population = sum(Population))
  
  death %>%
    filter(CountyState == "Philadelphia" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2015) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Asian/Pacific Islander
  pop17 %>%
    filter(County == "Berks County, PA" &
             Race %in% c(
               "Asian",
               "Native Hawaiian or Other Pacific Islander"
             ) &
             Age.Group.Code == "35-39" &
             Year == 2015) %>%
    select(Gender, Population) %>%
    group_by(Gender) %>%
    summarize(Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Asian/Pacific Islander" &
             Age == "35 to 39" &
             Year == 2015) %>%
    select(Sex, Obs_Count, Population)
  
  
  # Hispanic
  hispop17 %>%
    filter(County == "Berks County, PA" &
             Five.Year.Age.Groups.Code == "35-39" &
             Year == 2015) %>%
    select(Gender, Deaths, Population) %>%
    group_by(Gender) %>%
    summarize(Deaths = sum(Deaths), Population = sum(Population))
  
  death %>%
    filter(CountyState == "Berks" &
             RaceEthnicity == "Hispanic" &
             Age == "35 to 39" &
             Year == 2015) %>%
    select(Sex, Obs_Count, Population)
}




# clean data -------------------------------------------------------------------

##############
# population #
##############




# race data
# 2018-2022
pop.race22 <- pop22 %>%
  select(
    County.Code, Five.Year.Age.Groups.Code, Single.Race.6,
    Gender.Code, Population, Year
  ) %>%
  filter(Single.Race.6 %in% c(
    "Black or African American",
    "White",
    "Asian",
    "Native Hawaiian or Other Pacific Islander"
  ))  %>%
  mutate(Single.Race.6 = ifelse(Single.Race.6 == "Asian", "AAPI",
    ifelse(Single.Race.6 == "Native Hawaiian or Other Pacific Islander", "AAPI",
      ifelse(Single.Race.6 == "Black or African American", "Black",
        Single.Race.6
      )
    )
  )) %>%
  group_by(
    County.Code, Single.Race.6,
    Gender.Code, Year, Five.Year.Age.Groups.Code
  ) %>%
  summarize(Population = sum(Population, na.rm = TRUE)) %>%
  rename(
    fips = County.Code,
    race = Single.Race.6,
    sex = Gender.Code,
    age = Five.Year.Age.Groups.Code,
    year = Year,
    pop = Population
  )

pop.race22$fips <- as.character(pop.race22$fips)






# race data
# 2013-2016
pop.race17 <- pop17 %>%
  select(
    County.Code, Age.Group.Code, Race,
    Gender, Population, Year
  ) %>%
  filter(Race %in% c(
    "Asian",
    "Black or African American",
    "Native Hawaiian or Other Pacific Islander",
    "White"
  )) %>%
  mutate(Race = ifelse(Race == "Asian", "AAPI",
    ifelse(Race == "Native Hawaiian or Other Pacific Islander", "AAPI",
      ifelse(Race == "Black or African American", "Black",
        Race
      )
    )
  )) %>%
  mutate(Gender = ifelse(Gender == "Female", "F", "M")) %>%
  group_by(
    County.Code, Race,
    Gender, Year, Age.Group.Code
  ) %>%
  summarize(Population = sum(Population, na.rm = TRUE)) %>%
  rename(
    fips = County.Code,
    race = Race,
    sex = Gender,
    year = Year,
    age = Age.Group.Code,
    pop = Population
  )

pop.race17$fips <- as.character(pop.race17$fips)





# ethnicity
# 2018-2022

pop.eth22 <- hispop22 %>%
  select(
    County.Code, Five.Year.Age.Groups.Code,
    Gender.Code, Population, Year
  ) %>%
  mutate(race = "Hispanic") %>%
  group_by(County.Code, race, Gender.Code, Year, Five.Year.Age.Groups.Code) %>%
  summarize(Population = sum(Population, na.rm = TRUE)) %>%
  rename(
    fips = County.Code,
    sex = Gender.Code,
    year = Year,
    age = Five.Year.Age.Groups.Code,
    pop = Population
  )

colnames(pop.race22)
colnames(pop.race17)
colnames(pop.eth22)
pop.eth22$fips <- as.character(pop.eth22$fips)





# ethnicity
# 2013-2017

pop.eth17 <- hispop17 %>%
  select(
    County.Code, Five.Year.Age.Groups.Code,
    Gender.Code, Population, Year
  ) %>%
  mutate(race = "Hispanic") %>%
  group_by(County.Code, race, Gender.Code, Year, Five.Year.Age.Groups.Code) %>%
  summarize(Population = sum(Population, na.rm = TRUE)) %>%
  rename(
    fips = County.Code,
    sex = Gender.Code,
    year = Year,
    age = Five.Year.Age.Groups.Code,
    pop = Population
  )

pop.eth17$fips <- as.character(pop.eth17$fips)






# combine
pop.all <- rbind(
  pop.race17, pop.race22,
  pop.eth17, pop.eth22
)


# keep only age groups we need

pop.all <- pop.all %>%
  filter(age %in% c("35-39", "40-44", "45-49", "50-54"))

# pop.all <- pop.all %>%
#   filter(age %in% c("35-39", "40-44", "45-49", "50-54", "55-59", "60-64"))

pop.all%>%
  group_by(fips, sex, year) %>%
  summarize(count = n()) %>%
  filter(count != 16)

#########
# death #
#########

# death data from EDDIE (one data set)

death$FIPSCode <- paste0(42, str_pad(death$FIPSCode, width = 3, pad = 0))
death.all <- death %>%
  filter(FIPSCode != "42000") %>%
  select(FIPSCode, CountyState, Age, RaceEthnicity, Sex, Obs_Count, Year) %>%
  filter(RaceEthnicity %in% c(
    "White",
    "Black",
    "Asian/Pacific Islander",
    "Hispanic"
  )) %>%
  mutate(RaceEthnicity = ifelse(RaceEthnicity == "Asian/Pacific Islander", "AAPI", RaceEthnicity)) %>%
  mutate(Sex = ifelse(Sex == "Female", "F", "M")) %>%
  group_by(FIPSCode, CountyState, RaceEthnicity, Sex, Age, Year) %>%
  summarize(Obs_Count = sum(Obs_Count, na.rm = TRUE)) %>%
  rename(
    fips = FIPSCode,
    name = CountyState,
    race = RaceEthnicity,
    sex = Sex,
    year = Year,
    age = Age,
    death = Obs_Count
  )

death.all$age <- gsub(" to ", "-", death.all$age)

death.all <- death.all %>%
  filter(age %in% unique(pop.all$age))

dim(pop.all)
dim(death.all)
# same number of rows




# merge data -------------------------------------------------------------------

pa.data <- left_join(pop.all, death.all) %>%
  group_by(fips, race, sex, year) %>%
  summarize(pop = sum(pop, na.rm = TRUE),
    death = sum(death, na.rm = TRUE)
  )

dim(pa.data)

pa.data$rate <- pa.data$death / pa.data$pop * 10^5

pa.data %>%
  group_by(fips, race, sex, year) %>%
  summarize(count = n()) %>%
  filter(count != 1)

pa.data %>%
  group_by(fips, race, sex) %>%
  summarize(count = n()) %>%
  filter(count != 10)

pa.data %>%
  group_by(fips, race, year) %>%
  summarize(count = n()) %>%
  filter(count != 2)

pa.data %>%
  group_by(fips, sex, year) %>%
  summarize(count = n()) %>%
  filter(count != 4)

pa.data %>%
  group_by(race, sex, year) %>%
  summarize(count = n()) %>%
  filter(count != 67)

# # download shapefiles ----------------------------------------------------------
# 
# 
# pa <- counties(state = "PA", cb = TRUE, class = "sf")
# 
# 
# pa.sf <- pa %>%
#   arrange(GEOID) %>%
#   select(GEOID, geometry) %>%
#   rename(fips = GEOID)
# 
# pa.info <- pa %>%
#   arrange(GEOID) %>%
#   select(GEOID, NAME)
# 
# Ns <- dim(pa.sf)[1]
# 
# # nb
# nb <- poly2nb(pa.sf)
# nb.info <- nb2WB(nb)
# num <- nb.info$num
# adj <- nb.info$adj
# weights <- nb.info$weights
# 
# # D (diagonal matrix with m_i)
# D <- diag(num)
# 
# # W (weight matrix or adjacency matrix)
# W <- matrix(0, nrow = Ns, ncol = Ns)
# 
# for (i in 1:Ns) {
#   W[i, adj[(max(sum(num[(i - 1):0]), 0) + 1):sum(num[0:i])]] <- 1
# }
# 
# # neigh
# neigh <- list()
# count <- 0
# for (i in 1:Ns) {
#   if (num[i] > 0) {
#     neigh[[i]] <- adj[count + 1:num[i]]
#     count <- count + num[i]
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# # maps -------------------------------------------------------------------------
# 
# pa.sf.data <- left_join(pa.sf, pa.data, by = "fips")
# 
# pa.sf.data %>%
#   filter(sex == "M") %>%
#   filter(death > 5) %>%
#   ggplot(aes(fill = rate)) +
#   geom_sf(linewidth = 0.1) +
#   scale_fill_viridis_c(
#     direction = -1,
#     na.value = "lightgrey"
#   ) +
#   facet_grid(race ~ year) +
#   labs(
#     title = "Heart Disease Mortality Rates",
#     subtitle = "Males",
#     fill = "Death Per 100,000"
#   ) +
#   theme_light()
# 
# 
# pa.sf.data %>%
#   filter(sex == "M") %>%
#   ggplot(aes(fill = death)) +
#   geom_sf(linewidth = 0.1) +
#   scale_fill_viridis_c(
#     direction = -1,
#     na.value = "lightgrey"
#   ) +
#   facet_grid(race ~ year) +
#   labs(
#     title = "Heart Disease Mortality Cases",
#     subtitle = "Males",
#     fill = "Count"
#   ) +
#   theme_light()
# 
# 
# pa.data %>%
#   filter(sex == "M") %>%
#   ggplot(aes(x = death)) +
#   geom_histogram() +
#   xlim(c(0, 30)) +
#   ylim(c(0, 15)) +
#   facet_grid(year ~ race)
#   
# pa.data %>%
#   filter(sex == "F") %>%
#   ggplot(aes(x = death)) +
#   geom_histogram() +
#   xlim(c(0, 30)) +
#   ylim(c(0, 15)) +
#   facet_grid(year ~ race)
# 


# reshape data -----------------------------------------------------------------

# create sex-race variable so that I can convert them into 3 dim array
pa.data <- pa.data %>%
  mutate(sex_race = paste0(sex, "-", race))

# each label and length
Ns <- length(unique(pa.data$fips))
label_county <- unique(pa.data$name)
Nrs <- length(unique(pa.data$sex_race))
label_sexrace <- unique(pa.data$sex_race)
Ng <- length(unique(pa.data$year))
label_year <- unique(pa.data$year)

# create Y (death)

Y <- array(NA, dim = c(Ns, Ng, Nrs))

# fill the 3D array
for (i in 1:Ns) {
  for (j in 1:Ng) {
    for (k in 1:Nrs) {
      # find the subset of the data frame that matches the current fips, year, and sex_race
      subset_data <- subset(
        pa.data,
        fips == unique(pa.data$fips)[i] &
          year == unique(pa.data$year)[j] &
          sex_race == unique(pa.data$sex_race)[k]
      )
      if (nrow(subset_data) > 0) {
        # Use the 'death' column for the values in the array
        Y[i, j, k] <- subset_data$death
      }
    }
  }
}

# create n (population)

n <- array(NA, dim = c(Ns, Ng, Nrs))

# fill the 3D array
for (i in 1:Ns) {
  for (j in 1:Ng) {
    for (k in 1:Nrs) {
      # find the subset of the data frame that matches the current fips, year, and sex_race
      subset_data <- subset(
        pa.data,
        fips == unique(pa.data$fips)[i] &
          year == unique(pa.data$year)[j] &
          sex_race == unique(pa.data$sex_race)[k]
      )
      if (nrow(subset_data) > 0) {
        # Use the 'death' column for the values in the array
        n[i, j, k] <- subset_data$pop
      }
    }
  }
}


# label both Y and n

dimnames(Y) <- list(fips = label_county, year = label_year, sex_race = label_sexrace)
dimnames(n) <- list(fips = label_county, year = label_year, sex_race = label_sexrace)


# save data --------------------------------------------------------------------

save(pa.data, Y, n,
  Ns, Ng, Nrs,
  label_county, label_year, label_sexrace,
  file = here("application/data/hd-35-54.Rdata")
)
