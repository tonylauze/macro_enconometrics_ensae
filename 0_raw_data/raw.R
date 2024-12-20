### IMPORT RAW DATA ###

### Import packages 

library(ggplot2)
library(dplyr)
library(lubridate)

install.packages("ecb")
library(ecb)

install.packages("fredr")
library(fredr)
fredr_set_key("83e9bf2aec7d75da94707f7303b10042")

library(readxl)
library(readr)


### Import data 

### Import GDP 
gdp_quarterly <- fredr(
  series_id = "CLVMEURSCAB1GQEA19",
  observation_start = as.Date("2000-01-01"),
  observation_end = as.Date("2024-12-31")
)

# Garder uniquement les colonnes 'date' et 'value', et renommer 'value' en 'gdp'
gdp_cleaned <- gdp_quarterly %>%
  select(date, value) %>%
  rename(gdp = value)

# Calculer le taux de croissance en glissement annuel en utilisant la différence logarithmique
gdp_cleaned <- gdp_cleaned %>%
  arrange(date) %>%
  mutate(gdp_growth_log_yoy = (log(gdp) - log(lag(gdp, 4))) * 100) %>%   # Calcul du taux de croissance log
  filter(year(date) >= 2002)

# Tracer l'évolution du taux de croissance log en glissement annuel
ggplot(gdp_cleaned, aes(x = date, y = gdp_growth_log_yoy)) +
  geom_line(color = "blue") +
  labs(title = "Taux de Croissance en Logarithme Annuel du PIB (Euro Area)",
       x = "Date",
       y = "Taux de Croissance (%)",
       caption = "Source: FRED") +
  theme_minimal()


# GDP Deflator from FED
gdp_deflator <- fredr(
  series_id = "CP0000EZ19M086NEST",
  observation_start = as.Date("2000-01-01"),  # Adjust the start date as needed
  observation_end = as.Date("2024-12-31")    # Adjust the end date as needed
)

# GSCPI
gscpi <- read_excel("0_raw_data/gscpi_data.xls")

gscpi_cleaned <- gscpi %>%
  slice(-c(1:4)) %>% 
  select(Date, GSCPI) %>%
  rename(date = Date,
         gscpi = GSCPI)

gscpi_cleaned$date <- dmy(gscpi_cleaned$date)

# Extraire l'année et le trimestre
gscpi_quarterly <- gscpi_cleaned %>%
  mutate(year = year(date),
         quarter = quarter(date)) %>%
  group_by(year, quarter) %>%
  summarise(gscpi_quarterly = mean(gscpi, na.rm = TRUE)) %>%
  filter(year>= 2002) %>% 
  ungroup()

# Afficher les premières lignes du dataframe trimestriel pour vérifier
head(gscpi_quarterly)

ggplot(gscpi_quarterly, aes(x = as.Date(paste(year, quarter * 3, "01", sep = "-")), y = gscpi_quarterly)) +
  geom_line(color = "blue") +
  labs(title = "Évolution trimestrielle du GSCPI",
       x = "Date",
       y = "GSCPI",
       caption = "Source: Données agrégées mensuelles") +
  theme_minimal()

# HICP
hicp_data <- get_data("ICP.M.U2.N.000000.4.ANR")

hicp_clean <- hicp_data %>% 
  select(obstime, obsvalue) %>% 
  rename(date = obstime,
         hicp = obsvalue)

hicp_clean$date <- ym(hicp_clean$date)  # Convertir "YYYY-MM" en type Date

# Extraire l'année et le trimestre
hicp_quarterly <- hicp_clean %>%
  mutate(year = year(date),
         quarter = quarter(date)) %>%
  group_by(year, quarter) %>%
  summarise(hicp_quarterly = mean(hicp, na.rm = TRUE)) %>%
  filter(year >= 2002) %>% 
  ungroup()

ggplot(hicp_quarterly, aes(x = as.Date(paste(year, quarter * 3, "01", sep = "-")), y = hicp_quarterly)) +
  geom_line(color = "blue") +
  labs(title = "Évolution trimestrielle de l'indice HICP",
       x = "Date",
       y = "HICP (%)",
       caption = "Source: Données agrégées mensuelles") +
  theme_minimal()


# Deficit
deficit <- get_data("GFS.Q.N.I9.W0.S13.S1._Z.B.B9P._Z._Z._Z.XDC_R_B1GQ_CY._Z.S.V.CY._T")

deficit_clean <- deficit %>% 
  select(obstime, obsvalue) %>% 
  rename(date=obstime, 
         deficit = obsvalue) %>% 
  mutate(deficit=-deficit)

deficit_clean$date <- case_when(
  grepl("-Q1", deficit_clean$date) ~ paste0(sub("-Q1", "", deficit_clean$date), "-01-01"),
  grepl("-Q2", deficit_clean$date) ~ paste0(sub("-Q2", "", deficit_clean$date), "-04-01"),
  grepl("-Q3", deficit_clean$date) ~ paste0(sub("-Q3", "", deficit_clean$date), "-07-01"),
  grepl("-Q4", deficit_clean$date) ~ paste0(sub("-Q4", "", deficit_clean$date), "-10-01")
)


# Convertir la colonne 'date' en format Date
deficit_clean$date <- as.Date(deficit_clean$date, format = "%Y-%m-%d")

ggplot(deficit_clean, aes(x = date, y = deficit)) +
  geom_line(color = "blue") +
  labs(title = "Évolution du déficit par trimestre",
       x = "Date",
       y = "Déficit",
       caption = "Source: Données trimestrielles") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))






# Real Oil Price
real_oil_price <- read_csv("0_raw_data/DCOILBRENTEU.csv")

install.packages("quantmod")
library(quantmod)

# Télécharger les taux de change journaliers depuis FRED
getSymbols("DEXUSEU", src = "FRED", from = "2002-01-01", to = "2024-12-31")

# Voir les premières lignes de données
head(DEXUSEU)


library(dplyr)
library(zoo)

# Transformer les données en dataframe
data <- data.frame(
  date = index(DEXUSEU),
  exchange_rate = coredata(DEXUSEU)
)

# Supprimer les valeurs NA
data <- na.omit(data)

# Agréger par trimestre (moyenne trimestrielle)
quarterly_data <- data %>%
  mutate(quarter = as.yearqtr(date)) %>%
  group_by(quarter) %>%
  summarise(avg_exchange_rate = mean(DEXUSEU, na.rm = TRUE))

# Afficher les données trimestrielles
print(quarterly_data)

real_oil_price



