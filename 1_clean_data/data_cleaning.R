### IMPORT DATA ###

# Load packages -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)

install.packages("ecb")
library(ecb)

install.packages("fredr")
library(fredr)
fredr_set_key("83e9bf2aec7d75da94707f7303b10042")

library(readxl)
library(readr)

install.packages("quantmod")
library(quantmod)
library(dplyr)
library(zoo)

# Import and clean data ---------------------------------------------------

### Import GDP ###
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
  filter(year(date) >= 2002 & year(date) < 2024) %>% 
  rename(gdp_growth_log = gdp_growth_log_yoy) %>% 
  mutate(date = as.yearqtr(date))

# Tracer l'évolution du taux de croissance log en glissement annuel
ggplot(gdp_cleaned, aes(x = date, y = gdp_growth_log)) +
  geom_line(color = "blue") +
  labs(title = "Taux de Croissance en Logarithme Annuel du PIB (Euro Area)",
       x = "Date",
       y = "Taux de Croissance (%)",
       caption = "Source: FRED") +
  theme_minimal()

gdp_cleaned

### GDP Deflator from FED ###
gdp_deflator <- fredr(
  series_id = "CP0000EZ19M086NEST",
  observation_start = as.Date("2000-01-01"),  # Adjust the start date as needed
  observation_end = as.Date("2024-12-31")    # Adjust the end date as needed
) %>% 
  select(date, value)
gdp_deflator

gdp_deflator_quarterly <- gdp_deflator %>%
  # Ajouter une colonne pour l'année et le trimestre
  mutate(year = year(date), quarter = quarter(date)) %>%
  # Filtrer les données pour ne garder que les années entre 2002 et 2023
  filter(year >= 2002 & year <= 2023) %>%
  # Regrouper par année et trimestre, puis calculer la moyenne
  group_by(year, quarter) %>%
  summarise(avg_gdp_deflator = mean(value, na.rm = TRUE)) %>%
  # Créer une colonne date pour l'affichage
  mutate(date = as.Date(paste(year, quarter * 3, "01", sep = "-"))) %>%
  select(date, avg_gdp_deflator) %>% 
  rename(gdp_deflator_fed = avg_gdp_deflator) %>% 
  mutate(date = as.yearqtr(date))


### GSCPI ###
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
  filter(year>= 2002 & year<2024) %>% 
  ungroup() %>% 
  rename(gscpi = gscpi_quarterly) %>% 
  mutate(date = as.yearqtr(paste(year, quarter, "01", sep = "-")))  # Créer la colonne 'date' en 'yearqtr'


# Afficher les premières lignes du dataframe trimestriel pour vérifier
head(gscpi_quarterly)

ggplot(gscpi_quarterly, aes(x = as.Date(paste(year, quarter * 3, "01", sep = "-")), y = gscpi)) +
  geom_line(color = "blue") +
  labs(title = "Évolution trimestrielle du GSCPI",
       x = "Date",
       y = "GSCPI",
       caption = "Source: Données agrégées mensuelles") +
  theme_minimal()

### HICP ###
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
  filter(year >= 2002 & year<2024) %>% 
  ungroup() %>% 
  rename(hicp_ecb = hicp_quarterly) %>% 
  mutate(date = as.yearqtr(paste(year, quarter, "01", sep = "-")))  # Créer la colonne 'date' en 'yearqtr'


ggplot(hicp_quarterly, aes(x = as.Date(paste(year, quarter * 3, "01", sep = "-")), y = hicp_ecb)) +
  geom_line(color = "blue") +
  labs(title = "Évolution trimestrielle de l'indice HICP",
       x = "Date",
       y = "HICP (%)",
       caption = "Source: Données agrégées mensuelles") +
  theme_minimal()

hicp_quarterly

### Deficit ###
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

deficit_clean <- deficit_clean %>% 
  filter(year(date)>2001 & year(date)<2024) %>% 
  mutate(date = as.yearqtr(date))


ggplot(deficit_clean, aes(x = date, y = deficit)) +
  geom_line(color = "blue") +
  labs(title = "Évolution du déficit par trimestre",
       x = "Date",
       y = "Déficit",
       caption = "Source: Données trimestrielles") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


### Real Oil Price ###
real_oil_price <- read_csv("0_raw_data/DCOILBRENTEU.csv") %>% 
  rename(date = observation_date,
         real_price = DCOILBRENTEU)

real_oil_price <- real_oil_price %>%
  mutate(quarter = as.yearqtr(date))

# Télécharger les taux de change journaliers depuis FRED
getSymbols("DEXUSEU", src = "FRED", from = "2002-01-01", to = "2024-12-31")

# Voir les premières lignes de données
head(DEXUSEU)

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

oil_price_euro <- real_oil_price %>%
  inner_join(quarterly_data, by = "quarter") %>%
  mutate(
    real_price_euro = real_price / avg_exchange_rate, # Conversion en euros
    year = year(date) # Extraction de l'année
  ) %>%
  filter(year >= 2002 & year <= 2023) %>% # Filtrage des années
  select(date, quarter, real_price_euro, year)

deflated_oil_prices <- oil_price_euro %>%
  mutate(
    quarter_num = as.numeric(format(quarter, "%q")) # Conversion du trimestre
  ) %>%
  left_join(hicp_quarterly, by = c("year", "quarter_num" = "quarter")) %>%
  mutate(
    deflated_price_euro = real_price_euro / (1 + hicp_ecb / 100) # Calcul du prix déflaté
  ) %>% 
  mutate(
    log_deflated_price = log(deflated_price_euro) * 100 # Logarithme et mise à l'échelle
  ) %>% 
  select(date=quarter, log_deflated_price)

# Print the resulting dataframe
ggplot(deflated_oil_prices, aes(x = date, y = log_deflated_price)) +
  geom_line(color = "blue") +
  labs(title = "Log-Transformed Deflated Oil Prices",
       x = "Date",
       y = "Log(Deflated Price) × 100",
       caption = "Source: Deflated prices based on HICP data") +
  theme_minimal()

### SHADOW RATE ###
shadow_rate <- read_excel("0_raw_data/shadow_krippner.xlsx")

shadow_rate_clean <- shadow_rate %>%
  # Renommer les colonnes pour plus de clarté
  rename(date = `34758`, shadow_rate = `5.3346193777554403`) %>%
  # Convertir la colonne date en format Date
  mutate(date = as.Date(date)) %>%
  # Ajouter une colonne pour l'année et le trimestre
  mutate(year = year(date), quarter = quarter(date)) %>%
  # Filtrer les données pour ne garder que les années entre 2002 et 2023
  filter(year >= 2002 & year <= 2023) %>%
  # Regrouper par année et trimestre, puis calculer la moyenne
  group_by(year, quarter) %>%
  summarise(avg_shadow_rate = mean(shadow_rate, na.rm = TRUE)) %>%
  # Créer une colonne date pour l'affichage
  mutate(date = as.Date(paste(year, quarter * 3, "01", sep = "-"))) %>%
  select(date, avg_shadow_rate) %>% 
  mutate(date = as.yearqtr(date))

ggplot(shadow_rate_clean, aes(x = date, y = avg_shadow_rate)) +
  geom_line(color = "blue") +
  labs(title = "Taux moyen trimestriel de Shadow Rate (2002-2023)",
       x = "Date",
       y = "Shadow Rate",
       caption = "Source: Données agrégées trimestrielles") +
  theme_minimal()


# FINAL SERIES ------------------------------------------------------------

data_clean <- shadow_rate_clean %>%
  left_join(deflated_oil_prices, by = "date") %>%
  left_join(deficit_clean, by = "date") %>%
  left_join(hicp_quarterly, by = "date") %>%
  left_join(gscpi_quarterly, by = "date") %>%
  left_join(gdp_deflator_quarterly, by = "date") %>%
  left_join(gdp_cleaned, by = "date") %>%
  select(date, avg_shadow_rate, log_deflated_price, deficit, hicp_ecb, gscpi, gdp_deflator_fed, gdp_growth_log)

saveRDS(data_clean, "~/work/macro_enconometrics_ensae/1_clean_data/data_clean.rds")


# Vérifier les 5 premières lignes du dataframe final
head(data_clean)

data_long <- data_clean %>%
  pivot_longer(
    cols = -c(date, gdp_deflator_fed), # Exclude 'deflator_fed' column
    names_to = "variable", # Name for the new variable column
    values_to = "value" # Name for the new value column
  )

# Generate Plot
data_long$variable <- factor(data_long$variable, levels = c(
  "gdp_growth_log", "hicp_ecb", "gscpi", 
  "avg_shadow_rate", "log_deflated_price", "deficit"
))

ggplot(data_long, aes(x = date, y = value)) +
  geom_line(color = "black", size = 0.8) + # Main data line in black
  facet_wrap(~ variable, scales = "free_y", ncol = 3, 
             labeller = as_labeller(c(
               "gdp_growth_log" = "GDP",
               "hicp_ecb" = "Inflation",
               "gscpi" = "GSCPI",
               "avg_shadow_rate" = "Shadow Rate",
               "log_deflated_price" = "Real Oil Price",
               "deficit" = "Deficit"
             ))) + # Custom labels for variables
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold"), # Panel titles
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels for better readability
    panel.spacing = unit(1, "lines") # Increase space between panels
  ) +
  labs(
    title = "Figure A.1: Data (without Deflator Fed and Trend Line)",
    x = NULL,
    y = NULL,
    caption = "Black line: Data"
  )

