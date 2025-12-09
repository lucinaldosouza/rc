install.packages("PNADcIBGE")

library("PNADcIBGE")
library("tidyverse")
library(survey)

# Para família
variaveis_selecionadas <- c("UPA", "V1008", "V1014", "V2003", "VD4019")

# Para indivíduo
variaveis_selecionadas <- c("UPA", "V1008", "V1014", "V2003", "VD4019", "VD4047", "VD5002",
                            "V5001A2", "V2007", "V2008", "V20081", "V20082")
#"VDI4047", "VDI5002"

# Obtenção dos dados para os 5 trimestres consecutivos
PNADc_2017_1 <- get_pnadc(year = 2017, quarter = 1, 
                          vars = variaveis_selecionadas,
                          deflator = TRUE, design = TRUE)

# PNADc_2017_2 <- get_pnadc(year = 2017, quarter = 2, 
#                           vars = variaveis_selecionadas,
#                           deflator = TRUE, design = TRUE)
# 
# PNADc_2017_3 <- get_pnadc(year = 2017, quarter = 3, 
#                           vars = variaveis_selecionadas,
#                           deflator = TRUE, design = TRUE)
# 
# PNADc_2017_4 <- get_pnadc(year = 2017, quarter = 4, 
#                           vars = variaveis_selecionadas,
#                           deflator = TRUE, design = TRUE)

PNADc_2018_1 <- get_pnadc(year = 2018, quarter = 1, 
                          vars = variaveis_selecionadas,
                          deflator = TRUE, design = TRUE)


# Criar a chave de domicílio para os dois trimestres
# PNADc_2017_1 <- update(PNADc_2017_1, 
#                             domicilio_id = paste(UPA, V1008, V1014, sep = "_"))
# PNADc_2017_2 <- update(PNADc_2017_2, 
#                             domicilio_id = paste(UPA, V1008, V1014, sep = "_"))
# PNADc_2017_3 <- update(PNADc_2017_3, 
#                               domicilio_id = paste(UPA, V1008, V1014, sep = "_"))
# PNADc_2017_4 <- update(PNADc_2017_4, 
#                               domicilio_id = paste(UPA, V1008, V1014, sep = "_"))
# PNADc_2018_1 <- update(PNADc_2018_1, 
#                               domicilio_id = paste(UPA, V1008, V1014, sep = "_"))

# Atualizar os objetos de desenho amostral diretamente, garantindo que as chaves sejam reconhecidas
PNADc_2017_1 <- update(PNADc_2017_1, 
                       domicilio_id = paste(UPA, V1008, V1014, sep = "_"),
                       individuo_id = paste(V2007, V2008, V20081, V20082, sep = "_"))

###-----------------------------
PNADc_2017 <- get_pnadc(year = 2017, interview = 5, 
                          vars = variaveis_selecionadas,
                          deflator = TRUE, design = FALSE)

# svyquantile(~VD4019, design = PNADc_2017, c(0.010,0.025,.05,0.075,.10,.15,.25,.5,.75),na.rm = TRUE)
# svyquantile(~VD4019, design = PNADc_2017, seq(0, 1, by = 0.1),na.rm = TRUE)

#media_transferencias = mean(VD4047, na.rm = TRUE)
#svymean(~VD4047,design = PNADc_2017, na.rm = TRUE)

#pnad_design2 <- pnad_design %>%
#  filter(!is.na(VD4019) & VD4019 > 0 & !is.na(VD4047) & VD4047 > 0) %>%  
#  mutate(soma_beneficios = svytotal(~VD4047, design = pnad_design2, na.rm = TRUE)) 

##

pnadc_design <- svydesign(ids = ~UPA, 
                          strata = ~Estrato, 
                          weights = ~V1032,  # Ajuste o nome do peso amostral se necessário
                          data = PNADc_2017,
                          nest = TRUE)

dadosPNADc <- pnadc_design(data_pnadc=PNADc_2017)

svytotal(x=~VD4019, design=dadosPNADc, na.rm=TRUE)
pnadc_design %>%
  summarise(Total = survey_total(x = VD4019, vartype = "ci", na.rm=T))

# Definição do desenho amostral
pnadc_design <- PNADc_2017$design

# Criar uma nova variável de faixa de renda na base de dados original
PNADc_2017$data <- PNADc_2017$data %>%
  mutate(faixa_renda = cut(VD4019,
                           breaks = quantile(VD4019, probs = seq(0, 1, 0.1), na.rm = TRUE),
                           labels = paste0("Decil ", 1:10),
                           include.lowest = TRUE
  ))

###-----------------------------

PNADc_2018_1 <- update(PNADc_2018_1, 
                       domicilio_id = paste(UPA, V1008, V1014, sep = "_"),
                       individuo_id = paste(V2007, V2008, V20081, V20082, sep = "_"))


# Extrair as variáveis dos dois trimestres para acessar as chaves de domicílio
# chaves_2017_T1 <- PNADc_2017_1$variables$domicilio_id
# chaves_2017_T2 <- PNADc_2017_2$variables$domicilio_id
# chaves_2017_T3 <- PNADc_2017_3$variables$domicilio_id
# chaves_2017_T4 <- PNADc_2017_4$variables$domicilio_id
# chaves_2018_T1 <- PNADc_2018_1$variables$domicilio_id

# Extrair as variáveis dos dois trimestres para acessar as chaves de domicílio
df_2017_1 <- PNADc_2017_1$variables
df_2018_1 <- PNADc_2018_1$variables


# Identificar os domicílios que estão em ambos os trimestres
# domicilios_comuns_5_trimestres <- Reduce(intersect, list(chaves_2017_T1, chaves_2017_T2, 
#                                                          chaves_2017_T3, chaves_2017_T4, 
#                                                          chaves_2018_T1))

# Identificar os indivíduos que estão nas mesmas famílias
# Primeiro, identificamos os domicílios em comum
domicilios_comuns <- intersect(df_2017_1$domicilio_id, df_2018_1$domicilio_id)


# Filtrar os dados para manter apenas os domicílios comuns
df_2017_1_comuns <- df_2017_1 %>%
  filter(domicilio_id %in% domicilios_comuns)

df_2018_1_comuns <- df_2018_1 %>%
  filter(domicilio_id %in% domicilios_comuns)


# Agora, identificar os indivíduos que estão presentes nos dois períodos
individuos_comuns <- intersect(df_2017_1_comuns$individuo_id, df_2018_1_comuns$individuo_id)


# Filtrar as bases de desenho amostral originais para manter apenas os indivíduos comuns
PNADc_2017_1_comuns <- subset(PNADc_2017_1, individuo_id %in% individuos_comuns)
PNADc_2018_1_comuns <- subset(PNADc_2018_1, individuo_id %in% individuos_comuns)

# Exemplo: calcular a média do rendimento efetivo (VD4019) em 2017 e 2018
media_rendimento_2017_1 <- svymean(~VD4019, design = PNADc_2017_1_comuns, na.rm = TRUE)
media_rendimento_2018_1 <- svymean(~VD4019, design = PNADc_2018_1_comuns, na.rm = TRUE)

# Exibir os resultados
media_rendimento_2017_1
media_rendimento_2018_1


# Define a function to calculate income percentiles
calculate_percentiles_unique <- function(data, variable, percentiles = seq(0, 1, 0.01)) {
  quantiles <- quantile(data[[variable]], percentiles, na.rm = TRUE)
  # Ensure breaks are unique by adding a small jitter if necessary
  quantiles <- jitter(quantiles, factor = 1e-7)
  return(quantiles)
}

# Calculate percentiles for VD4019 in 2017_1 with unique breaks
percentiles_2017_1 <- calculate_percentiles_unique(df_2017_1_comuns, "VD4019")

# Calculate percentiles for VD4019 in 2018_1 with unique breaks
percentiles_2018_1 <- calculate_percentiles_unique(df_2018_1_comuns, "VD4019")


# Classify families into percentile categories based on 2017_1
df_2017_1_comuns <- df_2017_1_comuns %>%
  mutate(percentile_2017_1 = cut(VD4019, breaks = percentiles_2017_1, labels = FALSE, include.lowest = TRUE))

# Classify families into percentile categories based on 2018_1
df_2018_1_comuns <- df_2018_1_comuns %>%
  mutate(percentile_2018_1 = cut(VD4019, breaks = percentiles_2018_1, labels = FALSE, include.lowest = TRUE))

# Merge the data from 2017_1 and 2018_1 based on family identifier (e.g., domicilio_id)
merged_data <- inner_join(df_2017_1_comuns, df_2018_1_comuns, by = "domicilio_id", suffix = c("_2017", "_2018"))

# Now, you can see the transitions between percentiles
# Example: How many families moved from one percentile to another?
table(merged_data$percentile_2017_1, merged_data$percentile_2018_1)



# 
##
df_combined <- bind_rows(
  df_2017_1_comuns %>% mutate(year = "2017"),
  df_2018_1_comuns %>% mutate(year = "2018")
)

######

# Compute deciles using svyquantile()
deciles_2017_1 <- svyquantile(~VD4019, design = PNADc_2017_1_comuns, 
                              quantiles = seq(0, 1, 0.1), na.rm = TRUE)

# Extract the quantile values for VD4019
decile_values_2017_1 <- deciles_2017_1$VD4019

# Step 2: Ensure breaks are unique
# Check for duplicated decile values
duplicated_deciles_2017_1 <- duplicated(decile_values_2017_1)
decile_values_2017_1 <- as.numeric(decile_values_2017_1)

duplicated_deciles_2017_1 <- duplicated(decile_values_2017_1)

# Step 3: Assign decile groups using adjusted breaks
df_2017_1_comuns <- df_2017_1_comuns %>%
  mutate(percentile_2017_1 = cut(VD4019, 
                                 breaks = decile_values_2017_1, 
                                 labels = FALSE, 
                                 include.lowest = TRUE))

######


# Divide income into percentiles (10 groups, i.e., deciles)
df_2017_1_comuns <- df_2017_1_comuns %>%
  mutate(percentile_2017_1 = ntile(VD4019, 10))  # ntile divides into deciles

df_2018_1_comuns <- df_2018_1_comuns %>%
  mutate(percentile_2018_1 = ntile(VD4019, 10))

# Merge the datasets based on 'domicilio_id' to track families between years
df_transition <- inner_join(df_2017_1_comuns, df_2018_1_comuns, by = "domicilio_id", suffix = c("_2017", "_2018"))

# Create a transition matrix: count how many moved between percentiles
transition_matrix <- table(df_transition$percentile_2017_1, df_transition$percentile_2018_1)

# Compute transition probabilities: normalize by row sums to get probabilities
transition_probabilities <- prop.table(transition_matrix, 1)  # Normalize by row (2017 percentiles)

# Print transition matrix and transition probabilities
transition_matrix
transition_probabilities


##########################################
##########################################
##########################################
## Até aqui são calculadas as transições #
##########################################
##########################################
##########################################


# Renda média por decil em 2017
mean_income_by_decile_2017 <- df_2017_1_comuns %>%
  filter(!is.na(percentile_2017_1)) %>%  # Exclude rows where the decile is NA
  group_by(percentile_2017_1) %>%
  summarise(mean_income_2017 = mean(VD4019, na.rm = TRUE))

# Renda média por decil em 2018
mean_income_by_decile_2018 <- df_2018_1_comuns %>%
  filter(!is.na(percentile_2018_1)) %>%  # Exclude rows where the decile is NA
  group_by(percentile_2018_1) %>%
  summarise(mean_income_2018 = mean(VD4019, na.rm = TRUE))

# View the mean income for each decile in 2017
mean_income_by_decile_2017

# View the mean income for each decile in 2018
mean_income_by_decile_2018

# Step 3: Extract the mean income of the 1st decile for both 2017 and 2018
first_decile_mean_2017 <- mean_income_by_decile_2017 %>% filter(percentile_2017_1 == 1) %>% pull(mean_income_2017)
first_decile_mean_2018 <- mean_income_by_decile_2018 %>% filter(percentile_2018_1 == 1) %>% pull(mean_income_2018)


# Step 4: Normalize the transition matrix for 2017 by the 1st decile income of 2017
normalized_matrix_2017 <- mean_income_by_decile_2017$mean_income_2017 / first_decile_mean_2017

# Step 5: Normalize the transition matrix for 2018 by the 1st decile income of 2018
normalized_matrix_2018 <- mean_income_by_decile_2018$mean_income_2018 / first_decile_mean_2018

# Usar dados normalizados de 2017 ou 2018? P/ mim, faz mais sentido de 2017, pois
# é o começo do período analisado.

# Não faria mais sentido normalizar pela renda maior? Nesse caso, os decís mais
# baixos seriam interpretados como um choque forte negativo, enquanto os decís
# mais altos como agentes que quase não sofrem choques

# Normalizando pelo último decil

last_decile_mean_2017 <- mean_income_by_decile_2017 %>% filter(percentile_2017_1 == 7) %>% pull(mean_income_2017)
last_decile_mean_2018 <- mean_income_by_decile_2018 %>% filter(percentile_2018_1 == 7) %>% pull(mean_income_2018)

# Step 4: Normalize the income for each decile in 2017 by the 10th decile
normalized_matrix_2017 <- mean_income_by_decile_2017$mean_income_2017 / last_decile_mean_2017

# Step 5: Normalize the income for each decile in 2018 by the 10th decile
normalized_matrix_2018 <- mean_income_by_decile_2018$mean_income_2018 / last_decile_mean_2018

##########################################
## Capital + transferências públicas
# Variáveis adicionando V5002A2 (transferências), V5007A2 e V5008A2 (alugueis e ren-
# dimentos)
variaveis_selecionadas_2 <- c("UPA", "V1008", "V1014", "V2003", "VD4019",  # Variáveis já utilizadas
                            "V2007", "V2008", "V20081", "V20082",       # Identificação dos indivíduos
                            "V5002A2", "V5007A2", "V5008A2")              # Novas variáveis


PNADc_2017_anual <- get_pnadc(year = 2017, interview = 5, 
                              vars = variaveis_selecionadas_2,
                              deflator = TRUE, design = TRUE)

# Criar a chave de indivíduo com base em V2007, V2008, V20081 e V20082
df_2017_anual <- PNADc_2017_anual$variables %>%
  mutate(individuo_id = paste(V2007, V2008, V20081, V20082, sep = "_"))

# Certifique-se de que as variáveis V5002A2, V5007A2, e V5008A2 estão como numéricas
df_2017_anual <- df_2017_anual %>%
  mutate(
    V5002A2 = as.numeric(as.character(V5002A2)),   # Transferências
    V5007A2 = as.numeric(as.character(V5007A2)),   # Aluguel
    V5008A2 = as.numeric(as.character(V5008A2))    # Outros rendimentos
  )


# Criar a renda de capital somando V5007A2 e V5008A2
df_2017_anual <- df_2017_anual %>%
  mutate(rendimento_capital = V5007A2 + V5008A2)

# Agrupar por indivíduo para calcular a média de transferências e de renda de capital
df_2017_anual_agrupado <- df_2017_anual %>%
  group_by(individuo_id) %>%
  summarise(
    media_transferencias = mean(V5002A2, na.rm = TRUE),       # Média de transferências
    media_rendimento_capital = mean(rendimento_capital, na.rm = TRUE),  # Média de renda de capital
    renda = mean(VD4019, na.rm = TRUE)
  )

# Calcular os decis de transferências e renda de capital a nível individual
decis_transferencias <- quantile(df_2017_anual_agrupado$media_transferencias, probs = seq(0, 1, 0.1), na.rm = TRUE)
decis_rendimento_capital <- quantile(df_2017_anual_agrupado$media_rendimento_capital, probs = seq(0, 1, 0.1), na.rm = TRUE)

## Calculado o número de indivíduos em cada decil
# Assign deciles to each individual
df_2017_anual_agrupado <- df_2017_anual_agrupado %>%
  mutate(
    decil_transferencias = cut(media_transferencias,
                               breaks = decis_transferencias,
                               include.lowest = TRUE,
                               labels = 1:10),
    decil_rendimento_capital = cut(media_rendimento_capital,
                                   breaks = decis_rendimento_capital,
                                   include.lowest = TRUE,
                                   labels = 1:10)
  )

# Compute the number of individuals in each decile for transfers
counts_transferencias <- df_2017_anual_agrupado %>%
  filter(!is.na(decil_transferencias)) %>%
  group_by(decil_transferencias) %>%
  summarise(numero_individuos = n())

# Compute the number of individuals in each decile for capital income
counts_rendimento_capital <- df_2017_anual_agrupado %>%
  filter(!is.na(decil_rendimento_capital)) %>%
  # group_by(decil_rendimento_capital) %>%
  summarise(numero_individuos = n())


mean(df_2017_anual_agrupado$media_transferencias, na.rm = TRUE)


## Proporção de indivíduos recebendo transferência
#
# Número total de indivíduos
total_individuals <- nrow(df_2017_anual_agrupado)

# N de indivíduos recebendo transferência
recipients <- df_2017_anual_agrupado %>%
  filter(!is.na(media_transferencias) & media_transferencias > 0)

num_recipients <- nrow(recipients)

proportion_recipients <- num_recipients / total_individuals

# Media e mediana
mean_transfer <- mean(recipients$media_transferencias, na.rm = TRUE)
median_transfer <- median(recipients$media_transferencias, na.rm = TRUE)


# Variancia e desvio padrão 
variance_transfer <- var(recipients$media_transferencias, na.rm = TRUE)
sd_transfer <- sd(recipients$media_transferencias, na.rm = TRUE)


## Montante transferido por Decis de renda
#
# Assuming you have income data, create income deciles
df_2017_anual_agrupado <- df_2017_anual_agrupado %>%
  mutate(
    decil_renda = ntile(renda, 10)  # VD4019 is the income variable
  )

# Calculate average transfer amount by income decile
avg_transfer_by_decile <- df_2017_anual_agrupado %>%
  filter(!is.na(media_transferencias) & media_transferencias > 0) %>%
  group_by(decil_renda) %>%
  summarise(
    mean_transfer = mean(media_transferencias, na.rm = TRUE),
    median_transfer = median(media_transferencias, na.rm = TRUE)
  )

print(avg_transfer_by_decile)

# Calculate probability of receiving transfers by income decile
prob_transfer_by_decile <- df_2017_anual_agrupado %>%
  group_by(decil_renda) %>%
  summarise(
    probability = mean(media_transferencias > 0, na.rm = TRUE)
  )

print(prob_transfer_by_decile)
#####

# Filter out observations with NA in decil_renda and media_transferencias
filtered_data <- df_2017_anual_agrupado %>%
  filter(!is.na(decil_renda)) %>%
  filter(!is.na(media_transferencias) & media_transferencias > 0)

# Total de beneficiarios
total_recipients <- df_2017_anual_agrupado %>%
  filter(!is.na(media_transferencias) & media_transferencias > 0) %>%
  nrow()

# Beneficiários em cada decil
recipients_by_decile <- filtered_data %>%
  group_by(decil_renda) %>%
  summarise(
    num_recipients = n()
  ) %>%
  ungroup() %>%
  mutate(
    percent_of_total_recipients = (num_recipients / total_recipients) * 100
  )

recipients_by_decile

library(dplyr)

df_filtered <- df_2017_anual_agrupado %>%
  filter(!is.na(renda) & !is.na(decil_renda))


# Compute the income intervals for each decile
income_intervals_by_decile <- df_filtered %>%
  group_by(decil_renda) %>%
  summarise(
    min_income = min(renda, na.rm = TRUE),
    max_income = max(renda, na.rm = TRUE),
    mean_income = mean(renda, na.rm = TRUE),
    median_income = median(renda, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(decil_renda)

counts_per_individual <- df_2017_anual %>%
  group_by(individuo_id) %>%
  summarise(
    num_observations = n()
  )

#PNADc_mais_atual <- get_pnadc(year = 2023, interview = 5, 
#                              vars = variaveis_selecionadas_2,
#                              deflator = TRUE, design = TRUE)


