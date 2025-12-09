library(PNADcIBGE)
library(tidyverse)
library(dplyr)
library(srvyr) # para uso com dplyr
library(survey)

variaveis_selecionadas <- c("UPA", "Estrato", "V1008", "V1014", "V2003","V1032", 
                            "VD4019", "VD4047", "VD5002", "V5001A2", "V5002A2",
                            "V5003A2","V2007", "V2008", "V20081", "V20082")

# Entrevista 5 de 2017
dadosPNADc <- get_pnadc(year=2017, interview=5, labels=TRUE, deflator=TRUE, design=FALSE, vars=variaveis_selecionadas)

# desenho amostral (automático)
# dadosPNADc <- pnadc_design(data_pnadc=dadosPNADc)
# class(dadosPNADc)

# desenho amostral manual, com a vantagem de uso do dplyr
pnad_design <- dadosPNADc %>%
  as_survey_design(
    ids = UPA,
    strata = Estrato,
    weights = V1032,
    nest = TRUE
  )

###############################################
### Inserindo código de cálculo dos quintis ###
# 1. Soma de transferências por domicílio e renda domiciliar
pnad_domicilios <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%
  group_by(ID_DOMICILIO) %>%
  summarise(
    renda_domicilio = sum(VD4019, na.rm = TRUE),
    soma_rendas_domicilio = sum(c_across(c(V5001A2, V5002A2, V5003A2)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(quintil_renda = ntile(renda_domicilio, 5))

# 2. Unindo os quintis aos dados individuais
pnad_individuos <- pnad_design %>%
  as_tibble() %>%
  left_join(pnad_domicilios, by = "ID_DOMICILIO") %>%
  mutate(
    soma_rendas_ind = rowSums(across(c(V5001A2, V5002A2, V5003A2)), na.rm = TRUE),
    recebeu_transferencia = soma_rendas_ind > 0
  ) %>%
  as_survey_design(ids = UPA, strata = Estrato, weights = V1032, nest = TRUE)


# 3. Proporção de indivíduos com transferência por quintil
resultado_quintis <- pnad_individuos %>%
  group_by(quintil_renda) %>%
  summarise(
    prop_recebeu = survey_mean(recebeu_transferencia, vartype = "ci", na.rm = TRUE)
  )

resultado_quintis
#
# Agora por renda domiciliar per capta
#

# 1. Utilizando o peso (v2003) (primeira pessoa em cada domicílio: menor V2003)
pesos_domicilio <- pnad_design %>%
  as_tibble() %>%
  group_by(ID_DOMICILIO) %>%
  filter(V2003 == min(V2003)) %>%
  ungroup() %>%
  select(ID_DOMICILIO, peso_dom = V1032)

# 2. Agregando por domicílio e calculando renda per capita e transferências
pnad_domicilios <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%
  group_by(ID_DOMICILIO) %>%
  summarise(
    renda_total = sum(VD4019, na.rm = TRUE),
    n_moradores = n(),
    renda_per_capita = renda_total / n_moradores,
    soma_transferencias = sum(c_across(c(V5001A2, V5002A2, V5003A2)), na.rm = TRUE),
    recebeu_transferencia = soma_transferencias > 0,
    .groups = "drop"
  ) %>%
  left_join(pesos_domicilio, by = "ID_DOMICILIO") %>%
  mutate(
    quintil_renda_pc = ntile(renda_per_capita, 5)
  )

# 3. Objeto de survey com peso da primeira pessoa do domicílio
pnad_domicilios_design <- pnad_domicilios %>%
  as_survey_design(ids = ID_DOMICILIO, weights = peso_dom)

# 4. Proporção de domicílios com transferência em cada quintil
resultado_domicilios <- pnad_domicilios_design %>%
  group_by(quintil_renda_pc) %>%
  summarise(
    prop_recebeu = survey_mean(recebeu_transferencia, vartype = "ci", na.rm = TRUE)
  )

resultado_domicilios


### Fim do código de cálculo dos quintis ###
############################################

############################################
# DECIS
# Calculo dos decis da renda do trabalho
pnad_design2 <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%  # Remover casos sem renda do trabalho
  mutate(decil_renda = ntile(VD4019, 10))  # Criando os 10 grupos de renda

# Proporção da renda dos 10% mais ricos
prop_renda_top10 <- pnad_design2 %>%
  summarise(
    renda_total = survey_total(VD4019),
    renda_top10 = survey_total(if_else(decil_renda == 10, VD4019, 0))
  ) %>%
  mutate(
    proporcao_top10 = renda_top10 / renda_total
  )

prop_renda_top10


# QUINTIS
# Pegando os quintis
pnad_design2 <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%  # Remover casos sem renda do trabalho
  mutate(quintil_renda = ntile(VD4019, 5))  # Criando os 10 grupos de renda

# resultado que cheguei mas não bate com os dados -> renda maior recebendo maior transferência
# como solucionar a seguinte


##### Soma das variáveis "V5001A2", "V5002A2", "V5003A2"
pnad_design3 <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%
  mutate(quintil_renda = ntile(VD4019, 5),
         soma_rendas = rowSums(across(c(V5001A2, V5002A2, V5003A2)), na.rm = TRUE))

# Agrupar por decil e calcular a soma ponderada da variável soma_rendas
resultado <- pnad_design3 %>%
  group_by(quintil_renda) %>%
  summarise(total_soma = survey_total(soma_rendas, na.rm = TRUE))


###
# 10 grupos de renda
pnad_design2 <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%
  mutate(decil_renda = ntile(VD4019, 10),
         soma_rendas = rowSums(across(c(V5001A2, V5002A2, V5003A2)), na.rm = TRUE))

resultado2 <- pnad_design2 %>%
  group_by(decil_renda) %>%
  summarise(total_soma = survey_total(soma_rendas, na.rm = TRUE))

library(scales)

resultado2 %>%
  mutate(
    total_soma = comma(total_soma),
    total_soma_se = comma(total_soma_se)
  ) %>%
  print()



# V5001A2: 445385 linhas com 439261 NA, 6124 linhas restantes
# V5002A2: 445385 linhas com 420981 NA, 24404 linhas restantes
# V5003A2: 445385 linhas com 443818 NA, 1567 linhas restantes

# length(dadosPNADc$V5001A2)
# sum(is.na(dadosPNADc$V5001A2)

# resultado
# A tibble: 10 × 3
# decil_renda total_soma total_soma_se
#       <int>       <dbl>         <dbl>
# 1         1 316.503.428     10618176.
# 2         2 175.705.999      6998828.
# 3         3  90.060.633      5725885.
# 4         4  69.130.246      6051775.
# 5         5  56.570.896      5134878.
# 6         6  29.655.130      4429855.
# 7         7  31.898.212      3974095.
# 8         8  14.046.206      2522055.
# 9         9   7.924.772      1496024.
# 10        10   4.105.926      1308500.


# Calculando a mediana de renda de cada decil
resultado_mediana <- pnad_design2 %>%
  group_by(decil_renda) %>%
  summarise(
    mediana_renda = survey_median(VD4019, na.rm = TRUE)
  )

# Com os resultados formatados
resultado_mediana %>%
  mutate(mediana_renda = scales::comma(mediana_renda)) %>%
  print()
# Tabela resultado_mediana 
# decil_renda mediana_renda mediana_renda_se
#       <int> <chr>                    <dbl>
# 1         1 200                     10.2  
# 2         2 500                     12.8  
# 3         3 900                      7.65 
# 4         4 937                      0.765
# 5         5 1,000                    3.83 
# 6         6 1,200                   11.5  
# 7         7 1,500                   12.8  
# 8         8 2,000                   25.5  
# 9         9 3,000                   25.5  
# 10        10 5,500                  153. 


# Agora, calcular a renda mediana da população
renda_mediana_populacao <- pnad_design %>%
  filter(!is.na(VD4019) & VD4019 > 0) %>%  # Filtramos apenas aqueles com renda positiva
  summarise(renda_mediana = survey_median(VD4019, na.rm = TRUE))

# Exibindo a renda mediana formatada
renda_mediana_populacao %>%
  mutate(renda_mediana = scales::comma(renda_mediana )) %>%
  print()
# mediana = 1,254

# Quantidade de pessoas em cada decil de renda
tamanho_decis <- pnad_design2 %>%
  group_by(decil_renda) %>%
  summarise(tamanho = survey_total(1, na.rm = TRUE))

# Exibir o resultado formatado
tamanho_decis %>%
  mutate(tamanho = scales::comma(tamanho)) %>%
  print()
# A tibble: 10 × 3
# decil_renda tamanho    tamanho_se
#       <int> <chr>           <dbl>
# 1         1 6,346,296      91372.
# 2         2 7,505,007      93764.
# 3         3 7,020,243      95262.
# 4         4 8,482,645     116895.
# 5         5 9,309,404     130529.
# 6         6 10,195,910    152646.
# 7         7 9,400,373     124270.
# 8         8 9,767,815     131855.
# 9         9 9,748,375     157181.
# 10       10 10,573,204    313509.

# Calculando T_i

# Df da soma
total_soma_df <- tibble(
  decil_renda = 1:10,
  total_soma = c(316503428, 175705999, 90060633, 69130246, 56570896, 
                 29655130, 31898212, 14046206, 7924772, 4105926)
)

# Df de n_i
tamanho_df <- tibble(
  decil_renda = 1:10,
  tamanho = c(6346296, 7505007, 7020243, 8482645, 9309404, 
              10195910, 9400373, 9767815, 9748375, 10573204)
)

# T_i
Ti_df <- total_soma_df %>%
  inner_join(tamanho_df, by = "decil_renda") %>%
  mutate(T_i = total_soma / tamanho)  # Cálculo de T_i

# Gráfico de T_i
library(ggplot2)

# Cálculo T_i = renda_total_i / tamanho_i onde i = decil
dados_Ti <- tibble(
  decil_renda = 1:10,
  T_i = c(49.872150, 23.411837, 12.828706, 8.149610, 6.076747, 
          2.908570, 3.392637, 1.438360, 0.812998, 0.388302)
)


# Gráfico
ggplot(dados_Ti, aes(x = decil_renda, y = T_i)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +  
  labs(
    title = "Valor Médio de Transferência Recebida por Decil",
    x = "Decil de Renda",
    y = expression(T[i]),
    caption = "Fonte: Cálculo próprio com dados da PNAD Contínua"
  ) +
  theme_minimal() 

# Para ti_df
ggplot(Ti_df, aes(x = decil_renda, y = T_i)) +
  geom_line(color = "blue", size = 1) +  # Linha azul conectando os pontos
  geom_point(color = "blue", size = 2) +  # Pontos vermelhos para destaque
  labs(
    title = "Valor médio de transferência recebida por decil",
    x = "Decil de renda",
    y = expression(T[i]),
    caption = "Fonte: Cálculo próprio com dados da PNAD Contínua"
  ) +
  theme_minimal()

#############################
## Cálculo da função perda ##
#############################

library(stats)

# Função T do modelo
T_modelo <- function(y, xi, m, ybar) {
  return (m * ybar * (2 * exp(-xi * (y / ybar))) / (1 + exp(-xi * (y / ybar))))
}

# Função de Soma dos Erros Quadráticos
sse <- function(params, Y, T_obs, ybar) {
  xi <- params[1]
  m  <- params[2]
  T_pred <- T_modelo(Y, xi, m, ybar)
  return(sum((T_obs - T_pred)^2))
}

# Definir chute inicial para o parâmetro (xi) (\zeta)
x0 <- c(xi = 1.0, m = 0.1)

Y_obs <- c(200,500,900,937,1000,1200,1500,2000,3000, 5500)

# Minimizando a SSE usando Nelder-Mead (Brent para unidimensional)
res <- optim(par = x0, fn = sse, Y = Y_obs, T_obs = Ti_df$T_i, ybar = 1254, method = "Nelder-Mead")

# Extrair os valores estimados de xi
xi_hat <- res$par[1]
m_hat <- res$par[2]

# Exibir os resultados
cat(sprintf("xi = %.4f", xi_hat))




# Verificando graficamente
xi_hat <- 3.66292653
m_hat  <- 0.05467039
ybar   <- 1254

# Ti normalizado sobre a mediana da renda
Ti_df <- Ti_df %>%
  mutate(T_i_normalizado = T_i / ybar)

# Vetor de renda média por decil (X):
Y_obs <- c(200, 500, 900, 937, 1000, 1200, 1500, 2000, 3000, 5500)
# Vetor de transferência média observada (Y):
T_obs <- Ti_df$T_i_normalizado

# Função predita
T_modelo <- function(y, xi, m, ybar) {
  m * ybar * (2 * exp(-xi * (y / ybar))) / (1 + exp(-xi * (y / ybar)))
}

xi_hat <- 3.66292653
m_hat  <- 0.05467039
ybar   <- 1254

# Crie um grid de renda para visualizar a curva
y_grid <- seq(200, 5500, length.out = 10)
T_pred_grid <- T_modelo(y_grid, xi_hat, m_hat, ybar)

# Plotar
plot(
  x = Y_obs, y = T_pred_grid / ybar,
  pch = 16, col = "blue",
  xlab = "Renda (y)",
  ylab = "Transferência T(y)",
  main = "Ajuste da função de transferência"
)
lines(
  x = y_grid,
  y = T_pred_grid / ybar,
  lwd = 2
)


ggplot(Ti_df, aes(x = decil_renda, y = Ti_df$T_i_normalizado)) +
  geom_line(color = "blue", size = 1) +  # Linha azul conectando os pontos
  geom_point(color = "red", size = 3) +  # Pontos vermelhos para destaque
  labs(
    title = "Ajuste do valor médio de transferência recebida por decil",
    x = "Decil de renda",
    y = expression(T[i]),
    caption = "Fonte: Cálculo próprio com dados da PNAD Contínua"
  ) +
  theme_minimal()

# P1: Ti como números maiores que 1

# Ti normalizado sobre a mediana da renda
Ti_df <- Ti_df %>%
  mutate(T_i_normalizado = T_i / ybar)



##########################
# Testes e experimentos  #
##########################

# Função T(y) normalizada
T_modelo_norm <- function(y, xi, m, ybar) {
  (m * ybar * (2 * exp(-xi * (y / ybar))) / (1 + exp(-xi * (y / ybar)))) / ybar
}

# Parâmetros calibrados
xi_hat <- 3.66292653
m_hat  <- 0.05467039
ybar   <- 1254

# Criar um grid mais denso de renda
y_grid_2 <- seq(0, 5 * ybar, length.out = 200)

# Curva 1: Calibrada
T_calibrada <- T_modelo_norm(y_grid_2, xi = xi_hat, m = m_hat, ybar = ybar)

# Curva 2: Lower level (mesmo xi, menor m)
T_lower <- T_modelo_norm(y_grid_2, xi = xi_hat, m = 0.03, ybar = ybar)

# Curva 3: Slower phase-out (mesmo m, menor xi)
T_slower <- T_modelo_norm(y_grid_2, xi = 1.0, m = m_hat, ybar = ybar)

# Plotar tudo no mesmo gráfico
plot(
  y_grid_2 / ybar, T_calibrada, type = "l", lwd = 2, col = "blue",
  xlab = expression(y / bar(y)), ylab = "Transfer over median income",
  ylim = c(0, 0.2),
  main = "Funções de Transferência: Calibração e Contrafactuais"
)
lines(y_grid_2 / ybar, T_lower, col = "darkorange", lwd = 2, lty = 2)
lines(y_grid_2 / ybar, T_slower, col = "goldenrod", lwd = 2, lty = 3)
abline(h = 0, col = "black", lwd = 1, lty = 3)

# Legenda
legend("topright", legend = c("Calibration", "Lower level", "Slower phase-out"),
       col = c("blue", "darkorange", "goldenrod"),
       lwd = 2, lty = c(1, 2, 3), bty = "n")
