# =============================================================
# Recrutamento de cupins em resposta a rompimentos do ninho
# Previsão 1: proximidade ao ninho -> tempo de recrutamento
# Previsão 2: volume do ninho -> nº de cupins e tempo de recrutamento
# =============================================================

# ---- 0. Pacotes ---------------------------------------------
pacotes <- c("glmmTMB", "DHARMa", "performance", "emmeans",
             "ggplot2", "dplyr", "tidyr", "bbmle")
novos <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(novos)) install.packages(novos)
invisible(lapply(pacotes, library, character.only = TRUE))

set.seed(123)

# ---- 1. Leitura e preparo -----------------------------------
dados <- read.csv2("p_final_cupins.csv",
                   fileEncoding = "UTF-8-BOM",
                   stringsAsFactors = FALSE)

# read.csv2 já assume ";" e vírgula decimal; aqui o decimal é ponto:
dados$volume <- as.numeric(as.character(dados$volume))
for (v in c("t_sold_s", "t_ope_s", "t_recrut", "n_cupim", "A_cm", "B_cm", "C_cm")) {
  dados[[v]] <- as.numeric(as.character(dados[[v]]))
}

dados <- dados %>%
  mutate(
    ninho    = factor(ninho),
    proximi  = factor(proximi, levels = c("p10", "p50", "p95")),
    # versão numérica (% de proximidade) para testar tendência linear.
    # ATENÇÃO à direção: p95 = rompimento MAIS PRÓXIMO do ninho,
    # p10 = MAIS DISTANTE. A escala já cresce com a proximidade.
    prox_num = as.numeric(sub("p", "", as.character(proximi))),
    # houve resposta? (t_recrut = 0 e n_cupim = 0 significam ausência de recrutamento)
    resposta = ifelse(t_recrut > 0 | n_cupim > 0, 1, 0),
    # tempo só faz sentido quando houve recrutamento
    t_rec    = ifelse(resposta == 1, t_recrut, NA),
    # volume é muito assimétrico (10 a 406 cm3) -> log + centragem
    log_vol  = as.numeric(scale(log(volume), scale = FALSE))
  )

str(dados)
table(dados$proximi, dados$resposta)
# ninhos sem nenhuma resposta (excluídos das análises de tempo):
dados %>% group_by(ninho) %>% summarise(n_resp = sum(resposta)) %>% filter(n_resp == 0)

# ---- 2. Exploração da distribuição das respostas -------------
par(mfrow = c(2, 2))
hist(dados$t_rec, main = "t_recrut (s)", xlab = "segundos")
hist(log(dados$t_rec), main = "log(t_recrut)", xlab = "log s")
hist(dados$n_cupim, main = "n_cupim", xlab = "nº de cupins")
plot(dados$volume, dados$n_cupim, xlab = "volume (cm3)", ylab = "n_cupim")
par(mfrow = c(1, 1))

# média x variância de n_cupim (indício de sobredispersão se var >> média)
c(media = mean(dados$n_cupim), variancia = var(dados$n_cupim))

# Observação de desenho: cada ninho foi medido nas 3 proximidades
# -> medidas repetidas -> (1 | ninho) como efeito aleatório.
# proximi varia DENTRO do ninho; volume varia ENTRE ninhos (n = nº de ninhos).

# =============================================================
# 3. PREVISÃO 1 - proximidade -> tempo de recrutamento
# =============================================================
# t_recrut é contínuo, positivo e assimétrico à direita:
# candidatos = Gamma(log) e lognormal. Poisson NÃO se aplica (não é contagem).

dt <- droplevels(subset(dados, !is.na(t_rec)))

m1_gamma <- glmmTMB(t_rec ~ proximi + (1 | ninho),
                    family = Gamma(link = "log"), data = dt)
m1_lnorm <- glmmTMB(t_rec ~ proximi + (1 | ninho),
                    family = lognormal(link = "log"), data = dt)
m1_null  <- glmmTMB(t_rec ~ 1 + (1 | ninho),
                    family = Gamma(link = "log"), data = dt)

AICtab(m1_gamma, m1_lnorm, base = TRUE)   # escolha da família pelo AIC

# --- diagnóstico dos resíduos (simulados, DHARMa) ---
res1 <- simulateResiduals(m1_gamma, n = 1000)
plot(res1)                       # QQ (KS, dispersão, outliers) + resíduos x preditos
testDispersion(res1)
testUniformity(res1)
plotResiduals(res1, form = dt$proximi)
# Se o painel do Gamma mostrar desvio e o lognormal não, use m1_lnorm daqui em diante.

# --- inferência ---
summary(m1_gamma)
anova(m1_null, m1_gamma)                     # teste do efeito de proximidade (LRT)
emmeans(m1_gamma, pairwise ~ proximi, type = "response")  # contrastes entre p10/p50/p95
r2(m1_gamma)

# versão com proximidade contínua (testa tendência monotônica, 1 gl)
m1_lin <- glmmTMB(t_rec ~ prox_num + (1 | ninho),
                  family = Gamma(link = "log"), data = dt)
summary(m1_lin)
# coeficiente NEGATIVO de prox_num = tempo diminui conforme aumenta a proximidade
# (= apoio à Previsão 1)

# =============================================================
# 4. PREVISÃO 2a - volume -> número de cupins patrulhando
# =============================================================
# n_cupim é contagem -> Poisson; conferir sobredispersão -> binomial negativa;
# há muitos zeros estruturais (ninhos que não responderam) -> testar hurdle/ZI.

m2_pois <- glmmTMB(n_cupim ~ log_vol + proximi + (1 | ninho),
                   family = poisson, data = dados)
check_overdispersion(m2_pois)

m2_nb1  <- glmmTMB(n_cupim ~ log_vol + proximi + (1 | ninho),
                   family = nbinom1, data = dados)
m2_nb2  <- glmmTMB(n_cupim ~ log_vol + proximi + (1 | ninho),
                   family = nbinom2, data = dados)
m2_zinb <- glmmTMB(n_cupim ~ log_vol + proximi + (1 | ninho),
                   ziformula = ~1, family = nbinom2, data = dados)

AICtab(m2_pois, m2_nb1, m2_nb2, m2_zinb, base = TRUE)

melhor2 <- m2_nb2          # troque pelo modelo com menor AIC
res2 <- simulateResiduals(melhor2, n = 1000)
plot(res2)
testDispersion(res2)
testZeroInflation(res2)
plotResiduals(res2, form = dados$log_vol)

summary(melhor2)
# coeficiente POSITIVO de log_vol = mais cupins em ninhos maiores (apoio à Previsão 2)

m2_semvol <- update(melhor2, . ~ . - log_vol)
anova(m2_semvol, melhor2)     # LRT para o efeito do volume
r2(melhor2)

# =============================================================
# 5. PREVISÃO 2b - volume -> tempo de recrutamento
# =============================================================
m3 <- glmmTMB(t_rec ~ log_vol + proximi + (1 | ninho),
              family = Gamma(link = "log"), data = dt)
res3 <- simulateResiduals(m3, n = 1000); plot(res3)
summary(m3)
anova(update(m3, . ~ . - log_vol), m3)
# coeficiente NEGATIVO de log_vol = ninhos maiores recrutam mais rápido

# --- extra: o volume afeta a PROBABILIDADE de haver resposta? ---
m4 <- glmmTMB(resposta ~ log_vol + proximi + (1 | ninho),
              family = binomial, data = dados)
summary(m4)
# (com poucos ninhos sem resposta este modelo pode ter separação completa;
#  se não convergir, interprete apenas descritivamente)

# =============================================================
# 6. Gráficos
# =============================================================
p1 <- ggplot(dt, aes(proximi, t_rec)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_line(aes(group = ninho), alpha = .3) +
  geom_point(aes(color = ninho), size = 2, show.legend = FALSE) +
  labs(x = "Proximidade do rompimento ao ninho",
       y = "Tempo até o recrutamento (s)") +
  theme_classic()

p2 <- ggplot(dados, aes(volume, n_cupim)) +
  geom_point(aes(shape = proximi), size = 2) +
  geom_smooth(method = "glm", method.args = list(family = "quasipoisson"),
              formula = y ~ log(x), color = "black") +
  labs(x = expression("Volume do ninho (cm"^3*")"),
       y = "Nº de cupins patrulhando") +
  theme_classic()

p3 <- ggplot(dt, aes(volume, t_rec)) +
  geom_point(aes(shape = proximi), size = 2) +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "lm", color = "black") +
  labs(x = expression("Volume do ninho (cm"^3*", log)"),
       y = "Tempo até o recrutamento (s, log)") +
  theme_classic()

print(p1); print(p2); print(p3)
# ggsave("fig1_proximidade.png", p1, width = 5, height = 4, dpi = 300)

# =============================================================
# 7. Tabela-resumo dos modelos finais
# =============================================================
lapply(list(Previsao1 = m1_gamma, Previsao2a = melhor2, Previsao2b = m3),
       function(m) round(summary(m)$coefficients$cond, 4))
