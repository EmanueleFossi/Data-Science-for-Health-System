pacchetti <- c("lhs", "ggplot2", "RColorBrewer","gridExtra","patchwork")
installati <- pacchetti %in% installed.packages()[, "Package"]

if (any(!installati)) {
  install.packages(pacchetti[!installati], dependencies = TRUE)
}

library(lhs)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Numero di pazienti
n_pazienti <- 10000

# Parametri medi per Meropenem

# Clearence di Creatinina serve per verificare quanto bene funzionano i reni della mia popolazione
# Litri/ora (Velocità di eliminazione) (I miei pazienti in media puliscono 14.6 L/h di sangue)
CL_pop_media <- 14.6  


# Volume di distribuione della popolazione, come il mio farmaco si distribuisce all' interno
# del corpo in base al BMI di una persona, se è più o meno grassa 
# Litri (Spazio in cui si distribuisce il farmaco)(il farmaco si distribuisce in media su 10.8L di sangue) 
V_pop_medio  <- 10.8


# Variabilità Inter-individuale (espressa come Coefficiente di Variazione - CV)
#Ci dice quanro sono distanti i pazienti l' uno dall' altro rispetto alla media,
# se avessi avuto 0% avevo tutti cloni uguali
cv_pk <- 0.30  # 30% di variabilità tra un paziente e l'altro

# Funzione per convertire Media e CV in parametri della distribuzione Lognormale
# La lognormale non usa media e varianza standard ma usa la meanlog ossia la media dei logaritmi dei dati
# e la sdlog la deviazione standard dei logaritmi dei dati
get_log_params <- function(media, cv) {
  sd_log <- sqrt(log(cv^2 + 1))
  mean_log <- log(media) - 0.5 * sd_log^2
  return(list(meanlog = mean_log, sdlog = sd_log))
}

# Calcoliamo i parametri logaritmici
params_cl <- get_log_params(CL_pop_media, cv_pk)
params_v  <- get_log_params(V_pop_medio, cv_pk)

# Generazione dei pazienti virtuali
set.seed(42)
pazienti <- data.frame(
  id = 1:n_pazienti,
  CL = rlnorm(n_pazienti, params_cl$meanlog, params_cl$sdlog),
  V  = rlnorm(n_pazienti, params_v$meanlog, params_v$sdlog)
)

p1 <- ggplot(pazienti, aes(x = CL)) +
  geom_density(fill = "#238b45", alpha = 0.5, color = "#238b45") +
  theme_minimal() +
  labs(title = "Distribuzione Log-normale della Clearance (N=10.000)",
       subtitle = "Rappresentazione dell'Incertezza Biologica",
       x = "Clearance (L/h)", 
       y = "Densità di Probabilità")


ggsave("Distribuzione_CL_Lognormale.png", plot = p1, width = 8, height = 5, dpi = 300)

p1 <- ggplot(pazienti, aes(x = V)) +
  geom_density(fill = "#238b45", alpha = 0.5, color = "#238b45") +
  theme_minimal() +
  labs(title = "Distribuzione Log-normale della Volume (N=10.000)",
       subtitle = "Rappresentazione dell'Incertezza Biologica",
       x = "Volume di distribuzione (L)", 
       y = "Densità di Probabilità")


ggsave("Distribuzione_V_Lognormale.png", plot = p1, width = 8, height = 5, dpi = 300)

print(head(pazienti))

#print(range(pazienti$CL))
#print(range(pazienti$V))

# --- Provo il mio modello matematico per vedere se funziona ---

# Definizione della funzione base
get_conc_base <- function(cl, v, d, ti, tc) {
  k <- cl / v
  r <- d / ti
  conc <- ifelse(tc <= ti,
                 (r / cl) * (1 - exp(-k * tc)), 
                 ((r / cl) * (1 - exp(-k * ti))) * exp(-k * (tc - ti)))
  return(conc)
}

get_conc <- Vectorize(get_conc_base, vectorize.args = c("cl", "v"))

#test
dose_test <- 2000
t_inf_test <- 3
p1_cl <- pazienti$CL[1]
p1_v  <- pazienti$V[1]

c0   <- get_conc(p1_cl, p1_v, dose_test, t_inf_test, 0)   # Inizio
c1.5 <- get_conc(p1_cl, p1_v, dose_test, t_inf_test, 1.5) # Metà
c3   <- get_conc(p1_cl, p1_v, dose_test, t_inf_test, 3)   # Picco
c8   <- get_conc(p1_cl, p1_v, dose_test, t_inf_test, 8)   # Valle

cat("--- TEST MODELLO PK (Paziente 1) ---\n")
cat("Concentrazione a 0h: ", round(c0, 2), "mg/L\n")
cat("Concentrazione a 3h (Fine Infusione): ", round(c1.5, 2), "mg/L\n")
cat("Concentrazione a 3h (Fine Infusione): ", round(c3, 2), "mg/L\n")
cat("Concentrazione a 8h (Fine Intervallo): ", round(c8, 2), "mg/L\n")

# SIMULAZIONE MONTECARLO 

# (6 x 6 = 36)
array_dosi  <- c(500, 1000, 1500, 2000, 2500, 3000)
array_tempi <- c(0.5, 1, 2, 3, 4, 5)

mic_target   <- 8                     # mg/L (Soglia di resistenza del mio virus + alta + resistente)
target_pd    <- 0.40                  # 40% del tempo (Meropenem, concentrazione sopra la MIC per almeno il 40% del tempo tra una dose e un altra)
tempo_critico <- 8 * target_pd        # 3.2 ore il paziente  ha ancora almeno 8mg/L del farmaco la terapia è un usccesso (8h x 0.4)

griglia_risultati <- expand.grid(Dose = array_dosi, T_infusione = array_tempi)
griglia_risultati$PTA <- NA

for(i in 1:nrow(griglia_risultati)) {
  
  d_scen  <- griglia_risultati$Dose[i]
  ti_scen <- griglia_risultati$T_infusione[i]
  
  conc_popolazione <- get_conc(pazienti$CL, pazienti$V, d_scen, ti_scen, tempo_critico)
  
  successi <- conc_popolazione >= mic_target
  
  griglia_risultati$PTA[i] <- mean(successi) * 100 # moltiplico per cento cosi ho una percentuale
}

scenari_ok <- subset(griglia_risultati, PTA >= 90)

if (nrow(scenari_ok) > 0) {
  
  tabella_report <- scenari_ok
  tabella_report$PTA <- round(as.numeric(tabella_report$PTA), 2)
  tabella_report <- tabella_report[order(tabella_report$Dose, tabella_report$T_inf), ]
  
  colnames(tabella_report) <- c("Dose_mg", "Durata_Infusione_h", "Successo_PTA_Percento")
  write.csv(tabella_report, "Protocolli_Ottimali_Meropenem_MIC8.csv", row.names = FALSE)
  
  cat("\n FILE GENERATO \n")
  print(tabella_report)
  
} else {
  cat("\n FILE NON GENERATO \n")

}

# Heatmap
p <- ggplot(griglia_risultati, aes(x = as.factor(T_infusione), y = as.factor(Dose), fill = PTA)) +
  geom_tile(color = "white", size = 0.5) +  # Aggiunge una griglia bianca tra i quadrati
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), # Rosso -> Giallo -> Verde
    limits = c(0, 100),
    name = "Successo %"
  ) +
  geom_text(aes(label = paste0(round(PTA, 0), "%")), color = "black", fontface = "bold") +
  labs(
    title = "Matrice di Successo Monte Carlo (PTA %)",
    subtitle = paste("Target: 40% del tempo sopra MIC =", mic_target, "mg/L | Popolazione: 10.000 pazienti"),
    x = "Durata Infusione (ore)", 
    y = "Dose (mg)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )
print(p)

ggsave("Risultato_MonteCarlo_Meropenem.png", plot = p, width = 10, height = 7, dpi = 300)

cat("\n--- SIMULAZIONE COMPLETATA ---\n")
cat("Il grafico è stato visualizzato e salvato come 'Risultato_MonteCarlo_Meropenem.png'\n")

# CONFRONTO MONTE CARLO STANDARD VS LHS

set.seed(42)
lhs_quantili <- randomLHS(100, 2) 

cl_lhs <- qlnorm(lhs_quantili[,1], params_cl$meanlog, params_cl$sdlog)
v_lhs  <- qlnorm(lhs_quantili[,2], params_v$meanlog, params_v$sdlog)

test_dose <- 2000
test_ti   <- 3
mic_test  <- 8
t_critico <- 3.2

conc_std <- get_conc(pazienti$CL, pazienti$V, test_dose, test_ti, t_critico)
pta_std  <- mean(conc_std >= mic_test) * 100

conc_lhs <- get_conc(cl_lhs, v_lhs, test_dose, test_ti, t_critico)
pta_lhs  <- mean(conc_lhs >= mic_test) * 100

cat("   CONFRONTO EFFICIENZA METODI MONTE CARLO    \n")
cat("Metodo Standard (10.000 pazienti): ", round(pta_std, 2), "%\n")
cat("Metodo LHS      (100 pazienti):     ", round(pta_lhs, 2), "%\n")
cat("Differenza:                        ", round(abs(pta_std - pta_lhs), 4), "%\n")


# Violin Plot delle Concentrazioni
tempi_test <- c(1, 3, 5) 
df_violin <- data.frame()

for(t in tempi_test) {
  conc_t <- mapply(get_conc, pazienti$CL, pazienti$V, 
                   MoreArgs = list(d = 2000, ti = 3, tc = t))
  df_violin <- rbind(df_violin, data.frame(Tempo = as.factor(t), Conc = conc_t))
}

p_violin <- ggplot(df_violin, aes(x = Tempo, y = Conc, fill = Tempo)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = mic_target, linetype = "dashed", color = "red") +
  labs(title = " Dose 2000mg @ 3h",
       x = "Tempo dalla somministrazione (ore)", y = "Concentrazione (mg/L)") +
  theme_light()
ggsave("ViolinPlot_Concentrazioni.png", plot = p_violin, width = 8, height = 6)

# Curva PK con Fascia di Incertezza (95% CI)
tempi_seq <- seq(0, 8, by = 0.2)
stats_pk <- sapply(tempi_seq, function(t) {
  c_vals <- mapply(get_conc, pazienti$CL, pazienti$V, MoreArgs = list(d = 2000, ti = 3, tc = t))
  return(c(Media = mean(c_vals), 
           L_Inf = as.numeric(quantile(c_vals, 0.025)), 
           L_Sup = as.numeric(quantile(c_vals, 0.975))))
})

df_pk <- as.data.frame(t(stats_pk))
df_pk$Tempo <- tempi_seq

p_curve <- ggplot(df_pk, aes(x = Tempo)) +
  geom_ribbon(aes(ymin = L_Inf, ymax = L_Sup), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = Media), color = "blue", linewidth = 1) +
  geom_hline(yintercept = mic_target, linetype = "dashed", color = "red") +
  labs(title = "Profilo Farmacocinetico di Popolazione",
       subtitle = "Area ombreggiata: variabilità del 95% dei pazienti (Metodo Monte Carlo)",
       x = "Tempo (ore)", y = "Concentrazione (mg/L)") +
  theme_minimal()
ggsave("CurvaPK_Incertezza.png", plot = p_curve, width = 10, height = 6)

print(p_violin)
print(p_curve)