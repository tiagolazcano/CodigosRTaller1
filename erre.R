library(dplyr)
library(lubridate)
library(writexl)
library(fields)

MUE_BIO_TRAMPA_TIERRA_2022 <- MUE_BIO_TRAMPA_TIERRA_2022 %>%
  mutate(
    FECHA_HORA_ZARPE = ymd_hms(FECHA_HORA_ZARPE),  # Convertir a formato fecha-hora
    FECHA = as.Date(FECHA_HORA_ZARPE),
    HORA = format(FECHA_HORA_ZARPE, "%H:%M:%S"),
    ANO = year(FECHA_HORA_ZARPE),
    MES = month(FECHA_HORA_ZARPE),
    DIA = day(FECHA_HORA_ZARPE)
  )

TRAMPA_TIERRA_2022 <- read_excel("D:/Users/Tiago Lazcano/Desktop/Taller 1/Asesorías+/Asesoría 1/TRAMPA_TIERRA_2022.xlsx")

TRAMPA_TIERRA_2022 <- TRAMPA_TIERRA_2022 %>%
  filter(COD_REGION == 10, COD_ESPECIE == 66)

attach(MUE_BIO_TRAMPA_TIERRA_2022)
table(MUE_BIO_TRAMPA_TIERRA_2022$COD_ESPECIE)

table(MUE_BIO_TRAMPA_TIERRA_2022$PUERTO_RECALADA)

table(COD_ESPECIE[COD_ESPECIE == 66 & PUERTO_RECALADA == 960 & SEXO == 2], MES[COD_ESPECIE == 66 & PUERTO_RECALADA == 960 & SEXO ==2])


table(SEXO[PUERTO_RECALADA == 947])

mean(ANCHO_MM[PUERTO_RECALADA == 947 & COD_ESPECIE == 66 & SEXO == 1] )


write_xlsx(TRAMPA_TIERRA_2022, "TRAMPA_TIERRA_2022.xlsx")





TRAMPA_BUCEO_2022 <- read_excel("D:/Users/Tiago Lazcano/Desktop/Taller 1/Asesorías+/Asesoría 1/MUE_BIO_BUCEO_TIERRA_2022.xlsx")

TRAMPA_BUCEO_2022 <- TRAMPA_BUCEO_2022 %>%
  filter(RECURSO == 66)

write_xlsx(TRAMPA_BUCEO_2022, "TRAMPA_BUCEO_MARMOLA_2022.xlsx")

############## BAYES IC #################

library(rjags)
library(readxl)
library(dplyr)

TRAMPAS_MARMOLA_2022 = read_excel("TRAMPAS_MARMOLA_2022.xlsx")
View(TRAMPAS_MARMOLA_2022)

datos = TRAMPAS_MARMOLA_2022 %>%
  filter(SEXO == 1, ARTE == 2, MES == 1)

plot(density(datos$ANCHO_MM))
max(datos$ANCHO_MM)

y = datos$ANCHO_MM
n = length(y)

talla_media_mes = function(arreglo, mes, sexo, puerto, arte){
  
  tabla = table(arreglo$ANCHO_MM[arreglo$SEXO == sexo & arreglo$MES == mes &
                                   arreglo$PUERTO_RECALADA == puerto &
                                   arreglo$ARTE == arte])
  prop = numeric(length(tabla))
  mult = numeric(length(tabla))
  talla <- as.numeric(row.names(tabla))
  
  tabla = as.matrix(tabla)
  suma = sum(tabla)
  
  for (i in 1:length(tabla)) {
    
    prop[i] = tabla[i]/suma
    
  }
  
  
  for (i in 1:length(tabla)) {
    
    mult[i] = talla[i]*prop[i]
    
  }
  
  sum(mult)
}


talla_media_mes(TRAMPAS_MARMOLA_2022, 1, 1, 947, 2)

results = matrix(0, nrow = 12, ncol = 3)

for(i in 1:12) {
  
  # Tabla de tallas por mes
  tab <- table(TRAMPAS_MARMOLA_2022$ANCHO_MM[
    TRAMPAS_MARMOLA_2022$MES == i &
      TRAMPAS_MARMOLA_2022$SEXO == 1 &
      TRAMPAS_MARMOLA_2022$PUERTO_RECALADA == 947 &
      TRAMPAS_MARMOLA_2022$ARTE == 1
  ])
  
  # Si el mes no tiene datos, pasa al siguiente
  if(length(tab) == 0) next
  
  talla  <- as.numeric(names(tab))   # t_j
  counts <- as.numeric(tab)          # n_j
  K <- length(talla)
  
  # Modelo: Dirichlet sobre pesos de tallas
  model_string <- "
  model{
    w[1:K] ~ ddirch(alpha[1:K])
    theta <- inprod(w[1:K], talla[1:K])   # media de talla (estimador)
  }
  "
  
  j <- jags.model(textConnection(model_string),
                  data = list(K=K, talla=talla, alpha=1 + counts),
                  n.chains = 2, quiet=TRUE)
  update(j, 1000)
  
  samp <- coda.samples(j, "theta", n.iter = 3000)
  th <- as.numeric(as.matrix(samp))
  
  # Resumen posterior
  media_post <- mean(th)
  IC95 <- quantile(th, c(0.025, 0.975))
  var_post <- var(th)
  
  # Guardar en la matriz
  results[i, 1] <- var_post
  results[i, 2] <- as.numeric(IC95[1])
  results[i, 3] <- as.numeric(IC95[2])
}
round(results,4)


#########

TRAMPAS_MARMOLA_2022 = read_excel("TRAMPAS_MARMOLA_2022.xlsx")
TRAMPAS_MARMOLA_2022$VIAJE = paste(TRAMPAS_MARMOLA_2022$COD_BARCO,
                                   TRAMPAS_MARMOLA_2022$FECHA_HORA_ARRIBO, sep = "")

modelo_alpha0 = "
model{

  # Prior estable para alpha0
  alpha0 ~ dgamma(2, 1)

  # Dirichlet base para p
  p_raw[1:K] ~ ddirch(rep1[])

  for(k in 1:K){
    p[k] <- max(p_raw[k], 1e-6)
    alpha_k[k] <- alpha0 * p[k]
  }

  # Nivel viajes
  for(i in 1:Nviajes){

    w[i,1:K] ~ ddirch(alpha_k[])

    counts_adj[i,1:K] ~ dmulti(w[i,], N_adj[i])

    theta_i[i] <- inprod(w[i,], talla[])
  }

  for(i in 1:Nviajes){
    weight[i] <- N_adj[i] / totalN_adj
  }

  theta_mes <- inprod(weight[], theta_i[])
}
"

estimacion_mes_alpha0 = function(df, mes, sexo, puerto, arte,
                                 n.chains = 3, n.burnin = 5000, n.iter = 10000,
                                 vars = c("theta_mes","alpha0")){

  df_mes = df %>%
    filter(MES == mes, SEXO == sexo, PUERTO_RECALADA == puerto, ARTE == arte)

  if(nrow(df_mes) == 0) 
    return(list(media = NA, var = NA, IC = c(NA,NA), msg = "sin datos"))

  tab = df_mes %>%
    group_by(VIAJE, talla = ANCHO_MM) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = talla, values_from = n, values_fill = 0)

  tallas = as.numeric(colnames(tab)[-1])
  counts = as.matrix(tab[,-1], row.names = NULL)

  Nviajes = nrow(counts)
  K = length(tallas)

  N_i = rowSums(counts)
  totalN = sum(N_i)

  # Pseudo-conteos
  counts_adj = counts + 1
  N_adj = rowSums(counts_adj)
  totalN_adj = sum(N_adj)

  # ----------------------------------------------------
  # 🔧 BLOQUE data_jags CORREGIDO
  # ----------------------------------------------------
  data_jags = list(
    counts_adj   = counts_adj,
    N_adj        = N_adj,
    totalN_adj   = totalN_adj,
    talla        = tallas,     # <- CORREGIDO
    K            = K,
    Nviajes      = Nviajes,
    rep1         = rep(1, K)
  )
  # ----------------------------------------------------

  j = jags.model(textConnection(modelo_alpha0),
                 data = data_jags,
                 n.chains = n.chains,
                 quiet = TRUE)

  update(j, n.burnin)

  samp = coda.samples(j, variable.names = vars, n.iter = n.iter)

  theta = as.numeric(as.matrix(samp[ , "theta_mes"]))

  out = list(
    media = mean(theta),
    var   = var(theta),
    IC    = quantile(theta, c(0.025, 0.975)),
    mcmc  = samp
  )

  return(out)
}

# Ejemplo:
out = estimacion_mes_alpha0(TRAMPAS_MARMOLA_2022,
                            mes = 1, sexo = 2, puerto = 947, arte = 2)

out$media; out$var; out$IC
