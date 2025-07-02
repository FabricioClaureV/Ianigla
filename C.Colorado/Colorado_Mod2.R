# 📦 Librerías necesarias
library(readr)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(lubridate)
library(HBV.IANIGLA)

message("✅ Librerías cargadas.")

# 📁 Rutas
# Asegúrate de que esta sea la ruta base correcta en tu sistema
ruta_base <- "C:\\Users\\fabri\\Desktop\\C.Colorado"
ruta_entradas <- file.path(ruta_base, "Entradas")
ruta_salidas <- file.path(ruta_base, "Salidas")
if (!dir.exists(ruta_salidas)) { dir.create(ruta_salidas, recursive = TRUE) }

# 📆 Fechas de simulación
fecha_ini <- as.Date("1979-01-01")
fecha_fin <- as.Date("2019-12-31")

## 1. Carga de Datos
message("➡️  1. Cargando datos de entrada...")

# 1.1 Topografía
topo <- read_csv(file.path(ruta_entradas, "topo_colorado_superficies.csv"), show_col_types = FALSE) %>%
  select(mean_elev, rel_area, surface)

# 1.2 Forzantes climáticos (de la estación base)
forzantes <- read_csv(file.path(ruta_entradas, "forcings_colorado.csv"), show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= fecha_ini & Date <= fecha_fin)

# 1.3 Caudales observados
qobs <- read_csv(file.path(ruta_entradas, "caudal_observado_colorado.csv"), show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= fecha_ini & Date <= fecha_fin)

# ####################################################################### #
# ##    SECCIÓN 1.4 MODIFICADA PARA UNA LECTURA DE PARÁMETROS ROBUSTA    ## #
# ####################################################################### #
message("➡️  1.4 Cargando y procesando parámetros desde 'par_colorado.txt'...")

ruta_parametros_txt <- file.path(ruta_entradas, "par_colorado.txt")

# Leer las líneas del archivo
param_lines <- tryCatch({
  readLines(ruta_parametros_txt)
}, error = function(e) {
  stop(paste("Error al leer el archivo de parámetros:", e$message,
             "\nAsegúrate de que la ruta sea correcta y el archivo exista en:", ruta_parametros_txt))
})

# Limpiar las líneas de comentarios y espacios en blanco
param_lines <- gsub("#.*", "", param_lines) # Eliminar comentarios
param_lines <- trimws(param_lines)          # Eliminar espacios al inicio/final
param_lines <- param_lines[param_lines != ""] # Eliminar líneas vacías

# Opciones para el manejo de parámetros duplicados
keep_first_duplicate <- TRUE  # TRUE: se conserva la primera aparición

# Crear listas para los nombres y valores de los parámetros
param_names <- character(0)
param_values <- numeric(0)

# Procesar cada línea para extraer el nombre y el valor de forma segura
# Este bucle ahora verifica que el valor sea un número válido antes de agregarlo
for (line in param_lines) {
  parts <- strsplit(line, "\\s+")[[1]] # Dividir por uno o más espacios

  if (length(parts) >= 2) {
    param_name <- parts[1]
    value_str  <- parts[2]

    # Intentar convertir a número y verificar si es un número válido (no NA)
    numeric_value <- suppressWarnings(as.numeric(value_str))

    if (!is.na(numeric_value)) {
      if (param_name %in% param_names) {
        warning(sprintf("Par\u00e1metro repetido: %s", param_name))
        if (keep_first_duplicate) {
          next
        } else {
          stop(sprintf("Nombres de par\u00e1metro duplicados encontrados: %s", param_name))
        }
      }

      param_names  <- c(param_names, param_name)
      param_values <- c(param_values, numeric_value)
    }
  }
}

# Crear el vector nombrado final
params_vector <- setNames(param_values, param_names)

# Verificar si se leyeron parámetros
if (length(params_vector) == 0) {
  stop("No se pudieron leer parámetros del archivo. Verifica que el archivo no esté vacío y tenga el formato correcto (NOMBRE VALOR).")
}

get_param <- function(param_name, source_vector, module_name) {
  if (!param_name %in% names(source_vector)) {
    stop(paste("Parámetro crítico faltante:", param_name, "para el módulo de", module_name,
               ". Por favor, añádelo a 'par_colorado.txt'."))
  }
  source_vector[param_name]
}

message("✅ Datos cargados.")

# --------------------------------------------------------------------------

## 2. Preparación de Datos de Entrada
message("➡️  2. Preparando matrices de entrada para el modelo...")

elev_estacion <- mean(topo$mean_elev)
grad_temp <- -0.6 / 100
grad_prec <- 0.05 / 100

mat_temp <- matrix(NA, nrow = nrow(forzantes), ncol = nrow(topo))
mat_prec <- matrix(NA, nrow = nrow(forzantes), ncol = nrow(topo))
mat_pet <- matrix(NA, nrow = nrow(forzantes), ncol = nrow(topo))

for (i in 1:nrow(topo)) {
  elev_banda <- topo$mean_elev[i]
  delta_elev <- elev_banda - elev_estacion
  mat_temp[, i] <- forzantes$T + (delta_elev * grad_temp)
  mat_prec[, i] <- forzantes$P * (1 + (delta_elev * grad_prec))
  mat_pet[, i] <- forzantes$PET
}

mat_prec[mat_prec < 0] <- 0
message("✅ Matrices de entrada creadas.")

# --------------------------------------------------------------------------

## 3. Definición de Parámetros por Módulo
message("➡️  3. Agrupando y validando parámetros por módulo...")

param_snow <- unlist(list(
  SFCF          = get_param("SFCF", params_vector, "Nieve/Hielo (SFCF)"),
  TR            = get_param("TT_1", params_vector, "Nieve/Hielo (TR)"),
  TT            = get_param("TT_1", params_vector, "Nieve/Hielo (TT)"),
  FM            = get_param("CFMAX_1", params_vector, "Nieve/Hielo (FM)"),
  CFR           = get_param("CFR", params_vector, "Nieve/Hielo (CFR)"),
  TT_glacier    = get_param("TT_2", params_vector, "Nieve/Hielo (TT_glacier)"),
  CFMAX_glacier = get_param("CFMAX_2", params_vector, "Nieve/Hielo (CFMAX_glacier)"),
  TT_soil       = get_param("TT_3", params_vector, "Nieve/Hielo (TT_soil)"),
  CFMAX_soil    = get_param("CFMAX_3", params_vector, "Nieve/Hielo (CFMAX_soil)")
))

param_soil <- unlist(list(
  FC   = get_param("FC", params_vector, "Suelo (FC)"),
  LP   = get_param("LP", params_vector, "Suelo (LP)"),
  BETA = get_param("BETA", params_vector, "Suelo (BETA)"),
  UZL  = get_param("UZL", params_vector, "Suelo (UZL)"),
  PERC = get_param("PERC", params_vector, "Suelo (PERC)")
))

param_route <- unlist(list(
  K0   = get_param("K0", params_vector, "Ruteo (K0)"),
  K1   = get_param("K1", params_vector, "Ruteo (K1)"),
  K2   = get_param("K2", params_vector, "Ruteo (K2)"),
  UZL  = get_param("UZL", params_vector, "Ruteo (UZL)"),
  PERC = get_param("PERC", params_vector, "Ruteo (PERC)")
))

param_tf <- unlist(list(
  MAXBAS = get_param("MAXBAS", params_vector, "Función de Transferencia (MAXBAS)")
))

message("✅ Parámetros agrupados y validados.")

# --------------------------------------------------------------------------

## 4. Definición de la Función del Modelo
message("➡️  4. Definiendo la función 'hydrological_hbv'...")

hydrological_hbv <- function(basin, tair, precip, pet,
                             param_snow, param_soil, param_route, param_tf,
                             init_snow = 20, init_soil = 0, init_routing = c(0, 0, 0)) {
  n_it <- nrow(basin)
  snow_module  <- list()
  soil_module  <- list()

  for(i in 1:n_it){
    snow_module[[i]] <- SnowGlacier_HBV(model = 1, inputData = cbind(tair[, i], precip[, i]),
                                        initCond =  c(init_snow, 2), param = param_snow)

    soil_module[[i]] <- Soil_HBV(model = 1, inputData = cbind(snow_module[[i]][, 5], pet[, i]),
                                 initCond = c(init_soil, basin$rel_area[i]), param = param_soil)
  }

  soil_disch <- lapply(X = 1:n_it, FUN = function(x) { soil_module[[x]][, 1] })
  soil_disch <- Reduce(f = `+`, x = soil_disch)

  route_module <- Routing_HBV(model = 1, lake = F, inputData = as.matrix(soil_disch),
                              initCond = init_routing, param = param_route)

  tf_module <- round(UH(model = 1, Qg = route_module[, 1], param = param_tf), 4)

  return(tf_module)
}
message("✅ Función definida.")

# --------------------------------------------------------------------------

## 5. Ejecución del Modelo
message("➡️  5. Ejecutando el modelo hidrológico...")

qsim_diario <- hydrological_hbv(
  basin = topo,
  tair = mat_temp,
  precip = mat_prec,
  pet = mat_pet,
  param_snow = param_snow,
  param_soil = param_soil,
  param_route = param_route,
  param_tf = param_tf
)

message("✅ Simulación completada.")

# --------------------------------------------------------------------------

## 6. Procesar Resultados y Evaluar
message("➡️  6. Procesando resultados y evaluando desempeño...")

df_sim <- data.frame(Date = forzantes$Date, Q_sim_m3s = qsim_diario)

qsim_mensual <- df_sim %>%
  mutate(year = year(Date), month = month(Date)) %>%
  group_by(year, month) %>%
  summarise(Q_sim_m3s = mean(Q_sim_m3s, na.rm = TRUE), .groups = "drop") %>%
  mutate(Date = make_date(year, month, 15))

qobs_mensual <- qobs %>%
  mutate(year = year(Date), month = month(Date)) %>%
  group_by(year, month) %>%
  summarise(Q_obs_m3s = mean(Q_obs_m3s, na.rm = TRUE), .groups = "drop") %>%
  mutate(Date = make_date(year, month, 15))

comparacion <- inner_join(qobs_mensual, qsim_mensual, by = "Date")

nse_val <- NA
kge_val <- NA

try({
  nse_val <- NSE(comparacion$Q_sim_m3s, comparacion$Q_obs_m3s)
  kge_val <- KGE(comparacion$Q_sim_m3s, comparacion$Q_obs_m3s)
  metricas <- data.frame(NSE = nse_val, KGE = kge_val)
  print("Métricas de desempeño:")
  print(metricas)
  write_csv(metricas, file.path(ruta_salidas, "metricas_colorado.csv"))
}, silent = TRUE)


g <- ggplot(comparacion, aes(x = Date)) +
  geom_line(aes(y = Q_obs_m3s, color = "Observado"), linewidth = 1) +
  geom_line(aes(y = Q_sim_m3s, color = "Simulado"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Caudal", values = c("Observado" = "blue", "Simulado" = "red")) +
  labs(title = "Caudal Mensual Observado vs. Simulado - Río Colorado",
       subtitle = paste("NSE:", round(nse_val, 2), " | KGE:", round(kge_val, 2)),
       y = "Caudal (m³/s)",
       x = "Fecha") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(g)
ggsave(file.path(ruta_salidas, "Qsim_vs_Qobs_colorado.png"), plot = g, width = 10, height = 6)

message("✅ Proceso finalizado. Revisa la carpeta de Salidas.")
