#1. Instalación y carga de paquetes
install.packages(c("HBV.IANIGLA", "hydroToolkit", "Evapotranspiration", "hydroGOF"))
library(HBV.IANIGLA)
library(hydroToolkit)
library(Evapotranspiration)
library(hydroGOF)




#2. Obtención y preparación de datos meteorológicos
# 2. Define período, coordenadas y características de la cuenca 5405001 (Río Los Riecillos)
start_date        <- "1980-01-01"
end_date          <- "2020-12-31"
lon               <- -70.32    # longitud en grados decimales
lat               <- -32.75    # latitud en grados decimales

# Datos de la cuenca
area_km2          <- 233.4     # km²
area_m2           <- area_km2 * 1e6  # convierte a m²
prec_media_anual  <- 393       # mm/año (CR2MET)
indice_aridez     <- 1.6

# Cotras de elevación (msnm)
z_max             <- 4844      # cota máxima
z_media           <- 3487      # cota media
z_salida          <- 1922      # cota en el punto de salida

# Mostrar por control
print(sprintf("Período: %s a %s", start_date, end_date))
print(sprintf("Estación (lon,lat) = (%.2f, %.2f)", lon, lat))
print(sprintf("Área = %.1f km²; Precip anual = %d mm; Índice aridez = %.1f",
              area_km2, prec_media_anual, indice_aridez))



