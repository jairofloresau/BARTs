# ============================================================================
# TESIS: Pronóstico e Identificación de Determinantes de los Retornos
#        Extremos del Tipo de Cambio en Perú: Un Enfoque BART
#
# SCRIPT 01: Construcción de la base de datos diaria
# ============================================================================
#
# Qué hace este script, en orden:
#   1. Carga paquetes y define rutas / fecha de inicio de la muestra
#   2. Define funciones auxiliares para leer cada "formato" de archivo
#      (BCRP, FRED, Bloomberg, NY Fed, CBOE, Yahoo)
#   3. Importa cada fuente de datos usando esas funciones
#   4. Usa el tipo de cambio (variable dependiente) como "ancla" de fechas
#      y une todas las demás series a esa columna de fechas
#   5. Homogeneiza frecuencias (forward-fill para NFCI que es semanal)
#   6. Calcula el tipo de cambio promedio y sus log-retornos
#   7. Exporta la base final en .csv y .rds
#
# ============================================================================


# ============================================================================
# 1. SETUP
# ============================================================================

# --- Paquetes ---
# readxl   : leer archivos .xlsx
# lubridate: manejo de fechas
# dplyr/tidyr: manipulación de datos
# quantmod : descargar DXY de Yahoo Finance
paquetes <- c("readxl", "writexl", "lubridate", "dplyr", "tidyr", "purrr", "stringr", "quantmod")
instalar_faltantes <- paquetes[!paquetes %in% installed.packages()[, "Package"]]
if (length(instalar_faltantes) > 0) install.packages(instalar_faltantes)
invisible(lapply(paquetes, library, character.only = TRUE))

# Aseguramos un locale UTF-8 para que R lea sin problemas nombres de archivo
# con tildes (ej. "Tasa de interés", "Embig perú"). Si tu sistema no tiene
# exactamente este locale, R usará el que tengas por defecto sin fallar.
intentos_locale <- c("en_US.UTF-8", "es_ES.UTF-8", "C.UTF-8", "Spanish")
for (loc in intentos_locale) {
  ok <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error = function(e) NA))
  if (!is.na(ok) && ok != "") break
}

# --- Ruta donde están todos tus archivos ---
# AJUSTA ESTO si mueves los archivos de carpeta
ruta_datos <- "D:/BARTs/data/datos/BASE DIARIO"

# --- Fecha de inicio de la muestra (definida en la conversación) ---
fecha_inicio <- as.Date("2000-01-01")
fecha_fin    <- Sys.Date()  # hasta hoy; ajusta si quieres una fecha fin fija


# ============================================================================
# 2. FUNCIONES AUXILIARES
# ============================================================================
# La idea: cada fuente de datos viene en un "formato de casa" distinto
# (BCRP, FRED, Bloomberg, etc.), así que escribimos UNA función por formato,
# y luego la reutilizamos para cada archivo de ese tipo. Esto evita repetir
# código y hace que si el formato cambia, solo arreglamos un lugar.

## ---- 2.1 Lector de archivos BCRP (Metadatos + Diarias) ----
# Las fechas del BCRP vienen como texto "02Ene97" (día + mes abreviado en
# español + año de 2 dígitos). R no reconoce "Ene" por defecto, así que
# armamos nuestro propio diccionario de meses.
meses_es <- c(
  "Ene" = "01", "Feb" = "02", "Mar" = "03", "Abr" = "04",
  "May" = "05", "Jun" = "06", "Jul" = "07", "Ago" = "08",
  "Set" = "09", "Sep" = "09", "Oct" = "10", "Nov" = "11", "Dic" = "12"
)

parsear_fecha_bcrp <- function(x) {
  # x es un vector tipo "02Ene97", "15Dic23", etc.
  dia <- str_sub(x, 1, 2)
  mes_abr <- str_sub(x, 3, 5)
  anio2 <- str_sub(x, 6, 7)
  mes <- meses_es[mes_abr]
  anio4 <- ifelse(as.integer(anio2) <= 30,
                  paste0("20", anio2),  # 00-30 -> 2000-2030
                  paste0("19", anio2))  # 31-99 -> 1931-1999
  as.Date(paste(anio4, mes, dia, sep = "-"))
}

leer_bcrp <- function(path, nombre_valor) {
  # Lee la hoja "Diarias" saltando las 2 filas de encabezado (código + nombre serie)
  df <- read_excel(path, sheet = "Diarias", col_names = FALSE, skip = 2)
  colnames(df) <- c("fecha_txt", "valor")
  df <- df %>%
    mutate(
      fecha = parsear_fecha_bcrp(fecha_txt),
      # "n.d." (no disponible) y otros textos no numéricos -> NA
      valor = suppressWarnings(as.numeric(valor))
    ) %>%
    filter(!is.na(fecha)) %>%
    select(fecha, !!nombre_valor := valor) %>%
    arrange(fecha)
  df
}

## ---- 2.2 Lector de archivos FRED (README + Daily/Weekly) ----
leer_fred <- function(path, hoja, nombre_valor) {
  df <- read_excel(path, sheet = hoja)
  colnames(df) <- c("fecha", "valor")
  df %>%
    mutate(fecha = as.Date(fecha), valor = as.numeric(valor)) %>%
    filter(!is.na(fecha)) %>%
    select(fecha, !!nombre_valor := valor) %>%
    arrange(fecha)
}

## ---- 2.3 Lector genérico de "bloques" Bloomberg ----
# Los archivos de Bloomberg tienen filas de metadata arriba (ticker, "Fecha
# inicial", "#NAME?", etc.) antes de que empiecen los datos reales. En vez de
# adivinar en qué fila exacta empieza cada bloque, esta función escanea la
# columna de fechas y se queda solo con las filas donde SÍ hay una fecha
# válida. Así es robusta sin importar cuántas filas de metadata haya.
# Las fechas en estos bloques de Bloomberg a veces vienen como fecha real y
# a veces como numero serial de Excel guardado como texto (ej. "35797").
# Ademas, R (a diferencia de Excel/Python) descarta columnas totalmente
# vacias al leer, lo que puede desfasar los indices de columna esperados -
# por eso conviene revisar la hoja cruda (col_names = FALSE) antes de fijar
# col_fecha/col_valor para un archivo nuevo.
parsear_fecha_flexible <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXct")) return(as.Date(x))
  x_chr <- as.character(x)
  out <- as.Date(rep(NA_character_, length(x_chr)))
  es_iso <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", x_chr) & !is.na(x_chr)
  es_num <- grepl("^[0-9]+(\\.[0-9]+)?$", x_chr) & !is.na(x_chr) & !es_iso
  out[es_iso] <- as.Date(substr(x_chr[es_iso], 1, 10))
  out[es_num] <- as.Date(as.numeric(x_chr[es_num]), origin = "1899-12-30")
  out
}

leer_bloque_bloomberg <- function(path, hoja, col_fecha, col_valor, nombre_valor) {
  df <- read_excel(path, sheet = hoja, col_names = FALSE)
  fechas <- parsear_fecha_flexible(df[[col_fecha]])
  valores <- suppressWarnings(as.numeric(df[[col_valor]]))
  
  tibble(fecha = fechas, valor = valores) %>%
    filter(!is.na(fecha)) %>%
    select(fecha, !!nombre_valor := valor) %>%
    arrange(fecha)
}

## ---- 2.4 Lector de fechas seriales de Excel (High Yield - Bloomberg VALORES) ----
leer_high_yield <- function(path) {
  df <- read_excel(path, sheet = "VALORES", col_names = FALSE, skip = 6)
  colnames(df) <- c("fecha_serial", "yield_to_worst")
  df %>%
    mutate(
      fecha = as.Date(as.numeric(fecha_serial), origin = "1899-12-30"),
      yield_to_worst = as.numeric(yield_to_worst)
    ) %>%
    filter(!is.na(fecha)) %>%
    select(fecha, yield_to_worst) %>%
    arrange(fecha)
}


# ============================================================================
# 3. IMPORTAR CADA FUENTE
# ============================================================================

## ---- 3.1 Tipo de cambio interbancario (variable dependiente) ----
tc_compra <- leer_bcrp(file.path(ruta_datos, "TC-Interbancario compra (diario).xlsx"), "tc_compra")
tc_venta  <- leer_bcrp(file.path(ruta_datos, "TC-Interbancario venta (diario).xlsx"),  "tc_venta")

tc <- full_join(tc_compra, tc_venta, by = "fecha") %>%
  mutate(tc = (tc_compra + tc_venta) / 2) %>%   # promedio compra-venta (mid-rate)
  select(fecha, tc) %>%
  arrange(fecha)

## ---- 3.2 Determinantes domésticos (Perú) ----
rin        <- leer_bcrp(file.path(ruta_datos, "RIN netas (diarias).xlsx"), "rin")
tasa_peru  <- leer_bcrp(file.path(ruta_datos, "Tasa de interés (diarias).xlsx"), "tasa_peru")
embig_peru <- leer_bcrp(file.path(ruta_datos, "Embig perú.xlsx"), "embig_peru")

# Operaciones cambiarias BCRP: "n.d." se convierte en NA en leer_bcrp().
# Se trata como dato genuinamente faltante (no como "no hubo intervención").
operaciones_bcrp <- leer_bcrp(file.path(ruta_datos, "Operaciones cambiarias BCRP (diarias).xlsx"), "interv_bcrp")

## ---- 3.3 Determinantes globales / condiciones financieras ----
tasa_usa <- leer_fred(file.path(ruta_datos, "Interest Rate USA (diarias).xlsx"), "Daily", "tasa_usa")
baa10y   <- leer_fred(file.path(ruta_datos, "BAA10Y - spread crediticio corporativo de EE.UU (diarias).xlsx"), "Daily", "baa10y")

# NFCI es SEMANAL (viernes) a pesar del nombre del archivo. Se homogeniza
# a diario más abajo con forward-fill, igual que el resto de variables
# de baja frecuencia.
nfci <- leer_fred(file.path(ruta_datos, "NFCI - Chicago (diario).xlsx"), "Weekly, Ending Friday", "nfci")

high_yield <- leer_high_yield(file.path(ruta_datos, "high yield.xlsx"))

## ---- 3.4 Term Premia (ACM) - plazos 1, 5 y 10 años ----
term_premia <- read_excel(file.path(ruta_datos, "Term Premia.xls"), sheet = "ACM Daily") %>%
  transmute(
    fecha = dmy(DATE),   # formato "14-Jun-1961"
    tp_1y  = ACMTP01,
    tp_5y  = ACMTP05,
    tp_10y = ACMTP10
  ) %>%
  filter(!is.na(fecha)) %>%
  arrange(fecha)

## ---- 3.5 Metales: Cobre y Zinc LME, Oro spot ----
metales_lme <- read_excel(file.path(ruta_datos, "Precios de metales.xlsx"), sheet = "Hoja1",
                          col_names = FALSE, skip = 9) %>%
  transmute(
    fecha = as.Date(...1),
    oro   = as.numeric(...2),
    cobre = as.numeric(...3),
    zinc  = as.numeric(...4)
  ) %>%
  filter(!is.na(fecha)) %>%
  arrange(fecha)

## ---- 3.6 Plata SPOT y WTI (de TESIS EXTRA DATA) ----
ruta_extra <- file.path(ruta_datos, "TESIS EXTRA DATA - new (1).xlsx")

plata <- leer_bloque_bloomberg(ruta_extra, "Precios de metales históricos",
                               col_fecha = 1, col_valor = 2, nombre_valor = "plata")
wti <- leer_bloque_bloomberg(ruta_extra, "Precios de metales históricos",
                             col_fecha = 6, col_valor = 7, nombre_valor = "wti")

## ---- 3.7 Bolsa: SPX, BVL Perú, MXEF ----
spx <- leer_bloque_bloomberg(ruta_extra, "Bolsa", col_fecha = 1, col_valor = 2, nombre_valor = "spx")
mxef <- leer_bloque_bloomberg(ruta_extra, "Bolsa", col_fecha = 6, col_valor = 7, nombre_valor = "mxef")

# NOTA - BVL Perú EXCLUIDA:
# La columna 3 de la hoja "Bolsa" tiene DOS tickers de Bloomberg mezclados en
# las celdas de encabezado: "SPBLPADV Index" (fila 4) y "S&P/BVL Peru Select
# Index TR" (fila 5). Los valores que trae (0-29, promedio ~9) corresponden
# al ticker SPBLPADV, que es un indicador de amplitud de mercado
# (advance/decline: cuántas acciones subieron vs bajaron ese día), NO al
# nivel del índice bursátil. Es un error de la descarga original en
# Bloomberg. Hasta que se vuelva a descargar con el ticker correcto del
# índice, esta variable queda excluida del pipeline.
# bvl <- leer_bloque_bloomberg(ruta_extra, "Bolsa", col_fecha = 1, col_valor = 3, nombre_valor = "bvl_peru")

## ---- 3.8 Otras monedas frente al dólar (JPY, CHF, EUR, GBP, CLP, COP, MXN, BRL, CNY) ----
# Todas comparten una sola columna de fecha (col 1) en la hoja "otros TC"
nombres_monedas <- c("jpy", "chf", "eur", "gbp", "clp", "cop", "mxn", "brl", "cny")
cols_monedas <- 2:10  # columnas 2 a 10 en la hoja, en el mismo orden que nombres_monedas

otras_monedas <- map2(cols_monedas, nombres_monedas, function(col, nombre) {
  leer_bloque_bloomberg(ruta_extra, "otros TC", col_fecha = 1, col_valor = col, nombre_valor = nombre)
}) %>%
  reduce(full_join, by = "fecha")

## ---- 3.9 VIX (se descarga automáticamente desde CBOE, ya que el Excel venía dañado) ----
vix <- read.csv("https://cdn.cboe.com/api/global/us_indices/daily_prices/VIX_History.csv") %>%
  transmute(
    fecha = mdy(DATE),   # formato MM/DD/YYYY
    vix = as.numeric(CLOSE)
  ) %>%
  filter(!is.na(fecha)) %>%
  arrange(fecha)

## ---- 3.10 DXY (se descarga automáticamente desde Yahoo Finance) ----
getSymbols("DX-Y.NYB", src = "yahoo", from = as.character(fecha_inicio), auto.assign = TRUE)
dxy <- data.frame(fecha = index(`DX-Y.NYB`), dxy = as.numeric(`DX-Y.NYB`[, "DX-Y.NYB.Close"])) %>%
  filter(!is.na(dxy)) %>%
  arrange(fecha)


# ============================================================================
# 4. UNIR TODO USANDO EL TIPO DE CAMBIO COMO ANCLA DE FECHAS
# ============================================================================
# El TC diario del BCRP es nuestra variable dependiente, así que sus fechas
# (días hábiles del mercado cambiario peruano) definen el calendario de la
# base final. Todo lo demás se une (left_join) a ese calendario.

base <- tc %>%
  filter(fecha >= fecha_inicio, fecha <= fecha_fin) %>%
  left_join(rin,              by = "fecha") %>%
  left_join(tasa_peru,        by = "fecha") %>%
  left_join(embig_peru,       by = "fecha") %>%
  left_join(operaciones_bcrp, by = "fecha") %>%
  left_join(tasa_usa,         by = "fecha") %>%
  left_join(baa10y,           by = "fecha") %>%
  left_join(nfci,             by = "fecha") %>%
  left_join(high_yield,       by = "fecha") %>%
  left_join(term_premia,      by = "fecha") %>%
  left_join(metales_lme,      by = "fecha") %>%
  left_join(plata,            by = "fecha") %>%
  left_join(wti,              by = "fecha") %>%
  left_join(spx,              by = "fecha") %>%
  left_join(mxef,             by = "fecha") %>%
  left_join(otras_monedas,    by = "fecha") %>%
  left_join(vix,              by = "fecha") %>%
  left_join(dxy,              by = "fecha") %>%
  arrange(fecha)


# ============================================================================
# 5. HOMOGENEIZAR FRECUENCIAS (forward-fill)
# ============================================================================
#  El NFCI es semanal; todo lo demás es diario pero puede tener huecos
# puntuales por feriados que no coinciden entre países (ej. feriado en EE.UU.
# pero no en Perú). Rellenamos esos huecos con el último valor observado
# (Last Observation Carried Forward), el estándar en la literatura para
# este tipo de series.
#
# EXCEPCIÓN: interv_bcrp es una variable de FLUJO (cuánto intervino el BCRP
# ESE día específico), no de nivel/stock. Rellenar sus NA con el valor del
# día anterior implicaría asumir que la intervención se repitió, lo cual no
# es un supuesto razonable. Por eso queda explícitamente fuera del
# forward-fill y conserva sus NA reales.
columnas_a_rellenar <- setdiff(names(base), c("fecha", "interv_bcrp"))

base <- base %>%
  fill(all_of(columnas_a_rellenar), .direction = "down")


# ============================================================================
# 6. VARIABLE DEPENDIENTE: LOG-RETORNOS DEL TIPO DE CAMBIO
# ============================================================================
base <- base %>%
  mutate(retorno_tc = c(NA, diff(log(tc))))

# NOTA: la identificación de "retornos extremos" (percentiles, EVT/POT, etc.)
# todavía no la hemos definido en la conversación - queda pendiente como
# siguiente paso, una vez que confirmes el método a usar.


# ============================================================================
# 7. EXPORTAR BASE FINAL
# ============================================================================
write.csv(base, file.path(ruta_datos, "base_final_diaria.csv"), row.names = FALSE)
saveRDS(base, file.path(ruta_datos, "base_final_diaria.rds"))
write_xlsx(base, file.path(ruta_datos, "base_final_diaria.xlsx"))

cat("Base final construida:", nrow(base), "filas x", ncol(base), "columnas\n")
cat("Rango de fechas:", as.character(min(base$fecha)), "a", as.character(max(base$fecha)), "\n")
cat("Guardado en:\n")
cat(" -", file.path(ruta_datos, "base_final_diaria.csv"), "\n")
cat(" -", file.path(ruta_datos, "base_final_diaria.rds"), "\n")
cat(" -", file.path(ruta_datos, "base_final_diaria.xlsx"), "\n")