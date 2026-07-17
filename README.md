# Proyecto de Investigacion

Estructura base para organizar un proyecto de investigacion reproducible.

## Estructura

- `data/`: datos (raw, interim, processed, external)
- `notebooks/`: exploracion y experimentos
- `src/`: codigo fuente del proyecto
- `reports/`: figuras y tablas para resultados
- `references/`: articulos y material bibliografico
- `docs/`: documentacion
- `config/`: archivos de configuracion
- `outputs/`: salidas de ejecucion
- `tests/`: pruebas

## Siguientes pasos

1. Agregar tus datos en `data/`.
2. Desarrollar notebooks en `notebooks/`.
3. Migrar logica estable a `src/`.
4. Documentar hallazgos en `reports/` y `docs/`.

## Flujo QBART

El flujo principal esta en `QBART_riesgo_de_cola.R`. Usa por defecto el
kernel C++ oficial de BayesQArt, tuning separado por percentil, predictores
causales de volatilidad, ventana movil y calibracion rolling.

```powershell
# Instalar dependencias, incluido BayesQArt oficial en un commit fijado
Rscript --vanilla scripts/install_qbart_deps.R

# Pruebas unitarias y de interfaz
Rscript --vanilla tests/test_qbart_components.R

# Simulation gate rapido
$env:QBART_SIM_MODE="quick"
Rscript --vanilla tests/qbart_simulation_gate.R

# Corrida de humo sobre los datos reales
$env:QBART_MODE="prueba"
Rscript --vanilla QBART_riesgo_de_cola.R
```

No se debe ejecutar `QBART_MODE=completa` hasta que el simulation gate completo
cumpla cobertura y desempeno relativo frente al QRF.
