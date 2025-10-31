# MonitoreoContinuo2

MonitoreoContinuo2 es un script de MATLAB pensado para la captura, procesamiento, visualización e interpretación de señales 5G NR recibidas usando un SDR compatible (como USRP). El objetivo principal es detectar y analizar la ráfaga SSB más fuerte y visualizar la evolución temporal de los resource grids.

## Funcionalidad principal

- Configuración automática de parámetros del SDR, banda de operación, frecuencias y OFDM acorde a la banda seleccionada.
- Captura secuencial de la señal 5G en intervalos regulares.
- Procesamiento y demodulación OFDM para cada ráfaga capturada.
- Detección y localización automática de la ráfaga (SSB) más fuerte en cada ventana temporal.
- Cálculo y almacenamiento de métricas de calidad: potencia media, SNR y CellID detectado.
- Visualización interactiva final (con slider temporal) de los resource grids, resaltando la ubicación de la SSB más potente.

## Estructura del script

- **Sección 1:** Configuración de SDR, parámetros de banda y numerología 5G NR.
- **Sección 2:** Definición de duración y parámetros de las capturas sucesivas.
- **Sección 3:** Inicialización de variables para guardar formas de onda, grids y métricas asociadas.
- **Sección 4:** Bucle principal de captura y procesamiento por ráfaga, con almacenamiento de resultados.
- **Sección 5:** Visualizador interactivo que permite explorar la evolución temporal de la señal y la SSB más fuerte.

## Archivos y funciones auxiliares principales

- `findSSB`: Detección y análisis avanzado de la ráfaga SSB (Strongest Synchronization Signal Block).
- `visualizeResourceGridsOverlay`: Visualización dinámica y comparativa de los resource grids.
- `estimateSNR`: Cálculo sencillo de SNR para cada ráfaga.

## Requisitos

- MATLAB R2020a o superior (recomendado).
- Toolboxes: Communications Toolbox, soporte SDR.
- SDR compatible y configurado.
- Scripts auxiliares del ecosistema 5G NR de MathWorks.

## Uso

1. Ejecutar el script `MonitoreoContinuo2.m` en MATLAB.
2. Ajustar parámetros de SDR, banda y frecuencias si es necesario.
3. Explorar resultados y métricas en el workspace y desde el visualizador interactivo.


