#  EXAFS Analyzer

<p align="center">
<img width="50%" alt="icon" src="https://github.com/user-attachments/assets/7ace0080-c1c7-4968-8e1c-18dd4ff39ec3" />
</p>


[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
--
![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)

**EXAFS Analyzer** es una interfaz gráfica (GUI) avanzada y de código abierto construida en Python (PyQt5) para el procesamiento, análisis y ajuste de datos de espectroscopía de absorción de rayos X (XAS / EXAFS). 
Este software actúa como un entorno visual unificado ("frontend") para la potente librería científica **Larch** y el motor de dispersión cuántica **FEFF8**, simplificando el flujo de trabajo desde la importación de datos crudos hasta la obtención de parámetros estructurales en la primera capa de coordinación.

**ACTUALMENTE PROGRAMADO PARA LEER FICHEROS DE ADVANCED PHOTON SOURCE (APS) (Archivos "setup-full-XXXXX")**

---
> **Install the required dependencies** (run this command in the folder containing the script):
> ```bash
> pip install -r requirements.txt
> ```
>  **Run the application by typing the following in your terminal (inside the script folder):**
> ```bash
> python "EXAFS_Analyzer.py"
> ```
> 
>  **Create a Standalone Executable (.exe) (Run by typing in your terminal inside the script folder)**:
> ```bash
> pyinstaller --onefile --noconsole --icon=icon.ico --exclude-module PyQt6 "EXAFS_Analyzer.py"
> ```

---


##  Características Principales

El programa está diseñado con un flujo de trabajo en 6 pestañas secuenciales y un lienzo de visualización global que se actualiza en tiempo real:

* **1. Explorador Universal de Datos:** Carga múltiple de archivos `.txt`, `.csv` o `.xlsx`.
  * Detección inteligente de columnas y trazado visual instantáneo.
  * Cálculo de medias e interpolación para promediar espectros.
  * Gráficos de diferencia entre muestras.
* **2. Normalización (Pre-edge):** Ajuste interactivo de E0, pre-edge y post-edge para obtener la señal normalizada.
* **3. Extracción del Background (Autobk):** Ajuste de Spline optimizado mediante el parámetro `Rbkg`.
  * Extracción visual de las oscilaciones de EXAFS puras con diferentes pesos (k-weight).
* **4. Transformada de Fourier (Espacio R):** Paso al espacio de distancias aplicando funciones de ventana y suavizado de bordes.
* **5. Motor FEFF Integrado:** Generador automático de archivos `feff.inp` a partir de estructuras cristalográficas (`.cif`).
  * Ajuste avanzado de parámetros termodinámicos, tamaño del clúster y criterios de dispersión.
  * Ejecución nativa de FEFF y extracción automática de las rutas de dispersión (paths) en una tabla interactiva.
* **6. Ajuste Dinámico Multipath:** Ajuste simultáneo en el espacio R y k.
  * **Autodetección inteligente** de la primera capa de coordinación según el rango R seleccionado.
  * Ajuste de parámetros globales y locales para cada path (N, Delta R, Sigma2).
* **Utilidades Globales:** Exportación rápida gráfica a un archivo `.txt` estructurado en columnas, listo para graficar en OriginLab o Excel.

---

## Guía Rápida de Uso
* Datos: En la pestaña 1, carga tus datos del sincrotrón. Selecciona tu columna de Energía en el "Eje X" y la de Absorción en el "Eje Y".
  
<img width="1392" height="879" alt="image" src="https://github.com/user-attachments/assets/e29d34b2-403f-4751-b1a6-f1dd7584281b" />

* Procesamiento: Avanza por las pestañas 2, 3 y 4 ajustando los parámetros para aislar correctamente tu señal EXAFS y obtener una transformada de Fourier limpia.
  
<img width="1399" height="884" alt="image" src="https://github.com/user-attachments/assets/f3311542-a714-4e09-a327-751f929070c0" />
<img width="1387" height="886" alt="image" src="https://github.com/user-attachments/assets/2b726210-6638-48f8-91bb-a9093643e003" />
<img width="1393" height="878" alt="image" src="https://github.com/user-attachments/assets/efe26c89-82b4-4436-9e62-e65ed3e35acc" />

* Paths Teóricos: Ve a la pestaña 5. Carga un archivo .cif de tu molécula teórica, selecciona tu átomo central y pulsa en "Ejecutar FEFF".
<img width="1393" height="882" alt="image" src="https://github.com/user-attachments/assets/f8091313-26b1-4cea-8187-390acfdc68e6" />

* Fiteo: En la pestaña 6, define tu ventana R. Haz clic en  Autodetectar 1ª Capa para cargar los paths de FEFF correspondientes, y haz clic en Ejecutar para ver los resultados de tu ajuste estructural.
  
<img width="1439" height="853" alt="image" src="https://github.com/user-attachments/assets/718f5ea3-491c-4146-9929-42f4fcf6f71e" />

* Exportar: Utiliza el botón verde 💾 Exportar Datos (TXT) en el panel derecho en cualquier momento para guardar tus curvas.

