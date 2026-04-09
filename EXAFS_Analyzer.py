import os
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as style
from matplotlib.figure import Figure

from larch import Group
from larch.xafs import pre_edge, autobk, xftf, feffrunner
from larch.fitting import param
from larch.xafs import (FeffPathGroup, feffit_transform, 
                        feffit_dataset, feffit, feffit_report)

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QGridLayout, QSplitter, QFrame, QScrollArea, QTabWidget, 
    QPushButton, QFileDialog, QLabel, QComboBox, QListWidget, 
    QAbstractItemView, QMessageBox, QCheckBox, QDoubleSpinBox, 
    QSpinBox, QPlainTextEdit, QTableWidget, QTableWidgetItem, 
    QHeaderView, QFormLayout
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class AnalizadorTotal(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EXAFS Analyzer")
        self.setGeometry(100, 100, 1400, 850) 
        
        self.file_paths = []
        self.columnas_disponibles = []
        self.data_cache = {} 
        self.datos_promedio_actual = None 
        self.datos_xas = None 
        self.bloquear_plot = False 

        style.use('ggplot')
        self.initUI()

    def initUI(self):
           
           central_widget = QWidget()
           self.setCentralWidget(central_widget)
           main_layout = QHBoxLayout(central_widget)
   
           # ==========================================
           # PANEL IZQUIERDO 
           # ==========================================
           # Creamos un widget contenedor para el panel izquierdo
           left_container = QWidget()
           left_panel = QVBoxLayout(left_container)
           
           self.tabs = QTabWidget()
           
           # Pestañas (Igual que antes)
           self.tab_explorador = QWidget()
           self.init_explorador_tab()
           self.tabs.addTab(self.tab_explorador, "1. Datos")
   
           self.tab_exafs = QWidget()
           self.init_exafs_tab()
           self.tabs.addTab(self.tab_exafs, "2. Pre-edge")
   
           self.tab_autobk = QWidget()
           self.init_autobk_tab()
           self.tabs.addTab(self.tab_autobk, "3. Background")
           
           self.tab_ft = QWidget()
           self.init_ft_tab()
           self.tabs.addTab(self.tab_ft, "4. Espacio R (FT)")
           
           self.tab_feff = QWidget()
           self.init_feff_tab()
           self.tabs.addTab(self.tab_feff, "5. FEFF (Paths)")
           
           self.tab_fit = QWidget()
           self.init_fit_tab()
           self.tabs.addTab(self.tab_fit, "6. Fit")
           
           self.tabs.currentChanged.connect(self.actualizar_vista_segun_pestana)
   
           left_panel.addWidget(self.tabs)
           # Eliminamos el addStretch() para que el QScrollArea de los paths 
           # dentro de las pestañas pueda ocupar todo el espacio disponible.
   
            # ==========================================
           # PANEL DERECHO 
           # ==========================================
           right_container = QWidget()
           right_panel = QVBoxLayout(right_container)
           
           # --- EL TRUCO ESTÁ AQUÍ: Un layout horizontal ---
           h_tools = QHBoxLayout()
           
           self.fig = Figure()
           self.canvas = FigureCanvas(self.fig)
           
           # 1. Metemos la barra típica de Matplotlib (la que tiene el disquete)
           self.toolbar = NavigationToolbar(self.canvas, self)
           
           # 2. Creamos nuestro súper-botón de exportar a TXT
           self.btn_export = QPushButton("💾 Exportar Datos (TXT)")
           self.btn_export.setStyleSheet("background-color: #27ae60; color: white; font-weight: bold; padding: 6px 15px; border-radius: 4px;")
           self.btn_export.clicked.connect(self.exportar_grafica_txt)
           
           # 3. Los ponemos en fila india: Toolbar -> Espacio vacío -> Botón verde
           h_tools.addWidget(self.toolbar)
           h_tools.addStretch() 
           h_tools.addWidget(self.btn_export)
           
           # 4. Añadimos esta fila superior al panel, y justo debajo la gráfica gigante
           right_panel.addLayout(h_tools)
           right_panel.addWidget(self.canvas)
           
           # --- La Tabla de Paths (Oculta por defecto) ---
           self.tabla_paths = QTableWidget()
           self.tabla_paths.setColumnCount(5)
           self.tabla_paths.setHorizontalHeaderLabels(["Archivo", "Amplitud (%)", "Degen (N)", "Distancia (Å)", "Scattering"])
           self.tabla_paths.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
           self.tabla_paths.setAlternatingRowColors(True)
           self.tabla_paths.setStyleSheet("QTableWidget { background-color: #fdfdfd; font-size: 12px; } QHeaderView::section { background-color: #34495e; color: white; font-weight: bold; padding: 4px; }")
           self.tabla_paths.setMaximumHeight(220)
           self.tabla_paths.setVisible(False)
           
           right_panel.addWidget(self.tabla_paths)
   
           # ==========================================
           # SPLITTER (EL DIVISOR MÓVIL)
           # ==========================================
           splitter = QSplitter(Qt.Horizontal)
           splitter.addWidget(left_container)
           splitter.addWidget(right_container)
           
           # Definimos los tamaños iniciales en píxeles. 
           # Ajusta 500 para la izquierda y 1000 para la derecha según tu pantalla.
           splitter.setSizes([500, 1000]) 
           
           # Añadimos el splitter al layout principal de la ventana
           main_layout.addWidget(splitter)
    # ---------------------------------------------------------
    # INICIALIZACIÓN DE PESTAÑAS
    # ---------------------------------------------------------
    def init_explorador_tab(self):
        layout = QVBoxLayout(self.tab_explorador)
        self.btn_cargar = QPushButton("Cargar Archivos")
        self.btn_cargar.setStyleSheet("font-weight: bold; padding: 10px; background-color: #4CAF50; color: white;")
        self.btn_cargar.clicked.connect(self.cargar_archivos)
        layout.addWidget(self.btn_cargar)
        self.lbl_archivos = QLabel("Archivos en memoria: 0")
        layout.addWidget(self.lbl_archivos)

        layout.addWidget(QLabel("<b>Eje X:</b>")); self.combo_x = QComboBox()
        self.combo_x.currentIndexChanged.connect(self.procesar_y_plotear_explorador)
        layout.addWidget(self.combo_x)

        layout.addWidget(QLabel("<b>Ejes Y (múltiple):</b>")); self.lista_y = QListWidget()
        self.lista_y.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.lista_y.itemSelectionChanged.connect(self.procesar_y_plotear_explorador)
        layout.addWidget(self.lista_y)

        self.check_promedio = QCheckBox("Calcular promedio e interpolar"); self.check_promedio.setChecked(True)
        self.check_promedio.stateChanged.connect(self.procesar_y_plotear_explorador)
        layout.addWidget(self.check_promedio)

        self.check_diferencia = QCheckBox("Mostrar gráfico de diferencia (2 columnas)"); self.check_diferencia.setChecked(False)
        self.check_diferencia.stateChanged.connect(self.procesar_y_plotear_explorador)
        layout.addWidget(self.check_diferencia)

    def init_exafs_tab(self):
        layout = QVBoxLayout(self.tab_exafs)
        form_layout = QFormLayout()

        self.spin_e0 = QDoubleSpinBox(); self.spin_e0.setRange(0, 100000); self.spin_e0.setValue(0)
        self.spin_pre1 = QDoubleSpinBox(); self.spin_pre1.setRange(-500, 0); self.spin_pre1.setValue(-150)
        self.spin_pre2 = QDoubleSpinBox(); self.spin_pre2.setRange(-200, 0); self.spin_pre2.setValue(-30)
        self.spin_norm1 = QDoubleSpinBox(); self.spin_norm1.setRange(0, 1000); self.spin_norm1.setValue(150)
        self.spin_norm2 = QDoubleSpinBox(); self.spin_norm2.setRange(0, 2000); self.spin_norm2.setValue(800)

        for spin in [self.spin_e0, self.spin_pre1, self.spin_pre2, self.spin_norm1, self.spin_norm2]:
            spin.valueChanged.connect(self.procesar_xas_pipeline)

        form_layout.addRow("E0 (eV, 0=Auto):", self.spin_e0)
        form_layout.addRow("Pre-edge start:", self.spin_pre1)
        form_layout.addRow("Pre-edge end:", self.spin_pre2)
        form_layout.addRow("Norm start:", self.spin_norm1)
        form_layout.addRow("Norm end:", self.spin_norm2)
        layout.addLayout(form_layout); layout.addStretch()

    def init_autobk_tab(self):
        """Controles para el Background (Autobk)"""
        layout = QVBoxLayout(self.tab_autobk)
        form_layout = QFormLayout()

        # Rbkg: Define la "rigidez" del spline. Valores bajos = spline flexible.
        self.spin_rbkg = QDoubleSpinBox()
        self.spin_rbkg.setRange(0.1, 5.0); self.spin_rbkg.setValue(1.0); self.spin_rbkg.setSingleStep(0.1)
        self.spin_rbkg.valueChanged.connect(self.procesar_xas_pipeline)

        # Rango de k para el ajuste
        self.spin_kmin = QDoubleSpinBox()
        self.spin_kmin.setRange(0.0, 10.0); self.spin_kmin.setValue(0.0)
        self.spin_kmin.valueChanged.connect(self.procesar_xas_pipeline)

        self.spin_kmax = QDoubleSpinBox()
        self.spin_kmax.setRange(2.0, 30.0); self.spin_kmax.setValue(15.0)
        self.spin_kmax.valueChanged.connect(self.procesar_xas_pipeline)

        # Peso en k (k-weight)
        self.spin_kweight = QDoubleSpinBox()
        self.spin_kweight.setRange(0, 3); self.spin_kweight.setValue(2); self.spin_kweight.setDecimals(0)
        self.spin_kweight.valueChanged.connect(self.procesar_xas_pipeline)

        form_layout.addRow("Rbkg (Å):", self.spin_rbkg)
        form_layout.addRow("k min (Å⁻¹):", self.spin_kmin)
        form_layout.addRow("k max (Å⁻¹):", self.spin_kmax)
        form_layout.addRow("k-weight:", self.spin_kweight)
        layout.addLayout(form_layout)
        
        lbl_info = QLabel("<i>El Rbkg determina qué oscilaciones son 'background' (bajas frecuencias) y cuáles son señal EXAFS.</i>")
        lbl_info.setWordWrap(True)
        layout.addWidget(lbl_info)
        layout.addStretch()
        
    def init_ft_tab(self):
            """Controles para la Transformada de Fourier (Paso a espacio R)"""
            layout = QVBoxLayout(self.tab_ft)
            form_layout = QFormLayout()
    
            # Selector de ventana (Hanning es el estándar en EXAFS)
            self.combo_window = QComboBox()
            self.combo_window.addItems(['hanning', 'kaiser', 'welch', 'parzen', 'sine'])
            self.combo_window.currentTextChanged.connect(self.procesar_xas_pipeline)
    
            # Controles de los límites de la ventana en k
            self.spin_ft_kmin = QDoubleSpinBox()
            self.spin_ft_kmin.setRange(0.0, 10.0); self.spin_ft_kmin.setValue(3.0)
            self.spin_ft_kmin.valueChanged.connect(self.procesar_xas_pipeline)
    
            self.spin_ft_kmax = QDoubleSpinBox()
            self.spin_ft_kmax.setRange(2.0, 30.0); self.spin_ft_kmax.setValue(12.0)
            self.spin_ft_kmax.valueChanged.connect(self.procesar_xas_pipeline)
    
            # dk: El "suavizado" de los bordes de la ventana
            self.spin_ft_dk = QDoubleSpinBox()
            self.spin_ft_dk.setRange(0.0, 5.0); self.spin_ft_dk.setValue(1.0); self.spin_ft_dk.setSingleStep(0.1)
            self.spin_ft_dk.valueChanged.connect(self.procesar_xas_pipeline)
    
            form_layout.addRow("Tipo de Ventana:", self.combo_window)
            form_layout.addRow("k min (Inicio ventana):", self.spin_ft_kmin)
            form_layout.addRow("k max (Fin ventana):", self.spin_ft_kmax)
            form_layout.addRow("dk (Suavizado de bordes):", self.spin_ft_dk)
    
            layout.addLayout(form_layout)
            
            lbl_info = QLabel("<i>Ajusta k max hasta donde tus datos tengan buena relación señal/ruido. El parámetro dk apaga los extremos para evitar 'ripples' en el espacio R.</i>")
            lbl_info.setWordWrap(True)
            layout.addWidget(lbl_info)
            layout.addStretch()
            
    def init_feff_tab(self):
            from PyQt5.QtWidgets import QLineEdit
            
            # 1. Asignamos el layout directamente a la pestaña (sin QSplitter ni left_widgets)
            layout = QVBoxLayout(self.tab_feff)
            layout.setContentsMargins(10, 10, 10, 10)
    
            # 2. Botones de Cargar Archivos (CIF o INP)
            h_load = QHBoxLayout()
            
            self.btn_load_cif = QPushButton("Cargar .CIF (DFT)")
            self.btn_load_cif.setStyleSheet("background-color: #8e44ad; color: white; font-weight: bold; padding: 10px; border-radius: 4px;")
            self.btn_load_cif.clicked.connect(self.cargar_cif)
            
            self.btn_load_inp = QPushButton("Cargar feff.inp")
            self.btn_load_inp.setStyleSheet("background-color: #2980b9; color: white; font-weight: bold; padding: 10px; border-radius: 4px;")
            self.btn_load_inp.clicked.connect(self.seleccionar_feff_inp)
            
            h_load.addWidget(self.btn_load_cif)
            h_load.addWidget(self.btn_load_inp)
            layout.addLayout(h_load)
            
            self.lbl_cif_path = QLabel("Ningún archivo cargado")
            self.lbl_cif_path.setStyleSheet("color: #7f8c8d; font-style: italic; margin-bottom: 5px;")
            self.lbl_cif_path.setWordWrap(True)
            layout.addWidget(self.lbl_cif_path)
    
            # 3. Parámetros Avanzados FEFF
            group_params = QFrame()
            group_params.setStyleSheet("QFrame { border: 1px solid #bdc3c7; border-radius: 6px; background-color: #f8f9fa; }")
            grid_params = QGridLayout(group_params)
            
            # Básico
            grid_params.addWidget(QLabel("<b>Átomo Target (Índice):</b>"), 0, 0)
            self.spin_target = QSpinBox(); self.spin_target.setRange(0, 999)
            self.spin_target.valueChanged.connect(self.actualizar_visor_3d)
            grid_params.addWidget(self.spin_target, 0, 1)
    
            grid_params.addWidget(QLabel("Borde (HOLE, 1=K):"), 1, 0)
            self.spin_edge = QSpinBox(); self.spin_edge.setValue(1); self.spin_edge.setRange(1, 4)
            grid_params.addWidget(self.spin_edge, 1, 1)
    
            grid_params.addWidget(QLabel("Radio Clúster (RPATH Å):"), 2, 0)
            self.spin_rpath = QDoubleSpinBox(); self.spin_rpath.setValue(6.0); self.spin_rpath.setSingleStep(0.5)
            grid_params.addWidget(self.spin_rpath, 2, 1)
    
            # Avanzado
            grid_params.addWidget(QLabel("Amplitud (S02):"), 3, 0)
            self.spin_feff_s02 = QDoubleSpinBox(); self.spin_feff_s02.setRange(0.1, 2.0); self.spin_feff_s02.setValue(1.0); self.spin_feff_s02.setSingleStep(0.1)
            grid_params.addWidget(self.spin_feff_s02, 3, 1)
    
            grid_params.addWidget(QLabel("Max scattering (NLEG):"), 4, 0)
            self.spin_nleg = QSpinBox(); self.spin_nleg.setRange(2, 8); self.spin_nleg.setValue(4)
            grid_params.addWidget(self.spin_nleg, 4, 1)
    
            grid_params.addWidget(QLabel("Filtro Paths (CRITERIA %):"), 5, 0)
            self.spin_crit = QDoubleSpinBox(); self.spin_crit.setRange(0.0, 100.0); self.spin_crit.setValue(4.0)
            grid_params.addWidget(self.spin_crit, 5, 1)
    
            grid_params.addWidget(QLabel("Límite EXAFS (k-max):"), 6, 0)
            self.spin_exafs_k = QDoubleSpinBox(); self.spin_exafs_k.setRange(5.0, 40.0); self.spin_exafs_k.setValue(20.0)
            grid_params.addWidget(self.spin_exafs_k, 6, 1)
    
            grid_params.addWidget(QLabel("DEBYE (T_exp, Theta_D):"), 7, 0)
            h_debye = QHBoxLayout()
            self.spin_temp = QDoubleSpinBox(); self.spin_temp.setRange(0, 1000); self.spin_temp.setValue(300); self.spin_temp.setDecimals(0)
            self.spin_theta = QDoubleSpinBox(); self.spin_theta.setRange(0, 1000); self.spin_theta.setValue(343); self.spin_theta.setDecimals(0)
            h_debye.addWidget(self.spin_temp); h_debye.addWidget(self.spin_theta)
            grid_params.addLayout(h_debye, 7, 1)
    
            grid_params.addWidget(QLabel("Self-Consistent (SCF):"), 8, 0)
            self.line_scf = QLineEdit("7.0 1 100 0.2 1") 
            grid_params.addWidget(self.line_scf, 8, 1)
    
            layout.addWidget(group_params)
    
            # 4. Acciones
            self.btn_gen_inp = QPushButton("2. Generar feff.inp")
            self.btn_gen_inp.setStyleSheet("background-color: #f39c12; color: white; font-weight: bold; padding: 10px; border-radius: 4px;")
            self.btn_gen_inp.clicked.connect(self.generar_feff_inp)
            self.btn_gen_inp.setEnabled(False)
            layout.addWidget(self.btn_gen_inp)
    
            self.btn_run_feff = QPushButton("3. Ejecutar FEFF (Calcular Paths)")
            self.btn_run_feff.setStyleSheet("background-color: #e74c3c; color: white; font-weight: bold; padding: 10px; border-radius: 4px;")
            self.btn_run_feff.clicked.connect(self.ejecutar_feff) 
            self.btn_run_feff.setEnabled(False)
            layout.addWidget(self.btn_run_feff)
    
            # 5. Consola
            self.consola_feff = QPlainTextEdit()
            self.consola_feff.setReadOnly(True)
            self.consola_feff.setStyleSheet("background-color: #1e1e1e; color: #2ecc71; font-family: Consolas, monospace; margin-top: 5px;")
            layout.addWidget(self.consola_feff)
    
            # Variables internas
            self.atomos_cif = []
            self.caja_a = 1.0

    def cargar_cif(self):
        archivo, _ = QFileDialog.getOpenFileName(self, "Selecciona archivo .cif", "", "Archivos CIF (*.cif);;Todos (*.*)")
        if not archivo: return
        
        self.lbl_cif_path.setText(os.path.basename(archivo))
        self.atomos_cif = []
        
        leyendo_atomos = False
        with open(archivo, 'r') as f:
            for linea in f:
                partes = linea.split()
                if not partes: continue
                
                if partes[0] == '_cell_length_a':
                    self.caja_a = float(partes[1])
                
                if partes[0] == '_atom_site_type_symbol':
                    leyendo_atomos = True
                    continue
                
                if leyendo_atomos and len(partes) >= 6:
                    try:
                        simbolo = partes[0]
                        elemento = ''.join([i for i in simbolo if not i.isdigit()])
                        
                        x = float(partes[2]) * self.caja_a
                        y = float(partes[3]) * self.caja_a
                        z = float(partes[4]) * self.caja_a
                        
                        self.atomos_cif.append({'elemento': elemento, 'x': x, 'y': y, 'z': z})
                    except ValueError:
                        break 
        
        if self.atomos_cif:
            self.consola_feff.setPlainText(f"CIF cargado correctamente.\nSe han detectado {len(self.atomos_cif)} átomos.\nElige tu átomo target en el panel de control.")
            self.spin_target.setRange(0, len(self.atomos_cif)-1)
            self.btn_gen_inp.setEnabled(True)
            self.feff_folder = os.path.dirname(archivo)
            self.actualizar_visor_3d()
            
    def leer_y_mostrar_paths_feff(self):
            # 1. Extraer Amplitudes (cw ratio) desde files.dat
            cw_ratios = {}
            ruta_files = os.path.join(self.feff_folder, 'files.dat')
            if os.path.exists(ruta_files):
                with open(ruta_files, 'r') as f:
                    for linea in f:
                        if 'feff' in linea and '.dat' in linea:
                            partes = linea.split()
                            if len(partes) >= 3:
                                idx_str = partes[0].replace('feff', '').replace('.dat', '')
                                if idx_str.isdigit():
                                    cw_ratios[int(idx_str)] = partes[2]
    
            # 2. Extraer los caminos de dispersión desde paths.dat (¡LA CLAVE ESTABA AQUÍ!)
            rutas_scattering = {}
            ruta_paths = os.path.join(self.feff_folder, 'paths.dat')
            if os.path.exists(ruta_paths):
                with open(ruta_paths, 'r') as f:
                    current_path_idx = None
                    leyendo_atomos = False
                    for linea in f:
                        # Buscar el inicio de un path (ej: "  1    2   1.000  index, nleg...")
                        if "index," in linea and "nleg," in linea:
                            partes = linea.split()
                            if partes[0].isdigit():
                                current_path_idx = int(partes[0])
                                rutas_scattering[current_path_idx] = []
                                leyendo_atomos = False
                            continue
    
                        if current_path_idx is not None:
                            # Buscar los encabezados de la tabla de átomos
                            if "ipot" in linea and "label" in linea:
                                leyendo_atomos = True
                                continue
    
                            if leyendo_atomos:
                                partes = linea.split()
                                if len(partes) >= 5:
                                    try:
                                        # Validamos que la primera columna sea un número (coordenada X)
                                        float(partes[0]) 
                                        # Extraemos el átomo quitando las comillas que pone FEFF (Ej: "'N   '" -> "N")
                                        atomo = partes[4].replace("'", "").strip()
                                        if atomo:
                                            rutas_scattering[current_path_idx].append(atomo)
                                    except ValueError:
                                        # Si no es un número, se acabó la tabla de átomos
                                        leyendo_atomos = False
                                else:
                                    leyendo_atomos = False
    
            # 3. Leer feffNNNN.dat para sacar distancia exacta y poblar la tabla
            archivos_feff = glob.glob(os.path.join(self.feff_folder, "feff[0-9][0-9][0-9][0-9].dat"))
            archivos_feff.sort()
            self.tabla_paths.setRowCount(0)
    
            for archivo in archivos_feff:
                nombre = os.path.basename(archivo)
                idx_path = int(nombre.replace("feff", "").replace(".dat", "")) if "feff" in nombre else -1
                
                porcentaje = cw_ratios.get(idx_path, "N/A")
                scattering_list = rutas_scattering.get(idx_path, [])
                
                # Formatear el string al estilo Artemis (Absorbedor -> Scatterer -> Absorbedor)
                if scattering_list:
                    absorbedor = scattering_list[-1] # En paths.dat el último átomo siempre es el retorno al central
                    ruta_str = f"{absorbedor} \u2192 " + " \u2192 ".join(scattering_list)
                else:
                    ruta_str = "Desconocido"
    
                degen, reff = "N/A", "N/A"
    
                with open(archivo, 'r') as f:
                    for linea in f:
                        if "nleg," in linea and "deg," in linea and "reff," in linea:
                            partes = linea.split()
                            try:
                                degen = partes[1]
                                reff = partes[2]
                            except IndexError: pass
                            break # Optimización: Paramos de leer para no tragarnos todo el archivo numérico
    
                row = self.tabla_paths.rowCount()
                self.tabla_paths.insertRow(row)
                self.tabla_paths.setItem(row, 0, QTableWidgetItem(nombre))
                self.tabla_paths.setItem(row, 1, QTableWidgetItem(f"{porcentaje} %"))
                self.tabla_paths.setItem(row, 2, QTableWidgetItem(degen))
                self.tabla_paths.setItem(row, 3, QTableWidgetItem(reff))
                self.tabla_paths.setItem(row, 4, QTableWidgetItem(ruta_str))
                
    def exportar_grafica_txt(self):
            # 1. Comprobar si hay gráficas dibujadas en el lienzo
            if not self.fig.axes:
                QMessageBox.warning(self, "Aviso", "No hay ninguna gráfica en pantalla para exportar.")
                return
    
            # 2. Pedir al usuario dónde guardar el archivo
            ruta, _ = QFileDialog.getSaveFileName(self, "Guardar Datos como TXT", "Datos_Exportados.txt", "Archivo de Texto (*.txt);;Archivo CSV (*.csv);;Todos (*.*)")
            if not ruta:
                return
    
            # 3. Extraer todas las curvas directamente de Matplotlib
            df_export = pd.DataFrame()
            lineas_encontradas = 0
            
            for ax in self.fig.axes:
                for linea in ax.lines:
                    label = linea.get_label()
                    
                    # Ignoramos líneas decorativas ocultas de Matplotlib (suelen empezar por un guión bajo)
                    if label.startswith('_'):
                        continue
                        
                    # Si la curva no tiene nombre, le inventamos uno
                    if not label:
                        label = f"Curva_{lineas_encontradas}"
    
                    # Extraemos los datos físicos
                    x_data = linea.get_xdata()
                    y_data = linea.get_ydata()
    
                    # Usamos pd.Series para evitar errores si las curvas tienen diferentes longitudes
                    df_export[f"{label}_X"] = pd.Series(x_data)
                    df_export[f"{label}_Y"] = pd.Series(y_data)
                    lineas_encontradas += 1
    
            if lineas_encontradas == 0:
                QMessageBox.information(self, "Aviso", "La gráfica actual no contiene curvas exportables.")
                return
    
            # 4. Guardar a disco usando tabulaciones (ideal para pegar en Origin o Excel)
            try:
                df_export.to_csv(ruta, sep='\t', index=False)
                QMessageBox.information(self, "Éxito", f"Se han exportado {lineas_encontradas} curvas correctamente en:\n{ruta}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Fallo al guardar el archivo:\n{e}")  
                
    def actualizar_visor_3d(self):
            if not self.atomos_cif: return
            
            self.fig.clear() # ¡Borramos las gráficas 2D!
            self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            ax = self.fig.add_subplot(111, projection='3d')
            ax.set_facecolor('#2c3e50')
            
            target_idx = self.spin_target.value()
            
            colores_cpk = {'Cu': '#d35400', 'N': '#2980b9', 'O': '#c0392b', 'C': '#95a5a6', 'H': '#ecf0f1'}
            tamanos_cpk = {'Cu': 500, 'N': 250, 'O': 250, 'C': 250, 'H': 100}
            
            coords = np.array([[a['x'], a['y'], a['z']] for a in self.atomos_cif])
            
            dist_corte_normal = 1.7 
            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    dist = np.linalg.norm(coords[i] - coords[j])
                    es_cobre = self.atomos_cif[i]['elemento'] == 'Cu' or self.atomos_cif[j]['elemento'] == 'Cu'
                    limite = 2.5 if es_cobre else dist_corte_normal
                    
                    if 0.5 < dist < limite:
                        ax.plot([coords[i, 0], coords[j, 0]], [coords[i, 1], coords[j, 1]], [coords[i, 2], coords[j, 2]], 
                                color='#bdc3c7', linewidth=4, alpha=1.0, zorder=1, solid_capstyle='round')
    
            xs, ys, zs, colores, sizes = [], [], [], [], []
            for i, atomo in enumerate(self.atomos_cif):
                xs.append(atomo['x'])
                ys.append(atomo['y'])
                zs.append(atomo['z'])
                if i == target_idx:
                    colores.append('#2ecc71')
                    sizes.append(700) 
                    ax.text(atomo['x'], atomo['y'], atomo['z'] + 0.4, f" Target: {atomo['elemento']}{i}", 
                            color='white', fontweight='bold', fontsize=12, zorder=10)
                else:
                    colores.append(colores_cpk.get(atomo['elemento'], '#7f8c8d'))
                    sizes.append(tamanos_cpk.get(atomo['elemento'], 150))
    
            ax.scatter(xs, ys, zs, c=colores, s=sizes, edgecolors='#1a252f', linewidths=1.5, alpha=1.0, depthshade=True, zorder=2)
    
            max_range = np.array([coords[:,0].max()-coords[:,0].min(), coords[:,1].max()-coords[:,1].min(), coords[:,2].max()-coords[:,2].min()]).max() / 2.0
            mid_x = (coords[:,0].max()+coords[:,0].min()) * 0.5
            mid_y = (coords[:,1].max()+coords[:,1].min()) * 0.5
            mid_z = (coords[:,2].max()+coords[:,2].min()) * 0.5
            
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)
            
            try: ax.set_box_aspect((1, 1, 1))
            except AttributeError: pass
    
            ax.axis('off')
            self.canvas.draw()
    def generar_feff_inp(self):
            if not self.atomos_cif: return
            
            target_idx = self.spin_target.value()
            target_atom = self.atomos_cif[target_idx]
            
            # Recogemos todas las variables de la interfaz
            radio = self.spin_rpath.value()
            hole = self.spin_edge.value()
            s02 = self.spin_feff_s02.value()
            nleg = self.spin_nleg.value()
            crit = self.spin_crit.value()
            kmax = self.spin_exafs_k.value()
            temp = int(self.spin_temp.value())
            theta = int(self.spin_theta.value())
            scf = self.line_scf.text().strip()
            
            elementos_unicos = []
            for a in self.atomos_cif:
                if a['elemento'] not in elementos_unicos and a['elemento'] != target_atom['elemento']:
                    elementos_unicos.append(a['elemento'])
                    
            potenciales = {target_atom['elemento']: 0}
            for i, el in enumerate(elementos_unicos):
                potenciales[el] = i + 1
    
            atomos_exportar = []
            for i, a in enumerate(self.atomos_cif):
                dist = math.sqrt((a['x'] - target_atom['x'])**2 + (a['y'] - target_atom['y'])**2 + (a['z'] - target_atom['z'])**2)
                if dist <= radio:
                    atomos_exportar.append({
                        'x': a['x'], 'y': a['y'], 'z': a['z'],
                        'ipot': 0 if i == target_idx else potenciales[a['elemento']],
                        'elemento': a['elemento'],
                        'dist': dist
                    })
                    
            atomos_exportar.sort(key=lambda x: x['dist'])
    
            self.feff_filename = "feff.inp"
            ruta_guardado = os.path.join(self.feff_folder, self.feff_filename)
            
            # Escribimos el archivo feff.inp con los parámetros avanzados
            with open(ruta_guardado, 'w') as f:
                f.write(f"TITLE Generado por Python - Absorbedor: {target_atom['elemento']}{target_idx}\n\n")
                
                f.write("* --- Configuracion del Espectro ---\n")
                f.write(f"HOLE {hole}   1.0\n")
                f.write(f"S02 {s02}\n\n")
                
                f.write("CONTROL 1 1 1 1 1 1\n")
                f.write("PRINT   1 0 0 0 0 3\n\n")
                
                f.write("* --- Parametros de EXAFS ---\n")
                f.write(f"NLEG {nleg}\n")
                # CRITERIA suele usar dos valores (curved wave y plane wave). Ponemos el mismo a ambos.
                f.write(f"CRITERIA {crit} {crit}\n") 
                
                if scf: # Si la caja de texto no está vacía, escribimos SCF
                    f.write(f"SCF {scf}\n")
                    
                f.write(f"DEBYE {temp} {theta} 0\n")
                f.write(f"EXAFS {kmax}\n")
                f.write(f"RPATH {radio}\n\n")
                
                f.write("POTENTIALS\n")
                f.write(f"0  {self._get_z(target_atom['elemento'])}  {target_atom['elemento']}\n")
                for el in elementos_unicos:
                    f.write(f"{potenciales[el]}  {self._get_z(el)}  {el}\n")
                    
                f.write("\nATOMS\n")
                for a in atomos_exportar:
                    f.write(f"  {a['x']:>10.5f} {a['y']:>10.5f} {a['z']:>10.5f}   {a['ipot']}  {a['elemento']}\n")
    
            self.consola_feff.appendPlainText("\n--- feff.inp GENERADO ---")
            self.consola_feff.appendPlainText(f"Ruta: {ruta_guardado}")
            self.consola_feff.appendPlainText(f"Átomos calculados: {len(atomos_exportar)}")
            self.btn_run_feff.setEnabled(True)
    def _get_z(self, elemento):
        tabla = {'H':1, 'C':6, 'N':7, 'O':8, 'F':9, 'P':15, 'S':16, 'Cl':17, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30}
        return tabla.get(elemento, 0)
    def init_fit_tab(self):
            layout = QVBoxLayout(self.tab_fit)
            
            # --- 1. PARÁMETROS GLOBALES (Arriba) ---
            group_global = QFrame(); group_global.setStyleSheet("QFrame { background-color: #ecf0f1; border-radius: 5px; }")
            grid_global = QGridLayout(group_global)
            
            grid_global.addWidget(QLabel("<b>Parámetros Globales (Compartidos):</b>"), 0, 0, 1, 3)
            
            grid_global.addWidget(QLabel("S02 (Amplitud):"), 1, 0)
            self.spin_s02 = QDoubleSpinBox(); self.spin_s02.setRange(0.1, 2.0); self.spin_s02.setValue(0.85); self.spin_s02.setSingleStep(0.05)
            self.chk_s02 = QCheckBox("Fitear"); self.chk_s02.setChecked(True)
            grid_global.addWidget(self.spin_s02, 1, 1); grid_global.addWidget(self.chk_s02, 1, 2)
    
            grid_global.addWidget(QLabel("ΔE0 (Shift eV):"), 2, 0)
            self.spin_de0 = QDoubleSpinBox(); self.spin_de0.setRange(-15.0, 15.0); self.spin_de0.setValue(0.0)
            self.chk_de0 = QCheckBox("Fitear"); self.chk_de0.setChecked(True)
            grid_global.addWidget(self.spin_de0, 2, 1); grid_global.addWidget(self.chk_de0, 2, 2)
    
            grid_global.addWidget(QLabel("Rango Rmin - Rmax:"), 3, 0)
            h_r = QHBoxLayout()
            self.spin_rmin = QDoubleSpinBox(); self.spin_rmin.setValue(1.0)
            self.spin_rmax = QDoubleSpinBox(); self.spin_rmax.setValue(3.0)
            h_r.addWidget(self.spin_rmin); h_r.addWidget(QLabel("-")); h_r.addWidget(self.spin_rmax)
            grid_global.addLayout(h_r, 3, 1, 1, 2)
    
            layout.addWidget(group_global)
    
            # --- 2. GESTOR DINÁMICO DE PATHS (Medio) ---
            h_botones_path = QHBoxLayout()
            
            self.btn_add_path = QPushButton("+ Añadir Path Manual")
            self.btn_add_path.setStyleSheet("font-weight: bold; padding: 8px; background-color: #3498db; color: white; border-radius: 4px;")
            self.btn_add_path.clicked.connect(self.agregar_path_dinamico)
            
            self.btn_auto_shell = QPushButton(" Autodetectar 1ª Capa")
            self.btn_auto_shell.setStyleSheet("font-weight: bold; padding: 8px; background-color: #f1c40f; color: black; border-radius: 4px;")
            self.btn_auto_shell.clicked.connect(self.autodetectar_primera_capa)
            
            h_botones_path.addWidget(self.btn_add_path)
            h_botones_path.addWidget(self.btn_auto_shell)
            layout.addLayout(h_botones_path)
            # Aquí guardaremos la información de los widgets de cada path
            self.lista_paths_ui = []
    
            # Área con Scroll para meter infinitos paths
            self.scroll_area = QScrollArea()
            self.scroll_area.setWidgetResizable(True)
            self.scroll_widget = QWidget()
            self.scroll_layout = QVBoxLayout(self.scroll_widget)
            self.scroll_layout.setAlignment(Qt.AlignTop)
            self.scroll_area.setWidget(self.scroll_widget)
            layout.addWidget(self.scroll_area, stretch=1)
    
            # --- 3. EJECUCIÓN Y RESULTADOS (Abajo) ---
            self.btn_fitear = QPushButton("¡EJECUTAR fit MULTI-PATH!")
            self.btn_fitear.setStyleSheet("font-weight: bold; padding: 12px; background-color: #27ae60; color: white;")
            self.btn_fitear.clicked.connect(self.ejecutar_fit)
            self.btn_fitear.setEnabled(False) # Se activa al añadir el primer path
            layout.addWidget(self.btn_fitear)
    
            self.consola_fit = QPlainTextEdit()
            self.consola_fit.setReadOnly(True)
            self.consola_fit.setStyleSheet("background-color: #fdf6e3; color: #000000; font-family: monospace;")
            layout.addWidget(self.consola_fit, stretch=1)

    # --- FUNCIONES AUXILIARES PARA CREAR Y DESTRUIR PATHS ---
    def agregar_path_dinamico(self,archivo=None):
        
        if not isinstance(archivo, str): 
            archivo, _ = QFileDialog.getOpenFileName(self, "Selecciona Path", "", "Archivos DAT (feff*.dat);;Todos (*.*)")
        if not archivo: return

        # Crear una cajita (Frame) para este path
        frame = QFrame()
        frame.setStyleSheet("QFrame { border: 1px solid #bdc3c7; border-radius: 4px; margin-bottom: 5px; padding: 5px; }")
        f_layout = QGridLayout(frame)

        num_path = len(self.lista_paths_ui) + 1
        lbl_nombre = QLabel(f"<b>Path {num_path}:</b> {os.path.basename(archivo)}")
        lbl_nombre.setStyleSheet("color: #2980b9; border: none;")
        
        # Botón para eliminar este path específico
        btn_eliminar = QPushButton("X")
        btn_eliminar.setStyleSheet("background-color: #e74c3c; color: white; font-weight: bold; border: none; max-width: 30px;")

        # Controles locales (N, DR y Sigma2)
        spin_n = QDoubleSpinBox(); spin_n.setRange(0.1, 20.0); spin_n.setValue(1.0); spin_n.setSingleStep(0.5)
        chk_n = QCheckBox("Fitear"); chk_n.setChecked(True)
        
        spin_dr = QDoubleSpinBox(); spin_dr.setRange(-0.5, 0.5); spin_dr.setValue(0.0); spin_dr.setSingleStep(0.01)
        chk_dr = QCheckBox("Fitear"); chk_dr.setChecked(True)
        
        spin_sig2 = QDoubleSpinBox(); spin_sig2.setRange(0.001, 0.05); spin_sig2.setValue(0.005); spin_sig2.setDecimals(4)
        chk_sig2 = QCheckBox("Fitear"); chk_sig2.setChecked(True)

        # Organizar en el grid del frame
        f_layout.addWidget(lbl_nombre, 0, 0, 1, 2)
        f_layout.addWidget(btn_eliminar, 0, 2)
        
        # FILA 1: N (Coord)
        f_layout.addWidget(QLabel("N (Coord):"), 1, 0)
        f_layout.addWidget(spin_n, 1, 1); f_layout.addWidget(chk_n, 1, 2)

        # FILA 2: Delta R  <--- ¡Asegúrate de que este es un 2!
        f_layout.addWidget(QLabel("ΔR (Å):"), 2, 0)
        f_layout.addWidget(spin_dr, 2, 1); f_layout.addWidget(chk_dr, 2, 2)
        
        # FILA 3: Sigma2 <--- ¡Asegúrate de que este es un 3!
        f_layout.addWidget(QLabel("σ²:"), 3, 0)
        f_layout.addWidget(spin_sig2, 3, 1); f_layout.addWidget(chk_sig2, 3, 2)

        self.scroll_layout.addWidget(frame)

        # Guardar todas las referencias en un diccionario para leerlas al fitear
        path_data = {
                    'frame': frame,
                    'archivo': archivo,
                    'spin_n': spin_n, 'chk_n': chk_n,  # <-- NUEVO
                    'spin_dr': spin_dr, 'chk_dr': chk_dr,
                    'spin_sig2': spin_sig2, 'chk_sig2': chk_sig2
                }
        self.lista_paths_ui.append(path_data)
        
        # Conectar el botón de eliminar
        btn_eliminar.clicked.connect(lambda: self.eliminar_path_dinamico(path_data))
        
        self.btn_fitear.setEnabled(True)
    def autodetectar_primera_capa(self):
            # 1. Comprobar si hemos calculado FEFF antes
            if not hasattr(self, 'feff_folder') or not self.feff_folder:
                QMessageBox.warning(self, "Aviso", "Primero debes ejecutar FEFF (Pestaña 5) para generar los paths.")
                return
                
            # 2. Leemos la ventana de distancia de la interfaz
            rmin = self.spin_rmin.value()
            rmax = self.spin_rmax.value()
            
            archivos_feff = glob.glob(os.path.join(self.feff_folder, "feff[0-9][0-9][0-9][0-9].dat"))
            paths_encontrados = 0
            
            # 3. Escaneamos los archivos
            for archivo in archivos_feff:
                with open(archivo, 'r') as f:
                    for linea in f:
                        # Buscamos la línea de la cabecera que tiene la distancia "reff"
                        if "nleg," in linea and "reff," in linea:
                            partes = linea.split()
                            try:
                                reff = float(partes[2])
                                
                                # MAGIA: Si el path cae dentro de Rmin y Rmax, lo añadimos
                                if rmin <= reff <= rmax:
                                    self.agregar_path_dinamico(archivo)
                                    paths_encontrados += 1
                                    
                            except (ValueError, IndexError):
                                pass
                                
                            # Una vez leída la cabecera, pasamos al siguiente archivo para ir súper rápido
                            break 
                            
            # 4. Avisar al usuario del resultado
            if paths_encontrados > 0:
                self.consola_fiteo.appendPlainText(f"\n⚡ Autodetección: Se han añadido {paths_encontrados} paths en el rango de {rmin} a {rmax} Å.")
            else:
                QMessageBox.information(self, "Resultado", f"No se encontraron paths en el rango de {rmin} a {rmax} Å.\nPrueba a subir el Rmax o a generar los paths de FEFF de nuevo.")
    def eliminar_path_dinamico(self, path_data):
        # Borrar de la interfaz
        path_data['frame'].deleteLater()
        # Borrar de la lista lógica
        self.lista_paths_ui.remove(path_data)
        # Desactivar fit si no quedan paths
        if len(self.lista_paths_ui) == 0:
            self.btn_fitear.setEnabled(False) 
        
    # ---------------------------------------------------------
    # MOTOR DE CÁLCULO EXAFS (NUEVO PIPELINE CENTRALIZADO)
    # ---------------------------------------------------------
    def procesar_xas_pipeline(self):
            if not self.datos_promedio_actual: return 
    
            energia_original = np.array(self.datos_promedio_actual['x'], dtype=float)
            if np.max(energia_original) < 1000: energia_larch = energia_original * 1000.0
            else: energia_larch = energia_original
    
            self.datos_xas = Group(energy=energia_larch, mu=np.array(self.datos_promedio_actual['y'], dtype=float))
            
            e0 = self.spin_e0.value() if self.spin_e0.value() > 0 else None 
            try:
                pre_edge(self.datos_xas, e0=e0, pre1=self.spin_pre1.value(), pre2=self.spin_pre2.value(), 
                         norm1=self.spin_norm1.value(), norm2=self.spin_norm2.value())
                
                autobk(self.datos_xas, rbkg=self.spin_rbkg.value(), kmin=self.spin_kmin.value(), 
                       kmax=self.spin_kmax.value(), kweight=int(self.spin_kweight.value()))
                
                
                xftf(self.datos_xas, 
                     kweight=int(self.spin_kweight.value()),
                     kmin=self.spin_ft_kmin.value(),
                     kmax=self.spin_ft_kmax.value(),
                     dk=self.spin_ft_dk.value(),
                     window=self.combo_window.currentText())
                
            except Exception as e:
                print(f"Error en procesado Larch: {e}")
                return
    
            idx = self.tabs.currentIndex()
            if idx == 1: self.dibujar_pestana_exafs()
            elif idx == 2: self.dibujar_pestana_autobk()
            elif idx == 3: self.dibujar_pestana_ft() 
            
    def actualizar_vista_segun_pestana(self):
            idx = self.tabs.currentIndex()
            
            # 1. Control de la visibilidad de la tabla de FEFF
            if hasattr(self, 'tabla_paths'):

                self.tabla_paths.setVisible(idx == 4)
    
            if idx == 0:
                self.procesar_y_plotear_explorador()
            elif idx == 4:
                self.actualizar_visor_3d() 
            elif idx == 5:
                self.procesar_xas_pipeline()
            elif idx in [1, 2, 3]:
                self.procesar_xas_pipeline()

    # ---------------------------------------------------------
    # DIBUJO DE GRÁFICOS
    # ---------------------------------------------------------
    def dibujar_pestana_exafs(self):
        """Dibuja los resultados de la Normalización (Pestaña 2)"""
        self.fig.clear()
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)

        ax1.plot(self.datos_xas.energy, self.datos_xas.mu, label='Raw $\mu$', color='blue')
        ax1.plot(self.datos_xas.energy, self.datos_xas.pre_edge, label='Pre-edge', color='red', linestyle='--')
        ax1.plot(self.datos_xas.energy, self.datos_xas.post_edge, label='Post-edge', color='green', linestyle='--')
        ax1.axvline(self.datos_xas.e0, color='black', linestyle=':', label=f'E0 = {self.datos_xas.e0:.1f} eV')
        ax1.legend(); ax1.set_ylabel("$\mu(E)$")

        ax2.plot(self.datos_xas.energy, self.datos_xas.flat, label='Normalizado', color='black', lw=2)
        ax2.axhline(0, color='gray', linestyle='--'); ax2.axhline(1, color='gray', linestyle='--')
        ax2.set_xlabel("Energía (eV)"); ax2.set_ylabel("Norm $\mu(E)$")

        self.fig.tight_layout()
        self.canvas.draw()

    def dibujar_pestana_autobk(self):
            """Dibuja los resultados del Background y el espacio K (Pestaña 3)"""
            self.fig.clear()
            
            # Arriba: Señal RAW + Curva de Background RAW (Spline)
            ax1 = self.fig.add_subplot(211)
            
            # ¡CORRECCIÓN AQUÍ! Usamos 'mu' (crudo) en vez de 'flat' (normalizado)
            ax1.plot(self.datos_xas.energy, self.datos_xas.mu, label='Señal Cruda ($\mu$)', color='black', lw=2)
            ax1.plot(self.datos_xas.energy, self.datos_xas.bkg, label='Background (Spline)', color='red', lw=1.5)
            
            ax1.set_title(f"Ajuste de Background (Rbkg = {self.spin_rbkg.value()})")
            ax1.set_ylabel("$\mu(E)$ (Escala Cruda)")
            ax1.legend()
    
            # Abajo: La oscilación pura en el espacio k (Esto estaba bien)
            ax2 = self.fig.add_subplot(212)
            k = self.datos_xas.k
            chi = self.datos_xas.chi
            kw = int(self.spin_kweight.value())
            
            ax2.plot(k, chi * (k**kw), color='blue', lw=2, label=f'$k^{kw}\chi(k)$')
            ax2.axhline(0, color='gray', linestyle='--')
            ax2.set_title("Señal EXAFS (Espacio $k$)")
            ax2.set_xlabel("$k \ (\AA^{-1})$")
            ax2.set_ylabel(f"$k^{kw}\chi(k)$")
            ax2.legend()
    
            self.fig.tight_layout()
            self.canvas.draw()
            
    def dibujar_pestana_ft(self):
            """Dibuja los resultados de la Transformada de Fourier (Pestaña 4)"""
            self.fig.clear()
            
            # Gráfico Superior: Mostrar la ventana superpuesta en la señal k
            ax1 = self.fig.add_subplot(211)
            k = self.datos_xas.k
            chi = self.datos_xas.chi
            kw = int(self.spin_kweight.value())
            chi_kw = chi * (k**kw)
    
            ax1.plot(k, chi_kw, label=f'$k^{kw}\chi(k)$', color='blue', lw=1.5)
            
            # Dibujamos la "ventana" que Larch ha usado. 
            # La multiplicamos por el máximo de nuestra señal para que se vea en la misma escala visual
            ventana_visual = self.datos_xas.kwin * np.max(np.abs(chi_kw))
            ax1.plot(k, ventana_visual, label=f'Ventana ({self.combo_window.currentText()})', color='red', linestyle='--', lw=2)
            
            ax1.set_title("Espacio $k$ y Ventana de Fourier")
            ax1.set_xlabel("$k \ (\AA^{-1})$")
            ax1.set_ylabel(f"$k^{kw}\chi(k)$")
            ax1.legend(loc='upper right')
    
            # Gráfico Inferior: El Espacio R
            ax2 = self.fig.add_subplot(212)
            r = self.datos_xas.r
            chir_mag = self.datos_xas.chir_mag # La magnitud de la transformada
    
            ax2.plot(r, chir_mag, label='$|\chi(R)|$ (Magnitud)', color='black', lw=2)
            
            ax2.set_title("Transformada de Fourier (Espacio $R$)")
            ax2.set_xlabel("$R \ (\AA)$")
            ax2.set_ylabel(f"$|\chi(R)| \ (\AA^{{-{kw+1}}})$")
            ax2.set_xlim(0, 6) # Normalmente el EXAFS interesante está entre 0 y 6 Angstroms
            ax2.legend()
    
            self.fig.tight_layout()
            self.canvas.draw()
            
            # --- EL MOTOR DEL fit DINÁMICO ---
    def ejecutar_fit(self):
            if len(self.lista_paths_ui) == 0 or not self.datos_xas:
                QMessageBox.warning(self, "Error", "Faltan datos XAS o no has añadido ningún Path.")
                return
    
            self.consola_fit.setPlainText(f"Iniciando fit con {len(self.lista_paths_ui)} paths...\n")
    
            # 1. Parámetros Globales
            param_amp = param(self.spin_s02.value(), vary=self.chk_s02.isChecked())
            param_de0 = param(self.spin_de0.value(), vary=self.chk_de0.isChecked())
            
            dict_params = {'amp': param_amp, 'del_e0': param_de0}
            lista_larch_paths = []
    
            # 2. Bucle Dinámico de Paths
            for i, path_data in enumerate(self.lista_paths_ui):
                nombre_n = f'n_{i}'
                nombre_dr = f'del_r_{i}'
                nombre_sig2 = f'sig2_{i}'
    
                dict_params[nombre_n] = param(path_data['spin_n'].value(), vary=path_data['chk_n'].isChecked())
                dict_params[nombre_dr] = param(path_data['spin_dr'].value(), vary=path_data['chk_dr'].isChecked())
                dict_params[nombre_sig2] = param(path_data['spin_sig2'].value(), vary=path_data['chk_sig2'].isChecked())
    
                feff_path = FeffPathGroup(path_data['archivo'], 
                                          s02=f'{nombre_n} * amp',
                                          e0='del_e0',         
                                          deltar=nombre_dr,    
                                          sigma2=nombre_sig2)  
                lista_larch_paths.append(feff_path)
    
            params_fit = Group(**dict_params)
    
            # 3. Transformada y fit
            kw = int(self.spin_kweight.value())
            trans = feffit_transform(kmin=self.spin_ft_kmin.value(), kmax=self.spin_ft_kmax.value(), 
                                     kweight=kw, dk=self.spin_ft_dk.value(), 
                                     window=self.combo_window.currentText(),
                                     rmin=self.spin_rmin.value(), rmax=self.spin_rmax.value())
    
            dataset = feffit_dataset(data=self.datos_xas, pathlist=lista_larch_paths, transform=trans)
            
            try:
                resultado = feffit(params_fit, dataset)
            except Exception as e:
                self.consola_fit.appendPlainText(f"\nError matemático: {e}")
                return
    
            self.consola_fit.appendPlainText(feffit_report(resultado))
    
            # 4. GRÁFICOS ACTUALIZADOS
            self.fig.clear()
            
            # Gráfico Superior: Espacio R (Magnitud) - Se mantiene igual
            ax1 = self.fig.add_subplot(211)
            r = dataset.data.r
            ax1.plot(r, dataset.data.chir_mag, label='Datos Exp.', color='black', lw=2)
            ax1.plot(r, dataset.model.chir_mag, label='fit Total', color='red', lw=2, linestyle='--')
            ax1.axvspan(self.spin_rmin.value(), self.spin_rmax.value(), color='gray', alpha=0.15)
            ax1.set_title("fit en Espacio $R$ (Magnitud)")
            ax1.set_ylabel(r"$|\chi(R)|$")
            ax1.set_xlim(0, 6)
            ax1.legend(fontsize='small')
    
            # Gráfico Inferior: Espacio K (Oscilaciones) - ¡NUEVO!
            ax2 = self.fig.add_subplot(212)
            k = dataset.data.k
            # Aplicamos el k-weight a los datos y al modelo para que se vea como en la pestaña 3
            chi_exp = dataset.data.chi * (k**kw)
            chi_mod = dataset.model.chi * (k**kw)
            
            ax2.plot(k, chi_exp, label=f'Datos $k^{kw}\chi(k)$', color='blue', lw=1.5)
            ax2.plot(k, chi_mod, label='fit Modelo', color='red', lw=1.5, linestyle='--')
            
            # Sombreado de la ventana de fit en K
            ax2.axvspan(self.spin_ft_kmin.value(), self.spin_ft_kmax.value(), color='gray', alpha=0.15, label='Ventana Fourier')
            
            ax2.set_title(f"fit en Espacio $k$ (weight={kw})")
            ax2.set_xlabel(r"$k \ (\AA^{-1})$")
            ax2.set_ylabel(r"$k^{%d}\chi(k)$" % kw)
            ax2.set_xlim(0, max(k))
            ax2.axhline(0, color='black', lw=0.5, alpha=0.5)
            ax2.legend(fontsize='small', loc='upper right')
    
            self.fig.tight_layout()
            self.canvas.draw()
    # ---------------------------------------------------------
    # FUNCIONES DE LECTURA (Ocultas por brevedad, son idénticas)
    # ---------------------------------------------------------
    def leer_archivo_inteligente(self, filepath):
        ext = os.path.splitext(filepath)[1].lower()
        try:
            if ext in ['.xls', '.xlsx']: return pd.read_excel(filepath)
            elif ext == '.csv': return pd.read_csv(filepath)
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f.readlines():
                    if line.startswith('#L'):
                        return pd.read_csv(filepath, sep=r'\s+', comment='#', header=None, names=line.replace('#L', '').strip().split(), engine='python')
            return pd.read_csv(filepath, sep=None, engine='python')
        except: return None
        
    def seleccionar_feff_inp(self):
            archivo, _ = QFileDialog.getOpenFileName(self, "Selecciona archivo feff.inp", "", "Archivos INP (*.inp);;Todos (*.*)")
            if not archivo: return
    
            # Guardamos la ruta para cuando le demos al botón rojo de Ejecutar
            self.feff_folder = os.path.dirname(archivo)
            self.feff_filename = os.path.basename(archivo)
            
            self.lbl_cif_path.setText(self.feff_filename)
            self.consola_feff.setPlainText(f"--- ARCHIVO {self.feff_filename} CARGADO ---\nRuta: {archivo}\nListo para ejecutar FEFF.")
            
            # ¡Magia! Vamos a leer el INP para dibujarlo en el 3D
            self.atomos_cif = []
            leyendo_atomos = False
             
            try:
                with open(archivo, 'r') as f:
                    for linea in f:
                        linea = linea.strip()
                        # Ignorar líneas vacías o comentarios (en FEFF empiezan por *)
                        if not linea or linea.startswith('*'): 
                            continue
                            
                        if linea.upper().startswith('ATOMS'):
                            leyendo_atomos = True
                            continue
                            
                        if leyendo_atomos:
                            partes = linea.split()
                            # Un archivo FEFF tiene: x, y, z, ipot, tag (5 columnas min)
                            if len(partes) >= 5:
                                x, y, z = float(partes[0]), float(partes[1]), float(partes[2])
                                tag = partes[4]
                                # Limpiar números (ej: Cu1 -> Cu) para que el visor sepa el color
                                elemento = ''.join([i for i in tag if not i.isdigit()])
                                self.atomos_cif.append({'elemento': elemento, 'x': x, 'y': y, 'z': z})
            except Exception as e:
                self.consola_feff.appendPlainText(f"\nAviso: No se pudo dibujar en 3D ({e})")
    
            # Actualizar la interfaz
            if self.atomos_cif:
                self.spin_target.setRange(0, len(self.atomos_cif)-1)
                self.spin_target.setValue(0)
                self.actualizar_visor_3d()
    
            self.btn_run_feff.setEnabled(True)
            # Desactivamos el botón naranja porque ya tienes un INP hecho
            self.btn_gen_inp.setEnabled(False)
            
    def cargar_archivos(self):
        archivos, _ = QFileDialog.getOpenFileNames(self, "Selecciona", "", "Todos (*.*)")
        if not archivos: return
        self.bloquear_plot = True; self.data_cache = {}
        for p in archivos:
            df = self.leer_archivo_inteligente(p)
            if df is not None: self.data_cache[p] = df.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
        if not self.data_cache: return
        self.columnas_disponibles = list(list(self.data_cache.values())[0].columns)
        self.combo_x.clear(); self.lista_y.clear()
        self.combo_x.addItems(self.columnas_disponibles); self.lista_y.addItems(self.columnas_disponibles)
        self.bloquear_plot = False; self.procesar_y_plotear_explorador()

    def procesar_y_plotear_explorador(self):
            if self.bloquear_plot or not self.data_cache: return
            
            col_x = self.combo_x.currentText()
            cols_y = [item.text() for item in self.lista_y.selectedItems()]
            if not col_x or not cols_y: return
    
            # 1. Limpieza total y preparación del lienzo
            self.fig.clear()
            # Aseguramos que la tabla de FEFF esté oculta si volvemos a Datos
            if hasattr(self, 'tabla_paths'):
                self.tabla_paths.setVisible(False)
    
            # 2. Lógica de Subplots (Recuperamos la diferencia)
            if self.check_diferencia.isChecked() and len(cols_y) >= 2:
                self.ax_main = self.fig.add_subplot(211)
                self.ax_diff = self.fig.add_subplot(212, sharex=self.ax_main)
            else:
                self.ax_main = self.fig.add_subplot(111)
                self.ax_diff = None
    
            colores = plt.rcParams['axes.prop_cycle'].by_key()['color']
            res_medias = {}
            
            # 3. Dibujo de señales
            for i, y_t in enumerate(cols_y):
                all_y, x_ref = [], None
                c = colores[i % len(colores)]
                for path, df in self.data_cache.items():
                    if col_x not in df.columns or y_t not in df.columns: continue
                    t = df[[col_x, y_t]].dropna()
                    x_d, y_d = t[col_x].values, t[y_t].values
                    if col_x == "dutd": x_d = x_d * -1
                    idx = np.argsort(x_d)
                    x_d, y_d = x_d[idx], y_d[idx]
                    
                    if self.check_promedio.isChecked():
                        if x_ref is None: x_ref, all_y = x_d, [y_d]
                        else: all_y.append(np.interp(x_ref, x_d, y_d))
                        self.ax_main.plot(x_d, y_d, color=c, alpha=0.1)
                    else: 
                        self.ax_main.plot(x_d, y_d, alpha=0.7, label=f"{os.path.basename(path)} ({y_t})")
    
                if self.check_promedio.isChecked() and all_y:
                    ym = np.mean(all_y, axis=0)
                    self.ax_main.plot(x_ref, ym, label=f"Media: {y_t}", lw=2, color=c)
                    res_medias[y_t] = ym
                    if i == 0: self.datos_promedio_actual = {'x': x_ref, 'y': ym, 'name': y_t}
    
            # 4. Cálculo y dibujo de la DIFERENCIA (¡Restaurado!)
            if self.ax_diff is not None and len(res_medias) >= 2:
                claves = list(res_medias.keys())
                # Restamos la segunda columna seleccionada de la primera
                diff = res_medias[claves[0]] - res_medias[claves[1]]
                self.ax_diff.plot(x_ref, diff, color='black', lw=1.5, label=f"Diff: {claves[0]}-{claves[1]}")
                self.ax_diff.axhline(0, color='red', linestyle='--', alpha=0.5)
                self.ax_diff.set_ylabel("Diferencia")
                self.ax_diff.legend(fontsize='small')
    
            self.ax_main.legend(fontsize='small')
            self.ax_main.set_ylabel("Intensidad")
            self.fig.tight_layout()
            self.canvas.draw()
            
    def autodetectar_primera_capa(self):
        # 1. Comprobar si hemos calculado FEFF antes
        if not hasattr(self, 'feff_folder') or not self.feff_folder:
            QMessageBox.warning(self, "Aviso", "Primero debes ejecutar FEFF (Pestaña 5) para generar los paths.")
            return
            
        # 2. Leemos la ventana de distancia de la interfaz
        rmin = self.spin_rmin.value()
        rmax = self.spin_rmax.value()
        
        archivos_feff = glob.glob(os.path.join(self.feff_folder, "feff[0-9][0-9][0-9][0-9].dat"))
        paths_encontrados = 0
        
        # 3. Escaneamos los archivos
        for archivo in archivos_feff:
            with open(archivo, 'r') as f:
                for linea in f:
                    # Buscamos la línea de la cabecera que tiene la distancia "reff"
                    if "nleg," in linea and "reff," in linea:
                        partes = linea.split()
                        try:
                            reff = float(partes[2])
                            
                            # MAGIA: Si el path cae dentro de Rmin y Rmax, lo añadimos
                            if rmin <= reff <= rmax:
                                self.agregar_path_dinamico(archivo)
                                paths_encontrados += 1
                                
                        except (ValueError, IndexError):
                            pass
                            
                        # Una vez leída la cabecera, pasamos al siguiente archivo para ir súper rápido
                        break 
                        
        # 4. Avisar al usuario del resultado
        if paths_encontrados > 0:
            self.consola_fiteo.appendPlainText(f"\n⚡ Autodetección: Se han añadido {paths_encontrados} paths en el rango de {rmin} a {rmax} Å.")
        else:
            QMessageBox.information(self, "Resultado", f"No se encontraron paths en el rango de {rmin} a {rmax} Å.\nPrueba a subir el Rmax o a generar los paths de FEFF de nuevo.")        
    def ejecutar_feff(self):
            if not self.feff_folder or not self.feff_filename: return
    
            self.consola_feff.appendPlainText(f"\n--- INICIANDO FEFF EN: {self.feff_folder} ---")
            self.consola_feff.appendPlainText("Calculando... (La interfaz puede congelarse unos segundos, es normal).")
            QApplication.processEvents()
    
            try:
                runner = feffrunner(folder=self.feff_folder, feffinp=self.feff_filename)
                runner.run()
    
                # CORRECCIÓN: Buscamos paths.dat o files.dat
                ruta_paths = os.path.join(self.feff_folder, 'paths.dat')
                ruta_files = os.path.join(self.feff_folder, 'files.dat')
                
                if os.path.exists(ruta_paths) or os.path.exists(ruta_files):
                    self.consola_feff.appendPlainText("\n--- FEFF TERMINADO CON ÉXITO ---")
                    
                    # ¡LLAMAMOS A LA FUNCIÓN QUE RELLENA LA TABLA!
                    self.leer_y_mostrar_paths_feff()
                    
                    num_paths = self.tabla_paths.rowCount()
                    QMessageBox.information(self, "Éxito", f"Cálculo completado.\nSe han cargado {num_paths} paths en la tabla.")
                else:
                    self.consola_feff.appendPlainText("\nERROR: FEFF terminó pero no generó 'paths.dat'. Revisa la consola para errores.")
    
            except Exception as e:
                self.consola_feff.appendPlainText(f"\nERROR FATAL AL EJECUTAR FEFF: {str(e)}")
                QMessageBox.critical(self, "Error", f"Fallo en la ejecución:\n{str(e)}")

if __name__ == "__main__":
    # Si ya hay una app corriendo (como en Spyder), la usamos. Si no, creamos una nueva.
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    app.setStyle("Fusion") 
    
    ventana = AnalizadorTotal()
    ventana.show()
    sys.exit(app.exec_())