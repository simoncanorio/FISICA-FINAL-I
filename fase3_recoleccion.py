import numpy as np
import pandas as pd
import os

def calcular_parametro_orden(velocidades, velocidad_escalar):
    """
    Calcula el parámetro de orden polar (Phi) de forma eficiente y vectorizada.
    
    Args:
        velocidades (np.ndarray): Arreglo de forma (N, d) con los vectores de velocidad 
                                  de cada partícula (d es la dimensión, ej. 2 o 3).
        velocidad_escalar (float): Módulo constante de la velocidad (v0).
        
    Returns:
        float: Valor de Phi, normalizado entre [0, 1].
    """
    N = velocidades.shape[0]
    
    # 1. Suma vectorial a lo largo del eje de las partículas (eje 0)
    # Resultado: Un vector [suma_vx, suma_vy]
    suma_velocidades = np.sum(velocidades, axis=0)
    
    # 2. Magnitud del vector suma
    magnitud_resultante = np.linalg.norm(suma_velocidades)
    
    # 3. Normalización por número de partículas y su velocidad intrínseca
    phi = magnitud_resultante / (N * velocidad_escalar)
    
    return phi

def ejecutar_barrido(
    simulador, 
    n_particulas, 
    velocidad_escalar,
    ruido_max=5.0, 
    ruido_min=0.0, 
    pasos_ruido=20,
    t_termalizacion=1000, 
    t_muestreo=500,
    guardar_detalles=False,
    dir_salida="resultados_fase3"
):
    """
    Ejecuta el experimento automatizado realizando un barrido descendente de ruido.
    
    Args:
        simulador: Objeto de simulación (debe tener métodos `step()`, `get_velocities()`, `get_positions()`, y atributo `eta`).
        n_particulas (int): Número de agentes en el sistema (N).
        velocidad_escalar (float): Velocidad constante v0.
        ruido_max (float): Nivel máximo de ruido (inicio).
        ruido_min (float): Nivel mínimo de ruido (fin).
        pasos_ruido (int): Cantidad de niveles de ruido a evaluar.
        t_termalizacion (int): Pasos para alcanzar estado estacionario antes de medir.
        t_muestreo (int): Pasos durante los cuales se tomarán medidas para promediar.
        guardar_detalles (bool): Si es True, guarda posiciones/velocidades por cada paso.
        dir_salida (str): Carpeta destino para los archivos CSV.
    """
    os.makedirs(dir_salida, exist_ok=True)
    
    # Generamos los niveles de ruido de mayor a menor para simular un "enfriamiento"
    # y capturar correctamente la histéresis si la transición de fase es de primer orden.
    niveles_ruido = np.linspace(ruido_max, ruido_min, pasos_ruido)
    
    resultados_resumen = []
    
    print(f"=== Iniciando Barrido de Parámetros ===")
    print(f"Ruido: {ruido_max} -> {ruido_min} ({pasos_ruido} niveles)")
    print(f"Directorio de salida: {dir_salida}/")
    
    # Inicializamos el sistema con ruido máximo para destruir cualquier orden inicial
    simulador.eta = ruido_max
    
    for idx, eta in enumerate(niveles_ruido):
        print(f"[{idx+1}/{pasos_ruido}] Evaluando η = {eta:.3f} ... ", end="")
        simulador.eta = eta
        
        # --- Fase de Termalización ---
        # Dejamos evolucionar el sistema para que se adapte al nuevo nivel de ruido
        for _ in range(t_termalizacion):
            simulador.step()
            
        # --- Fase de Muestreo ---
        phis_temporales = np.zeros(t_muestreo)
        detalles_pasos = [] if guardar_detalles else None
        
        for t in range(t_muestreo):
            simulador.step()
            velocidades = simulador.get_velocities()
            
            # Cálculo vectorizado
            phi_t = calcular_parametro_orden(velocidades, velocidad_escalar)
            phis_temporales[t] = phi_t
            
            if guardar_detalles:
                posiciones = simulador.get_positions()
                # Construimos un dataframe rápido para el frame actual
                frame_data = np.column_stack((
                    np.full(n_particulas, t),      # timestep
                    np.arange(n_particulas),       # particula_id
                    posiciones,                    # x, y
                    velocidades                    # vx, vy
                ))
                detalles_pasos.append(frame_data)
                
        # --- Estadísticas del Nivel de Ruido ---
        phi_promedio = np.mean(phis_temporales)
        phi_varianza = np.var(phis_temporales)
        
        resultados_resumen.append({
            'ruido_eta': eta,
            'phi_promedio': phi_promedio,
            'phi_varianza': phi_varianza
        })
        print(f"Φ promedio = {phi_promedio:.4f}")
        
        # --- Exportación de Detalles (Opcional) ---
        if guardar_detalles:
            # Concatenamos la lista de arrays de 2D en un solo array de 2D
            datos_completos = np.vstack(detalles_pasos)
            columnas = ['timestep', 'particula_id', 'x', 'y', 'vx', 'vy']
            df_detalles = pd.DataFrame(datos_completos, columns=columnas)
            
            # Optimizamos tipos de datos para no hacer el CSV tan pesado
            df_detalles = df_detalles.astype({
                'timestep': 'int32', 
                'particula_id': 'int32'
            })
            
            archivo_detalles = os.path.join(dir_salida, f"trayectorias_eta_{eta:.3f}.csv")
            df_detalles.to_csv(archivo_detalles, index=False)

    # --- Exportación del Resumen Final ---
    df_resumen = pd.DataFrame(resultados_resumen)
    archivo_resumen = os.path.join(dir_salida, "resumen_transicion_fase.csv")
    df_resumen.to_csv(archivo_resumen, index=False)
    
    print("=== Barrido Completado ===")
    print(f"Resumen guardado en: {archivo_resumen}")
    return df_resumen

# ==========================================
# Ejemplo de uso / Adapter de Integración
# ==========================================
if __name__ == "__main__":
    # Importamos el simulador del main.py que el usuario ya tiene
    try:
        from main import VicsekSimulation
        
        # Creamos una clase Wrapper (Adapter) para asegurar que el simulador
        # tenga los métodos exactos que requiere nuestro recolector de datos.
        class SimuladorAdapter:
            def __init__(self, N, v0, eta):
                self.sim = VicsekSimulation(N=N, v0=v0, eta=eta)
                self.v0 = v0
                
            @property
            def eta(self):
                return self.sim.eta
                
            @eta.setter
            def eta(self, valor):
                self.sim.eta = valor
                
            def step(self):
                self.sim.step()
                
            def get_velocities(self):
                # El modelo actual de Vicsek (Fase 1) usa ángulos (theta).
                # Convertimos ángulos a vectores velocidad (vx, vy).
                vx = np.cos(self.sim.theta) * self.v0
                vy = np.sin(self.sim.theta) * self.v0
                return np.column_stack((vx, vy))
                
            def get_positions(self):
                return self.sim.pos

        # Parámetros del experimento
        N_PARTICULAS = 250
        VELOCIDAD_V0 = 0.1
        
        mi_simulador = SimuladorAdapter(N=N_PARTICULAS, v0=VELOCIDAD_V0, eta=5.0)
        
        # Ejecutamos un barrido rápido para prueba
        ejecutar_barrido(
            simulador=mi_simulador,
            n_particulas=N_PARTICULAS,
            velocidad_escalar=VELOCIDAD_V0,
            ruido_max=5.0,
            ruido_min=0.0,
            pasos_ruido=10,        # Solo 10 puntos de ruido para prueba rápida
            t_termalizacion=200,   # Reducido para prueba rápida
            t_muestreo=100,
            guardar_detalles=False # Cambiar a True si quieres las trayectorias completas
        )
    except ImportError:
        print("No se encontró main.py para ejecutar la prueba de integración.")
