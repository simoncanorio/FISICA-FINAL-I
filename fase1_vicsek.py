import numpy as np
from scipy.spatial import cKDTree
import time

class VicsekSimulation:
    def __init__(self, num_particles, box_size, interaction_radius, noise_amplitude, speed, delta_t=1.0):
        """
        Inicializa la simulación del Modelo Discreto de Vicsek.
        
        Args:
            num_particles (int): N, número de partículas.
            box_size (float): L, tamaño del lado del espacio cuadrado (topología de toroide).
            interaction_radius (float): R, radio de vecindad para la alineación.
            noise_amplitude (float): eta, amplitud máxima del ruido angular en radianes ([-eta/2, eta/2]).
            speed (float): v, velocidad escalar constante de las partículas.
            delta_t (float): Paso de tiempo para la actualización cinemática.
        """
        self.N = num_particles
        self.L = box_size
        self.R = interaction_radius
        self.eta = noise_amplitude
        self.v = speed
        self.dt = delta_t
        
        # Inicialización de posiciones (x, y) uniformemente distribuidas en [0, L)
        self.positions = np.random.uniform(0, self.L, size=(self.N, 2))
        
        # Inicialización de ángulos aleatorios en [-pi, pi)
        self.angles = np.random.uniform(-np.pi, np.pi, size=self.N)
        
    def step(self):
        """Avanza la simulación un paso discreto de tiempo (dt)."""
        # 1. Detección del Vecindario (Optimizado espacialmente)
        # cKDTree con el parámetro 'boxsize' implementa de forma nativa la métrica de 
        # distancia sobre una topología de toroide (límites periódicos).
        tree = cKDTree(self.positions, boxsize=self.L)
        
        # Buscamos todos los vecinos (incluyéndose a sí mismos) dentro del radio R
        # Retorna una lista donde el i-ésimo elemento es una lista de los índices de sus vecinos.
        neighbors_list = tree.query_ball_tree(tree, r=self.R)
        
        # 2. Regla de Alineación
        new_angles = np.zeros(self.N)
        
        # Para promediar ángulos correctamente y evitar problemas con el salto de pi a -pi,
        # promediamos los vectores de dirección (coseno y seno) en lugar de los ángulos directamente.
        cos_angles = np.cos(self.angles)
        sin_angles = np.sin(self.angles)
        
        for i, neighbors in enumerate(neighbors_list):
            # Suma de las componentes de velocidad de todos los vecinos
            sum_cos = np.sum(cos_angles[neighbors])
            sum_sin = np.sum(sin_angles[neighbors])
            
            # Ángulo promedio resultante (arctan2 determina el cuadrante correcto)
            avg_angle = np.arctan2(sum_sin, sum_cos)
            
            # Añadimos el ruido uniforme
            noise = np.random.uniform(-self.eta / 2, self.eta / 2)
            new_angles[i] = avg_angle + noise
            
        self.angles = new_angles
        
        # 3. Actualización Cinemática
        # Calculamos los desplazamientos en x e y
        vx = self.v * np.cos(self.angles)
        vy = self.v * np.sin(self.angles)
        
        self.positions[:, 0] += vx * self.dt
        self.positions[:, 1] += vy * self.dt
        
        # Aplicamos límites periódicos asegurando que x e y se mantengan en el rango [0, L).
        # El operador módulo (%) en Python maneja correctamente los números negativos.
        self.positions %= self.L
        
    def get_state(self):
        """
        Devuelve el estado actual de la simulación. 
        Este formato (listas/diccionario) es fácilmente serializable a JSON 
        para enviarlo a un dashboard interactivo web o API en el futuro.
        """
        return {
            'positions': self.positions.tolist(),
            'angles': self.angles.tolist()
        }

if __name__ == "__main__":
    # --- PRUEBA EN TERMINAL ---
    print("Iniciando prueba del Modelo de Vicsek (Fase 1)...")
    
    # Parámetros para una prueba pequeña
    N_test = 5
    L_test = 10.0
    
    sim = VicsekSimulation(
        num_particles=N_test, 
        box_size=L_test, 
        interaction_radius=3.0, 
        noise_amplitude=0.5, 
        speed=1.0
    )
    
    print("\n[ESTADO INICIAL]")
    print(f"Posiciones:\n{np.round(sim.positions, 2)}")
    print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
    
    start_time = time.time()
    
    # Avanzamos la simulación 3 pasos
    for step in range(1, 4):
        sim.step()
        print(f"\n--- Paso {step} ---")
        print(f"Posiciones:\n{np.round(sim.positions, 2)}")
        print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
        
    elapsed = time.time() - start_time
    print(f"\nSimulación de 3 pasos completada en {elapsed:.5f} segundos.")
