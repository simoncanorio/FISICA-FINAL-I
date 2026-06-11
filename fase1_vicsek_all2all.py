import numpy as np
import time

class VicsekSimulationAll2All:
    def __init__(self, num_particles, box_size, interaction_radius, noise_amplitude, speed, delta_t=1):
        """
        Inicializa la simulación del Modelo Discreto de Vicsek usando matriz de distancias (Todos contra Todos).
        
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
        """Avanza la simulación un paso discreto de tiempo (dt) comparando todos con todos."""
        # 1. Todos interactúan con todos sin importar la distancia (Campo Medio global)
        cos_angles = np.cos(self.angles)
        sin_angles = np.sin(self.angles)
        
        # Sumamos los vectores de dirección de TODAS las partículas
        sum_cos = np.sum(cos_angles)
        sum_sin = np.sum(sin_angles)
        
        # Ángulo promedio global resultante
        avg_angles = np.arctan2(sum_sin, sum_cos)
        
        # Añadimos el ruido uniforme
        noise = np.random.uniform(-self.eta / 2, self.eta / 2, size=self.N)
        self.angles = avg_angles + noise
        
        # 4. Actualización Cinemática
        self.positions[:, 0] += self.v * np.cos(self.angles) * self.dt
        self.positions[:, 1] += self.v * np.sin(self.angles) * self.dt
        
        # Aplicamos límites periódicos asegurando que x e y se mantengan en el rango [0, L)
        self.positions %= self.L
        
    def get_state(self):
        return {
            'positions': self.positions.tolist(),
            'angles': self.angles.tolist()
        }

if __name__ == "__main__":
    # --- PRUEBA EN TERMINAL ---
    print("Iniciando prueba del Modelo de Vicsek (All-To-All)...")
    
    # Parámetros para una prueba pequeña
    N_test = 500
    L_test = 20.0
    
    sim = VicsekSimulationAll2All(
        num_particles=N_test, 
        box_size=L_test, 
        interaction_radius=2.0, 
        noise_amplitude=0.5, 
        speed=1.0
    )
    
    start_time = time.time()
    
    # Avanzamos la simulación 100 pasos para ver el rendimiento
    pasos = 100
    for step in range(pasos):
        sim.step()
        
    elapsed = time.time() - start_time
    print(f"\nSimulación de {pasos} pasos con {N_test} partículas completada en {elapsed:.5f} segundos.")
