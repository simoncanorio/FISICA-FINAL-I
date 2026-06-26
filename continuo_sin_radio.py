import numpy as np
import time

class ContinuousSimulationSinRadio:
    def __init__(self, num_particles, box_size, 
                 mass=1.0, inertia=1.0, self_propulsion=1.0, 
                 linear_drag=1.0, rotational_drag=1.0, 
                 k_align=5.0, noise_amplitude=0.5, delta_t=0.02):
        """
        Inicializa la simulación del Modelo Continuo de Materia Activa sin radio de vecindario (Acoplamiento Global).
        
        Args:
            num_particles (int): N, número de partículas.
            box_size (float): L, tamaño del lado del espacio cuadrado (límites periódicos).
            mass (float): m, masa de la partícula (resistencia a cambiar velocidad lineal).
            inertia (float): I, momento de inercia (resistencia a cambiar dirección/giro).
            self_propulsion (float): f_0, fuerza constante de autopropulsión hacia adelante.
            linear_drag (float): gamma_lin, coeficiente de fricción lineal.
            rotational_drag (float): gamma_rot, coeficiente de fricción rotacional.
            k_align (float): Fuerza del torque de alineación global.
            noise_amplitude (float): eta, intensidad del ruido angular.
            delta_t (float): dt, incremento temporal minúsculo para la integración.
        """
        self.N = num_particles
        self.L = box_size
        self.m = mass
        self.I = inertia
        self.f0 = self_propulsion
        self.gamma_lin = linear_drag
        self.gamma_rot = rotational_drag
        self.k_align = k_align
        self.eta = noise_amplitude
        self.dt = delta_t
        
        # Inicialización de posiciones (x, y) uniformemente distribuidas en [0, L)
        self.positions = np.random.uniform(0, self.L, size=(self.N, 2))
        
        # Inicialización de ángulos aleatorios en [-pi, pi)
        self.angles = np.random.uniform(-np.pi, np.pi, size=self.N)
        
        # Inicialización de velocidad angular (omega) en cero
        self.omega = np.zeros(self.N)
        
        # Inicialización de velocidades lineales (v_x, v_y) alineadas con los ángulos iniciales
        v_init = self.f0 / self.gamma_lin if self.gamma_lin > 0 else 1.0
        self.velocities = np.zeros((self.N, 2))
        self.velocities[:, 0] = v_init * np.cos(self.angles)
        self.velocities[:, 1] = v_init * np.sin(self.angles)
        
    def step(self):
        """Avanza el estado físico de la simulación un paso delta_t sin buscar vecinos."""
        # 1. Fuerza de Autopropulsión lineal
        F_self = np.zeros((self.N, 2))
        F_self[:, 0] = self.f0 * np.cos(self.angles)
        F_self[:, 1] = self.f0 * np.sin(self.angles)
        
        # 2. Torque de Alineación Global (Campo Medio - Todos contra Todos)
        # En lugar de promediar vecinos en un radio, promediamos la dirección de TODAS las partículas del sistema.
        cos_angles = np.cos(self.angles)
        sin_angles = np.sin(self.angles)
        sum_cos = np.sum(cos_angles)
        sum_sin = np.sum(sin_angles)
        
        # Ángulo promedio global
        theta_avg = np.arctan2(sum_sin, sum_cos)
        
        # Torque global: tau = k_align * sin(theta_avg - theta)
        tau_align = self.k_align * np.sin(theta_avg - self.angles)
        
        # --- INTEGRACIÓN DE EDOs (EULER-MARUYAMA) ---
        
        # a. Actualización de Posición
        self.positions += self.velocities * self.dt
        self.positions %= self.L  # Límites periódicos
        
        # b. Actualización de Velocidad Lineal (sin fuerzas espaciales localizadas)
        # dv = (dt / m) * (F_self - gamma_lin * v)
        self.velocities += (self.dt / self.m) * (F_self - self.gamma_lin * self.velocities)
        
        # c. Actualización del Ángulo
        self.angles += self.omega * self.dt
        self.angles = (self.angles + np.pi) % (2 * np.pi) - np.pi
        
        # d. Actualización de Velocidad Angular con Ruido
        noise = self.eta * np.random.uniform(-0.5, 0.5, size=self.N) * np.sqrt(self.dt)
        self.omega += (self.dt / self.I) * (tau_align - self.gamma_rot * self.omega) + noise
        
    def get_state(self):
        """Devuelve el estado actual para compatibilidad y análisis."""
        return {
            'positions': self.positions.copy(),
            'velocities': self.velocities.copy(),
            'angles': self.angles.copy(),
            'omega': self.omega.copy()
        }

if __name__ == "__main__":
    # --- PRUEBA EN TERMINAL ---
    print("Iniciando prueba del Modelo Continuo Sin Radio (Acoplamiento Global)...")
    
    N_test = 5
    L_test = 10.0
    
    sim = ContinuousSimulationSinRadio(
        num_particles=N_test, 
        box_size=L_test,
        mass=1.0, inertia=1.0,
        self_propulsion=1.0, linear_drag=1.0, rotational_drag=1.0,
        k_align=5.0, noise_amplitude=0.2, delta_t=0.02
    )
    
    print("\n[ESTADO INICIAL]")
    print(f"Posiciones:\n{np.round(sim.positions, 2)}")
    print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
    
    for step in range(1, 4):
        sim.step()
        print(f"\n--- Paso {step} ---")
        print(f"Posiciones:\n{np.round(sim.positions, 2)}")
        print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
