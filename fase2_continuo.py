import numpy as np
import time

class ContinuousSimulation:
    def __init__(self, num_particles, box_size, 
                 mass=1.0, inertia=1.0, self_propulsion=1.0, 
                 linear_drag=1.0, rotational_drag=1.0, 
                 d_rep=1.0, k_rep=20.0, 
                 d_att=3.0, k_att=2.0, 
                 d_align=2.0, k_align=5.0, 
                 noise_amplitude=0.5, delta_t=0.02):
        """
        Inicializa la simulación del Modelo Continuo de Materia Activa (Fase 2).
        
        Args:
            num_particles (int): N, número de partículas.
            box_size (float): L, tamaño del lado del espacio cuadrado (límites periódicos).
            mass (float): m, masa de la partícula (resistencia a cambiar velocidad lineal).
            inertia (float): I, momento de inercia (resistencia a cambiar dirección/giro).
            self_propulsion (float): f_0, fuerza constante de autopropulsión hacia adelante.
            linear_drag (float): gamma_lin, coeficiente de fricción lineal.
            rotational_drag (float): gamma_rot, coeficiente de fricción rotacional.
            d_rep (float): Radio de repulsión de corto alcance.
            k_rep (float): Constante elástica de la fuerza de repulsión.
            d_att (float): Radio de atracción de mediano alcance.
            k_att (float): Constante elástica de la fuerza de atracción.
            d_align (float): Radio de vecindario para el torque de alineación.
            k_align (float): Fuerza del torque de alineación.
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
        self.d_rep = d_rep
        self.k_rep = k_rep
        self.d_att = d_att
        self.k_att = k_att
        self.d_align = d_align
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
        # Se asume una velocidad inicial típica en equilibrio con la propulsión y el arrastre (v = f0 / gamma_lin)
        v_init = self.f0 / self.gamma_lin if self.gamma_lin > 0 else 1.0
        self.velocities = np.zeros((self.N, 2))
        self.velocities[:, 0] = v_init * np.cos(self.angles)
        self.velocities[:, 1] = v_init * np.sin(self.angles)
        
    def step(self):
        """Avanza el estado físico de la simulación un paso delta_t usando Euler-Maruyama."""
        # 1. Cálculo de diferencias de posición con convención de mínima imagen (límites periódicos)
        # diff[i, j] representa el vector r_j - r_i
        diff = self.positions[np.newaxis, :, :] - self.positions[:, np.newaxis, :]
        diff = (diff + self.L / 2) % self.L - self.L / 2
        
        # Distancias euclidianas entre pares
        dist = np.sqrt(diff[:, :, 0]**2 + diff[:, :, 1]**2)
        
        # Para evitar división por cero en la diagonal (i == j)
        np.fill_diagonal(dist, 1.0)
        
        # Mascaras de vecindad (excluyendo la auto-interacción en la diagonal)
        eye_mask = ~np.eye(self.N, dtype=bool)
        rep_mask = (dist < self.d_rep) & eye_mask
        att_mask = (dist >= self.d_rep) & (dist < self.d_att) & eye_mask
        align_mask = (dist < self.d_align)  # Incluimos al propio elemento
        
        # 2. Fuerza de Repulsión de corto alcance
        # F_rep_ij = k_rep * (d_rep - r_ij) * (r_i - r_j) / r_ij
        # Como diff[i, j] = r_j - r_i, el vector unitario (r_i - r_j)/r_ij es -diff[i, j] / dist[i, j]
        rep_mag = np.zeros_like(dist)
        rep_mag[rep_mask] = -self.k_rep * (self.d_rep - dist[rep_mask]) / dist[rep_mask]
        F_rep = np.einsum('ij,ijo->io', rep_mag, diff)
        
        # 3. Fuerza de Atracción de alcance medio
        # F_att_ij = k_att * (r_ij - d_rep) * (r_j - r_i) / r_ij
        # El vector unitario (r_j - r_i)/r_ij es diff[i, j] / dist[i, j]
        att_mag = np.zeros_like(dist)
        att_mag[att_mask] = self.k_att * (dist[att_mask] - self.d_rep) / dist[att_mask]
        F_att = np.einsum('ij,ijo->io', att_mag, diff)
        
        # 4. Torque de Alineación
        # Promediamos la ORIENTACIÓN (theta) de los vecinos, no su velocidad instantánea.
        # La velocidad puede estar temporalmente distorsionada por colisiones (repulsión),
        # pero la orientación theta representa la "intención de movimiento" estable de cada partícula.
        # Usamos la técnica vectorial (cos/sin) para promediar ángulos correctamente,
        # igual que en el modelo de Vicsek (Fase 1).
        cos_angles = np.cos(self.angles)
        sin_angles = np.sin(self.angles)
        sum_cos = align_mask.dot(cos_angles)
        sum_sin = align_mask.dot(sin_angles)
        
        theta_avg = np.arctan2(sum_sin, sum_cos)
        
        # Torque: tau = k_align * sin(theta_avg - theta)
        tau_align = self.k_align * np.sin(theta_avg - self.angles)
        
        # 5. Fuerza de Autopropulsión lineal
        # F_self = f0 * (cos(theta), sin(theta))
        F_self = np.zeros((self.N, 2))
        F_self[:, 0] = self.f0 * np.cos(self.angles)
        F_self[:, 1] = self.f0 * np.sin(self.angles)
        
        # --- INTEGRACIÓN DE EDOs (EULER-MARUYAMA) ---
        
        # a. Actualización de Posición
        self.positions += self.velocities * self.dt
        self.positions %= self.L  # Límites periódicos
        
        # b. Actualización de Velocidad Lineal
        # dv = (dt / m) * (F_self + F_rep + F_att - gamma_lin * v)
        self.velocities += (self.dt / self.m) * (F_self + F_rep + F_att - self.gamma_lin * self.velocities)
        
        # c. Actualización del Ángulo
        self.angles += self.omega * self.dt
        # Mantenemos los ángulos en el intervalo [-pi, pi)
        self.angles = (self.angles + np.pi) % (2 * np.pi) - np.pi
        
        # d. Actualización de Velocidad Angular con Ruido
        # domega = (dt / I) * (tau_align - gamma_rot * omega) + eta * Uniform(-0.5, 0.5) * sqrt(dt)
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
    print("Iniciando prueba del Modelo Continuo de Materia Activa (Fase 2)...")
    
    # Parámetros para una prueba pequeña
    N_test = 5
    L_test = 10.0
    
    sim = ContinuousSimulation(
        num_particles=N_test, 
        box_size=L_test,
        mass=1.0, inertia=1.0,
        self_propulsion=1.0, linear_drag=1.0, rotational_drag=1.0,
        d_rep=1.0, k_rep=20.0,
        d_att=3.0, k_att=2.0,
        d_align=2.0, k_align=5.0,
        noise_amplitude=0.2, delta_t=0.02
    )
    
    print("\n[ESTADO INICIAL]")
    print(f"Posiciones:\n{np.round(sim.positions, 2)}")
    print(f"Velocidades:\n{np.round(sim.velocities, 2)}")
    print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
    print(f"Velocidad Angular (rad/s):\n{np.round(sim.omega, 2)}")
    
    start_time = time.time()
    
    # Avanzamos la simulación 3 pasos
    for step in range(1, 4):
        sim.step()
        print(f"\n--- Paso {step} ---")
        print(f"Posiciones:\n{np.round(sim.positions, 2)}")
        print(f"Velocidades:\n{np.round(sim.velocities, 2)}")
        print(f"Ángulos (rad):\n{np.round(sim.angles, 2)}")
        print(f"Velocidad Angular (rad/s):\n{np.round(sim.omega, 2)}")
        
    elapsed = time.time() - start_time
    print(f"\nSimulación de 3 pasos completada en {elapsed:.5f} segundos.")
