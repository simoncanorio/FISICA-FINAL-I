"""
Fase 3: Recolección de Datos — Barrido de Parámetros
=====================================================
Realiza un barrido automático sobre el nivel de ruido (eta) para ambos
modelos (Vicsek discreto y Materia Activa Continua), calculando en cada
paso el Parámetro de Orden Polar (Phi) y exportando los resultados a CSV.

Archivos generados:
    - resultados_vicsek.csv     : promedio de Phi vs eta (Modelo de Vicsek)
    - resultados_continuo.csv   : promedio de Phi vs eta (Modelo Continuo)
    - snapshots_vicsek.csv      : posiciones/ángulos al estabilizarse
    - snapshots_continuo.csv    : posiciones/velocidades/ángulos al estabilizarse
"""

import numpy as np
import csv
import time
import sys
from pathlib import Path

# Importamos los motores de simulación de las fases anteriores
from fase1_vicsek import VicsekSimulation
from fase2_continuo import ContinuousSimulation


# ===========================================================================
# FUNCIÓN DE CÁLCULO DEL PARÁMETRO DE ORDEN POLAR (Phi)
# ===========================================================================

def calcular_phi_vicsek(angles: np.ndarray) -> float:
    """
    Calcula el Parámetro de Orden Polar para el modelo de Vicsek (Fase 1).

    En el modelo discreto, todas las partículas tienen la misma rapidez v,
    por lo que el orden se mide directamente sobre los ángulos de orientación.

    Phi = (1/N) * |sum_i (cos(theta_i), sin(theta_i))|

    Valores:
        Phi ~ 1.0  →  alineación perfecta (enjambre ordenado)
        Phi ~ 0.0  →  movimiento aleatorio (caos total)

    Args:
        angles: Array de ángulos de orientación (rad) de las N partículas.

    Returns:
        Valor de Phi en el rango [0, 1].
    """
    N = len(angles)
    vx_sum = np.sum(np.cos(angles))
    vy_sum = np.sum(np.sin(angles))
    return np.sqrt(vx_sum**2 + vy_sum**2) / N


def calcular_phi_continuo(velocities: np.ndarray) -> float:
    """
    Calcula el Parámetro de Orden Polar para el modelo continuo (Fase 2).

    En el modelo continuo, la rapidez de cada partícula varía en el tiempo,
    por lo que normalizamos cada vector de velocidad antes de promediar.
    Esto evita que partículas rápidas (post-colisión) dominen el cálculo.

    Phi = (1/N) * |sum_i (v_i / |v_i|)|

    Args:
        velocities: Array (N, 2) de vectores de velocidad (vx, vy).

    Returns:
        Valor de Phi en el rango [0, 1].
    """
    N = len(velocities)
    speeds = np.linalg.norm(velocities, axis=1, keepdims=True)
    # Evitamos división por cero para partículas en reposo momentáneo
    speeds = np.where(speeds < 1e-9, 1.0, speeds)
    unit_vectors = velocities / speeds
    phi = np.linalg.norm(unit_vectors.sum(axis=0)) / N
    return float(phi)


# ===========================================================================
# BARRIDO DE PARÁMETROS — MODELO DE VICSEK (Fase 1)
# ===========================================================================

def barrer_vicsek(
    eta_values: np.ndarray,
    N: int = 150,
    L: float = 20.0,
    R: float = 2.0,
    v: float = 0.5,
    pasos_termal: int = 500,
    pasos_medicion: int = 300,
    archivo_resumen: str = "resultados_vicsek.csv",
    archivo_snapshots: str = "snapshots_vicsek.csv",
):
    """
    Barrido automático del parámetro de ruido eta para el Modelo de Vicsek.

    Para cada valor de eta, el sistema:
      1. Se inicializa aleatoriamente.
      2. Corre 'pasos_termal' pasos para alcanzar el estado estacionario.
      3. Mide Phi durante 'pasos_medicion' pasos y guarda el promedio.
      4. Exporta un snapshot del estado final.

    Args:
        eta_values      : Array de valores de ruido a explorar (de mayor a menor).
        N               : Número de partículas.
        L               : Tamaño de la caja.
        R               : Radio de interacción.
        v               : Velocidad escalar constante.
        pasos_termal    : Pasos de termalización (descartados del promedio).
        pasos_medicion  : Pasos de medición (promediados para obtener <Phi>).
        archivo_resumen : Nombre del CSV con eta vs <Phi>.
        archivo_snapshots: Nombre del CSV con el estado final de cada corrida.
    """
    print(f"\n{'='*60}")
    print(f"  BARRIDO VICSEK | N={N} | L={L} | R={R} | v={v}")
    print(f"  Pasos de termalización : {pasos_termal}")
    print(f"  Pasos de medición      : {pasos_medicion}")
    print(f"  Valores de eta         : {len(eta_values)}")
    print(f"{'='*60}")

    with (
        open(archivo_resumen, "w", newline="", encoding="utf-8") as f_res,
        open(archivo_snapshots, "w", newline="", encoding="utf-8") as f_snap,
    ):
        writer_res = csv.writer(f_res)
        writer_snap = csv.writer(f_snap)

        # Cabeceras
        writer_res.writerow(["eta", "phi_promedio", "phi_std", "phi_min", "phi_max",
                              "pasos_termal", "pasos_medicion", "N", "L", "R", "v"])
        writer_snap.writerow(["eta", "paso", "particula_id", "x", "y", "theta"])

        for idx, eta in enumerate(eta_values):
            t_inicio = time.time()

            sim = VicsekSimulation(
                num_particles=N,
                box_size=L,
                interaction_radius=R,
                noise_amplitude=float(eta),
                speed=v,
            )

            # --- Termalización (descartamos datos transitorios) ---
            for _ in range(pasos_termal):
                sim.step()

            # --- Medición del estado estacionario ---
            phi_series = np.zeros(pasos_medicion)
            for t in range(pasos_medicion):
                sim.step()
                phi_series[t] = calcular_phi_vicsek(sim.angles)

            phi_avg = float(np.mean(phi_series))
            phi_std = float(np.std(phi_series))
            phi_min = float(np.min(phi_series))
            phi_max = float(np.max(phi_series))

            # Guardar resumen
            writer_res.writerow([
                round(float(eta), 6), round(phi_avg, 6), round(phi_std, 6),
                round(phi_min, 6), round(phi_max, 6),
                pasos_termal, pasos_medicion, N, L, R, v,
            ])

            # Guardar snapshot del estado final (último paso de medición)
            paso_total = pasos_termal + pasos_medicion
            for pid in range(N):
                writer_snap.writerow([
                    round(float(eta), 6),
                    paso_total,
                    pid,
                    round(sim.positions[pid, 0], 6),
                    round(sim.positions[pid, 1], 6),
                    round(sim.angles[pid], 6),
                ])

            elapsed = time.time() - t_inicio
            bar = "#" * int(30 * (idx + 1) / len(eta_values))
            bar = bar.ljust(30)
            print(f"  [{bar}] eta={eta:.3f}  <Phi>={phi_avg:.4f} ± {phi_std:.4f}  ({elapsed:.1f}s)")

    print(f"\n  ✓ Resumen guardado en  : {archivo_resumen}")
    print(f"  ✓ Snapshots guardados en: {archivo_snapshots}")


# ===========================================================================
# BARRIDO DE PARÁMETROS — MODELO CONTINUO (Fase 2)
# ===========================================================================

def barrer_continuo(
    eta_values: np.ndarray,
    N: int = 150,
    L: float = 20.0,
    mass: float = 1.0,
    inertia: float = 0.3,
    self_prop: float = 0.5,
    lin_drag: float = 1.0,
    rot_drag: float = 3.0,
    d_rep: float = 0.4,
    k_rep: float = 40.0,
    d_att: float = 3.0,
    k_att: float = 0.15,
    d_align: float = 4.0,
    k_align: float = 8.0,
    dt: float = 0.03,
    pasos_termal: int = 3000,
    pasos_medicion: int = 1500,
    archivo_resumen: str = "resultados_continuo.csv",
    archivo_snapshots: str = "snapshots_continuo.csv",
):
    """
    Barrido automático del parámetro de ruido eta para el Modelo Continuo.

    Funciona igual que 'barrer_vicsek' pero usa ContinuousSimulation
    y registra también las velocidades instantáneas en el snapshot.

    Args:
        eta_values      : Array de valores de ruido a explorar (de mayor a menor).
        (resto)         : Parámetros físicos del modelo continuo (ver fase2_continuo.py).
        pasos_termal    : Pasos de termalización. Se recomienda un valor mayor que
                          Vicsek porque el dt del modelo continuo es mucho menor.
        pasos_medicion  : Pasos de medición promediados.
        archivo_resumen : Nombre del CSV con eta vs <Phi>.
        archivo_snapshots: Nombre del CSV con el estado final de cada corrida.
    """
    print(f"\n{'='*60}")
    print(f"  BARRIDO CONTINUO | N={N} | L={L} | dt={dt}")
    print(f"  Pasos de termalización : {pasos_termal}")
    print(f"  Pasos de medición      : {pasos_medicion}")
    print(f"  Valores de eta         : {len(eta_values)}")
    print(f"{'='*60}")

    with (
        open(archivo_resumen, "w", newline="", encoding="utf-8") as f_res,
        open(archivo_snapshots, "w", newline="", encoding="utf-8") as f_snap,
    ):
        writer_res = csv.writer(f_res)
        writer_snap = csv.writer(f_snap)

        # Cabeceras
        writer_res.writerow(["eta", "phi_promedio", "phi_std", "phi_min", "phi_max",
                              "pasos_termal", "pasos_medicion", "N", "L", "dt"])
        writer_snap.writerow(["eta", "paso", "particula_id",
                               "x", "y", "vx", "vy", "theta", "omega"])

        for idx, eta in enumerate(eta_values):
            t_inicio = time.time()

            sim = ContinuousSimulation(
                num_particles=N, box_size=L,
                mass=mass, inertia=inertia,
                self_propulsion=self_prop, linear_drag=lin_drag, rotational_drag=rot_drag,
                d_rep=d_rep, k_rep=k_rep,
                d_att=d_att, k_att=k_att,
                d_align=d_align, k_align=k_align,
                noise_amplitude=float(eta), delta_t=dt,
            )

            # --- Termalización ---
            for _ in range(pasos_termal):
                sim.step()

            # --- Medición ---
            phi_series = np.zeros(pasos_medicion)
            for t in range(pasos_medicion):
                sim.step()
                phi_series[t] = calcular_phi_continuo(sim.velocities)

            phi_avg = float(np.mean(phi_series))
            phi_std = float(np.std(phi_series))
            phi_min = float(np.min(phi_series))
            phi_max = float(np.max(phi_series))

            # Guardar resumen
            writer_res.writerow([
                round(float(eta), 6), round(phi_avg, 6), round(phi_std, 6),
                round(phi_min, 6), round(phi_max, 6),
                pasos_termal, pasos_medicion, N, L, dt,
            ])

            # Snapshot del estado final
            paso_total = pasos_termal + pasos_medicion
            for pid in range(N):
                writer_snap.writerow([
                    round(float(eta), 6),
                    paso_total,
                    pid,
                    round(sim.positions[pid, 0], 6),
                    round(sim.positions[pid, 1], 6),
                    round(sim.velocities[pid, 0], 6),
                    round(sim.velocities[pid, 1], 6),
                    round(sim.angles[pid], 6),
                    round(sim.omega[pid], 6),
                ])

            elapsed = time.time() - t_inicio
            bar = "#" * int(30 * (idx + 1) / len(eta_values))
            bar = bar.ljust(30)
            print(f"  [{bar}] eta={eta:.3f}  <Phi>={phi_avg:.4f} ± {phi_std:.4f}  ({elapsed:.1f}s)")

    print(f"\n  ✓ Resumen guardado en  : {archivo_resumen}")
    print(f"  ✓ Snapshots guardados en: {archivo_snapshots}")


# ===========================================================================
# PUNTO DE ENTRADA PRINCIPAL
# ===========================================================================

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Rango de ruido: de alto (caos) a bajo (orden), con 20 puntos.
    # Puedes modificar estos valores fácilmente aquí.
    # ------------------------------------------------------------------
    ETA_MAX   = 2.0    # Ruido máximo (movimiento casi aleatorio)
    ETA_MIN   = 0.05   # Ruido mínimo (alineación casi perfecta)
    N_PUNTOS  = 20     # Número de valores de eta a explorar

    eta_values = np.linspace(ETA_MAX, ETA_MIN, N_PUNTOS)

    print("\n╔══════════════════════════════════════════════════════════╗")
    print("║        FASE 3 — RECOLECCIÓN MASIVA DE DATOS             ║")
    print("║  Barrido de eta: de {:.2f} a {:.2f} ({} puntos)    ║".format(ETA_MAX, ETA_MIN, N_PUNTOS))
    print("╚══════════════════════════════════════════════════════════╝")

    t0_global = time.time()

    # --- Barrido del Modelo de Vicsek (Fase 1) ---
    barrer_vicsek(
        eta_values=eta_values,
        N=150, L=20.0, R=2.0, v=0.5,
        pasos_termal=500,
        pasos_medicion=300,
        archivo_resumen="resultados_vicsek.csv",
        archivo_snapshots="snapshots_vicsek.csv",
    )

    # --- Barrido del Modelo Continuo (Fase 2) ---
    barrer_continuo(
        eta_values=eta_values,
        N=150, L=20.0,
        pasos_termal=3000,
        pasos_medicion=1500,
        archivo_resumen="resultados_continuo.csv",
        archivo_snapshots="snapshots_continuo.csv",
    )

    total = time.time() - t0_global
    print(f"\n✓ Barrido completo finalizado en {total:.1f} segundos.")
    print("  Archivos generados:")
    print("    → resultados_vicsek.csv")
    print("    → resultados_continuo.csv")
    print("    → snapshots_vicsek.csv")
    print("    → snapshots_continuo.csv")
