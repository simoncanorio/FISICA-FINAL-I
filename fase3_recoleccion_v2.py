"""
Fase 3 (v2): Recolección de Datos — Barrido de Parámetros
==========================================================
Realiza un barrido automático sobre el nivel de ruido (eta) para los modelos
de acoplamiento global (campo medio):
  - Vicsek All-to-All    : todos interactúan con todos sin radio de corte.
  - Continuo Sin Radio   : modelo continuo con torque de alineación global.

Calcula el Parámetro de Orden Polar (Phi) y exporta resultados a CSV.

Archivos generados:
    - resultados_vicsek_all2all.csv       : promedio de Phi vs eta
    - resultados_continuo_sin_radio.csv   : promedio de Phi vs eta
    - snapshots_vicsek_all2all.csv        : posiciones/ángulos al estabilizarse
    - snapshots_continuo_sin_radio.csv    : posiciones/velocidades/ángulos al estabilizarse
"""

import numpy as np
import csv
import time
from pathlib import Path

# Importamos los motores de simulación de campo medio
from fase1_vicsek_all2all import VicsekSimulationAll2All
from continuo_sin_radio import ContinuousSimulationSinRadio


# ===========================================================================
# FUNCIÓN DE CÁLCULO DEL PARÁMETRO DE ORDEN POLAR (Phi)
# ===========================================================================

def calcular_phi_vicsek(angles: np.ndarray) -> float:
    """
    Calcula el Parámetro de Orden Polar para el modelo de Vicsek All-to-All.

    En este modelo todas las partículas tienen la misma rapidez v,
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
    Calcula el Parámetro de Orden Polar para el modelo continuo sin radio.

    La rapidez de cada partícula varía en el tiempo, por lo que normalizamos
    cada vector de velocidad antes de promediar. Esto evita que partículas
    rápidas dominen el cálculo.

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
# BARRIDO DE PARÁMETROS — MODELO DE VICSEK ALL-TO-ALL
# ===========================================================================

def barrer_vicsek_all2all(
    eta_values: np.ndarray,
    N: int = 150,
    L: float = 20.0,
    R: float = 2.0,
    v: float = 0.5,
    pasos_termal: int = 500,
    pasos_medicion: int = 300,
    archivo_resumen: str = "resultados_vicsek_all2all.csv",
    archivo_snapshots: str = "snapshots_vicsek_all2all.csv",
):
    """
    Barrido automático del parámetro de ruido eta para el Modelo de Vicsek All-to-All.

    En este modelo NO hay radio de interacción: cada partícula se alinea con
    el promedio global de TODAS las demás (campo medio / acoplamiento global).
    El parámetro R se conserva en la interfaz para compatibilidad pero no
    afecta la dinámica interna del modelo.

    Para cada valor de eta el sistema:
      1. Se inicializa aleatoriamente.
      2. Corre 'pasos_termal' pasos para alcanzar el estado estacionario.
      3. Mide Phi durante 'pasos_medicion' pasos y guarda el promedio.
      4. Exporta un snapshot del estado final.

    Args:
        eta_values      : Array de valores de ruido a explorar (de mayor a menor).
        N               : Número de partículas.
        L               : Tamaño de la caja.
        R               : Radio de interacción (solo informativo, no usado en dinámica).
        v               : Velocidad escalar constante.
        pasos_termal    : Pasos de termalización (descartados del promedio).
        pasos_medicion  : Pasos de medición (promediados para obtener <Phi>).
        archivo_resumen : Nombre del CSV con eta vs <Phi>.
        archivo_snapshots: Nombre del CSV con el estado final de cada corrida.
    """
    print(f"\n{'='*60}")
    print(f"  BARRIDO VICSEK ALL-TO-ALL | N={N} | L={L} | v={v}")
    print(f"  Modo            : Campo Medio Global (sin radio de corte)")
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
                              "pasos_termal", "pasos_medicion", "N", "L", "v"])
        writer_snap.writerow(["eta", "paso", "particula_id", "x", "y", "theta"])

        for idx, eta in enumerate(eta_values):
            t_inicio = time.time()

            sim = VicsekSimulationAll2All(
                num_particles=N,
                box_size=L,
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
                pasos_termal, pasos_medicion, N, L, v,
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

    print(f"\n  ✓ Resumen guardado en   : {archivo_resumen}")
    print(f"  ✓ Snapshots guardados en: {archivo_snapshots}")


# ===========================================================================
# BARRIDO DE PARÁMETROS — MODELO CONTINUO SIN RADIO
# ===========================================================================

def barrer_continuo_sin_radio(
    eta_values: np.ndarray,
    N: int = 150,
    L: float = 20.0,
    mass: float = 1.0,
    inertia: float = 1.0,
    self_prop: float = 1.0,
    lin_drag: float = 1.0,
    rot_drag: float = 1.0,
    k_align: float = 5.0,
    dt: float = 0.02,
    pasos_termal: int = 3000,
    pasos_medicion: int = 1500,
    archivo_resumen: str = "resultados_continuo_sin_radio.csv",
    archivo_snapshots: str = "snapshots_continuo_sin_radio.csv",
):
    """
    Barrido automático del parámetro de ruido eta para el Modelo Continuo Sin Radio.

    En este modelo el torque de alineación actúa sobre el promedio global de
    TODOS los ángulos del sistema (acoplamiento de campo medio continuo).
    No existe radio de vecindario local: todas las partículas se alinean entre sí.

    Funciona igual que 'barrer_vicsek_all2all' pero usa ContinuousSimulationSinRadio
    y registra también las velocidades y velocidad angular en el snapshot.

    Args:
        eta_values      : Array de valores de ruido a explorar (de mayor a menor).
        (resto)         : Parámetros físicos del modelo continuo (ver continuo_sin_radio.py).
        pasos_termal    : Pasos de termalización. Se recomienda mayor que Vicsek porque
                          el dt del modelo continuo es más pequeño.
        pasos_medicion  : Pasos de medición promediados.
        archivo_resumen : Nombre del CSV con eta vs <Phi>.
        archivo_snapshots: Nombre del CSV con el estado final de cada corrida.
    """
    print(f"\n{'='*60}")
    print(f"  BARRIDO CONTINUO SIN RADIO | N={N} | L={L} | dt={dt}")
    print(f"  Modo            : Acoplamiento Global (sin radio de vecindario)")
    print(f"  k_align         : {k_align}")
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
                              "pasos_termal", "pasos_medicion", "N", "L", "dt", "k_align"])
        writer_snap.writerow(["eta", "paso", "particula_id",
                               "x", "y", "vx", "vy", "theta", "omega"])

        for idx, eta in enumerate(eta_values):
            t_inicio = time.time()

            sim = ContinuousSimulationSinRadio(
                num_particles=N,
                box_size=L,
                mass=mass,
                inertia=inertia,
                self_propulsion=self_prop,
                linear_drag=lin_drag,
                rotational_drag=rot_drag,
                k_align=k_align,
                noise_amplitude=float(eta),
                delta_t=dt,
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
                pasos_termal, pasos_medicion, N, L, dt, k_align,
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

    print(f"\n  ✓ Resumen guardado en   : {archivo_resumen}")
    print(f"  ✓ Snapshots guardados en: {archivo_snapshots}")


# ===========================================================================
# PUNTO DE ENTRADA PRINCIPAL
# ===========================================================================

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Rango de ruido: de alto (caos) a bajo (orden), con 20 puntos.
    # Puedes modificar estos valores fácilmente aquí.
    # ------------------------------------------------------------------
    ETA_MAX  = 2.0    # Ruido máximo (movimiento casi aleatorio)
    ETA_MIN  = 0.05   # Ruido mínimo (alineación casi perfecta)
    N_PUNTOS = 20     # Número de valores de eta a explorar

    eta_values = np.linspace(ETA_MAX, ETA_MIN, N_PUNTOS)

    print("\n╔══════════════════════════════════════════════════════════╗")
    print("║    FASE 3 (v2) — RECOLECCIÓN MASIVA DE DATOS            ║")
    print("║    Modelos de Acoplamiento Global (Campo Medio)          ║")
    print("║  Barrido de eta: de {:.2f} a {:.2f} ({} puntos)    ║".format(ETA_MAX, ETA_MIN, N_PUNTOS))
    print("╚══════════════════════════════════════════════════════════╝")

    t0_global = time.time()

    # --- Barrido del Modelo de Vicsek All-to-All ---
    barrer_vicsek_all2all(
        eta_values=eta_values,
        N=150, L=20.0, R=2.0, v=0.5,
        pasos_termal=500,
        pasos_medicion=300,
        archivo_resumen="resultados_vicsek_all2all.csv",
        archivo_snapshots="snapshots_vicsek_all2all.csv",
    )

    # --- Barrido del Modelo Continuo Sin Radio ---
    barrer_continuo_sin_radio(
        eta_values=eta_values,
        N=150, L=20.0,
        k_align=5.0,
        pasos_termal=3000,
        pasos_medicion=1500,
        archivo_resumen="resultados_continuo_sin_radio.csv",
        archivo_snapshots="snapshots_continuo_sin_radio.csv",
    )

    total = time.time() - t0_global
    print(f"\n✓ Barrido completo finalizado en {total:.1f} segundos.")
    print("  Archivos generados:")
    print("    → resultados_vicsek_all2all.csv")
    print("    → resultados_continuo_sin_radio.csv")
    print("    → snapshots_vicsek_all2all.csv")
    print("    → snapshots_continuo_sin_radio.csv")
