import numpy as np
import csv
import time
import sys
from pathlib import Path
from fase1_vicsek import VicsekSimulation
from fase2_continuo import ContinuousSimulation

# ===========================================================================
# FUNCIÓN DE CÁLCULO DEL PARÁMETRO DE ORDEN POLAR (Phi)
# ===========================================================================

def calcular_phi_vicsek(angles: np.ndarray) -> float:
    N = len(angles)
    vx_sum = np.sum(np.cos(angles))
    vy_sum = np.sum(np.sin(angles))
    return np.sqrt(vx_sum**2 + vy_sum**2) / N


def calcular_phi_continuo(velocities: np.ndarray) -> float:
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

    print("\n")
    print("FASE 3 — RECOLECCIÓN MASIVA DE DATOS")
    print("Barrido de eta: de {:.2f} a {:.2f} ({} puntos)".format(ETA_MAX, ETA_MIN, N_PUNTOS))
    print("")

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
