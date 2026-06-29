# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

"""
Fase 4: Analisis Visual y Matematico
======================================
Analiza los datos generados por fase3_recoleccion.py  (modelos con radio local)
y por fase3_recoleccion_v2.py (modelos de campo medio / acoplamiento global).

Análisis realizados:
  1. Curvas de Transición de Fase  — Φ vs η para los 4 modelos.
  2. Susceptibilidad σ(Φ) vs η    — localización del punto crítico η_c.
  3. Correlación Espacial C(r)     — cómo se propaga el orden con la distancia.
  4. Ajuste del Decaimiento        — exponencial vs ley de potencias.

Archivos que necesita (generados por las Fases 3 y 3v2):
  - resultados_vicsek.csv
  - resultados_continuo.csv
  - resultados_vicsek_all2all.csv
  - resultados_continuo_sin_radio.csv
  - snapshots_vicsek.csv
  - snapshots_continuo.csv
  - snapshots_vicsek_all2all.csv
  - snapshots_continuo_sin_radio.csv

Archivos generados:
  - fase4_transicion_fase.png        : curvas Φ vs η  (4 modelos)
  - fase4_susceptibilidad.png        : σ(Φ) vs η      (4 modelos)
  - fase4_correlacion_espacial.png   : C(r) con ajustes en el punto crítico
  - fase4_resumen.txt                : tabla de resultados numéricos
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # backend sin GUI → funciona en cualquier entorno
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# PALETA Y ESTILO
# ─────────────────────────────────────────────────────────────────────────────

COLORES = {
    "vicsek":         "#4C9BE8",   # azul
    "continuo":       "#E85C4C",   # rojo
    "vicsek_all2all": "#5ECC8C",   # verde
    "continuo_sr":    "#F5A623",   # naranja
}

ETIQUETAS = {
    "vicsek":         "Vicsek Local  (R = 2.0)",
    "continuo":       "Continuo Local",
    "vicsek_all2all": "Vicsek All-to-All (campo medio)",
    "continuo_sr":    "Continuo Sin Radio (campo medio)",
}

plt.rcParams.update({
    "figure.facecolor":  "#12141C",
    "axes.facecolor":    "#1A1D2E",
    "axes.edgecolor":    "#3A3D5C",
    "axes.labelcolor":   "#E0E4FF",
    "axes.titlecolor":   "#E0E4FF",
    "xtick.color":       "#A0A4CC",
    "ytick.color":       "#A0A4CC",
    "text.color":        "#E0E4FF",
    "grid.color":        "#2A2D45",
    "grid.linewidth":    0.6,
    "legend.facecolor":  "#1A1D2E",
    "legend.edgecolor":  "#3A3D5C",
    "font.family":       "DejaVu Sans",
    "font.size":         11,
})


# ─────────────────────────────────────────────────────────────────────────────
# UTILIDADES
# ─────────────────────────────────────────────────────────────────────────────

def cargar_resultados(ruta: str) -> pd.DataFrame | None:
    """Carga un CSV de resultados. Devuelve None si el archivo no existe."""
    p = Path(ruta)
    if not p.exists():
        print(f"  [AVISO] No se encontró '{ruta}' — se omite este modelo.")
        return None
    df = pd.read_csv(ruta)
    df.sort_values("eta", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def cargar_snapshots(ruta: str) -> pd.DataFrame | None:
    """Carga un CSV de snapshots. Devuelve None si el archivo no existe."""
    p = Path(ruta)
    if not p.exists():
        print(f"  [AVISO] No se encontró '{ruta}' — sin correlación para este modelo.")
        return None
    return pd.read_csv(ruta)


def eta_critico_por_susceptibilidad(df: pd.DataFrame) -> float:
    """
    Devuelve el η en el que σ(Φ) es máxima.
    Es el estimador más robusto del punto crítico en sistemas finitos.
    """
    idx_max = df["phi_std"].idxmax()
    return float(df.loc[idx_max, "eta"])


def eta_critico_por_inflexion(df: pd.DataFrame) -> float:
    """
    Devuelve el η en el que d<Φ>/dη es máxima (inflexión de la curva de orden).
    Útil como verificación cruzada del punto crítico.
    """
    etas = df["eta"].values
    phis = df["phi_promedio"].values
    if len(etas) < 4:
        return float(etas[len(etas) // 2])
    try:
        phi_smooth = savgol_filter(phis, window_length=min(5, len(phis) - 2 | 1), polyorder=2)
    except Exception:
        phi_smooth = phis
    derivada = np.abs(np.gradient(phi_smooth, etas))
    return float(etas[np.argmax(derivada)])


# ─────────────────────────────────────────────────────────────────────────────
# CÁLCULO DE CORRELACIÓN ESPACIAL C(r)
# ─────────────────────────────────────────────────────────────────────────────

def vectores_unitarios_vicsek(snap_eta: pd.DataFrame) -> np.ndarray:
    """Devuelve vectores unitarios de velocidad para snapshots de Vicsek."""
    thetas = snap_eta["theta"].values
    return np.column_stack([np.cos(thetas), np.sin(thetas)])


def vectores_unitarios_continuo(snap_eta: pd.DataFrame) -> np.ndarray:
    """Devuelve vectores unitarios de velocidad para snapshots de modelo continuo."""
    vx = snap_eta["vx"].values
    vy = snap_eta["vy"].values
    speeds = np.sqrt(vx**2 + vy**2)
    speeds = np.where(speeds < 1e-9, 1.0, speeds)
    return np.column_stack([vx / speeds, vy / speeds])


def calcular_correlacion_espacial(
    posiciones: np.ndarray,   # (N, 2)
    direcciones: np.ndarray,  # (N, 2) vectores unitarios
    L: float,
    n_bins: int = 25,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcula la función de correlación espacial de velocidades:

        C(r) = < v̂_i · v̂_j >_{|r_ij| ≈ r}

    Usa condiciones de contorno periódicas (imagen mínima).

    Returns:
        r_centros  : centros de cada bin de distancia.
        C_r        : valor promedio de C(r) en cada bin.
        C_r_err    : error estándar de la media en cada bin.
    """
    N = len(posiciones)
    r_max = L / 2.0          # máxima distancia distinguible con PBC

    pares_r  = []
    pares_C  = []

    for i in range(N - 1):
        for j in range(i + 1, N):
            # distancia con imagen mínima
            dr = posiciones[j] - posiciones[i]
            dr -= L * np.round(dr / L)
            r  = np.linalg.norm(dr)
            if r < r_max:
                c_ij = float(np.dot(direcciones[i], direcciones[j]))
                pares_r.append(r)
                pares_C.append(c_ij)

    pares_r = np.array(pares_r)
    pares_C = np.array(pares_C)

    bins      = np.linspace(0, r_max, n_bins + 1)
    r_centros = 0.5 * (bins[:-1] + bins[1:])
    C_r       = np.zeros(n_bins)
    C_r_err   = np.zeros(n_bins)

    for k in range(n_bins):
        mask = (pares_r >= bins[k]) & (pares_r < bins[k + 1])
        vals = pares_C[mask]
        if len(vals) > 1:
            C_r[k]     = vals.mean()
            C_r_err[k] = vals.std() / np.sqrt(len(vals))
        elif len(vals) == 1:
            C_r[k]     = vals[0]
            C_r_err[k] = np.nan

    return r_centros, C_r, C_r_err


def obtener_snapshot_eta(snap_df: pd.DataFrame, eta_objetivo: float) -> pd.DataFrame:
    """Devuelve las filas del snapshot cuyo eta es el más cercano a eta_objetivo."""
    etas_disponibles = snap_df["eta"].unique()
    eta_mas_cercano  = etas_disponibles[np.argmin(np.abs(etas_disponibles - eta_objetivo))]
    return snap_df[snap_df["eta"] == eta_mas_cercano].copy()


# ─────────────────────────────────────────────────────────────────────────────
# MODELOS DE AJUSTE
# ─────────────────────────────────────────────────────────────────────────────

def modelo_exponencial(r, A, xi):
    """C(r) = A · exp(-r / ξ)   [decaimiento rápido, longitud de correlación ξ]"""
    return A * np.exp(-r / xi)


def modelo_ley_potencias(r, A, alpha):
    """C(r) = A · r^(-α)         [decaimiento lento, sistemas críticos]"""
    return A * np.power(r, -alpha)


def ajustar_correlacion(r: np.ndarray, C: np.ndarray) -> dict:
    """
    Ajusta C(r) con dos modelos:
      - Exponencial:   C(r) = A·exp(-r/ξ)
      - Ley potencias: C(r) = A·r^{-α}

    Retorna un dict con los parámetros y el R² de cada ajuste.
    """
    resultado = {}

    # Filtramos bins sin datos o con r ≈ 0
    mask = np.isfinite(C) & (r > 0.3)
    r_f, C_f = r[mask], C[mask]

    if len(r_f) < 3:
        return resultado

    # ── Ajuste exponencial ──────────────────────────────────────────────────
    try:
        popt_exp, _ = curve_fit(
            modelo_exponencial, r_f, C_f,
            p0=[C_f[0] if C_f[0] > 0 else 0.5, r_f.max() / 2],
            bounds=([0, 0.01], [2.0, r_f.max() * 3]),
            maxfev=5000,
        )
        C_pred_exp = modelo_exponencial(r_f, *popt_exp)
        ss_res = np.sum((C_f - C_pred_exp) ** 2)
        ss_tot = np.sum((C_f - C_f.mean()) ** 2)
        r2_exp = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        resultado["exp"] = {"A": popt_exp[0], "xi": popt_exp[1], "R2": r2_exp}
    except Exception:
        pass

    # ── Ajuste ley de potencias  (solo si C > 0) ────────────────────────────
    mask_pos = mask & (C > 1e-6)
    r_p, C_p = r[mask_pos], C[mask_pos]
    if len(r_p) >= 3:
        try:
            popt_pow, _ = curve_fit(
                modelo_ley_potencias, r_p, C_p,
                p0=[C_p[0], 0.5],
                bounds=([0, 0.01], [5.0, 5.0]),
                maxfev=5000,
            )
            C_pred_pow = modelo_ley_potencias(r_p, *popt_pow)
            ss_res = np.sum((C_p - C_pred_pow) ** 2)
            ss_tot = np.sum((C_p - C_p.mean()) ** 2)
            r2_pow = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            resultado["pow"] = {"A": popt_pow[0], "alpha": popt_pow[1], "R2": r2_pow}
        except Exception:
            pass

    return resultado


# ─────────────────────────────────────────────────────────────────────────────
# FIGURA 1: CURVAS DE TRANSICIÓN DE FASE
# ─────────────────────────────────────────────────────────────────────────────

def graficar_transicion_fase(datasets: dict, eta_criticos: dict) -> None:
    """
    Grafica Φ vs η para todos los modelos disponibles en un panel de 2x2.
    Los paneles superiores corresponden a los modelos con radio local (Fase 3 original)
    y los inferiores a los modelos de campo medio (Fase 3 v2).
    """
    print("\n  → Generando Figura 1: Curvas de Transición de Fase...")

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        "Transición de Fase: Parámetro de Orden Polar $\\langle\\Phi\\rangle$ vs Ruido $\\eta$",
        fontsize=15, fontweight="bold", y=0.98,
    )

    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    # Configuración de cada subplot: (clave, posición en gridspec)
    config = [
        ("vicsek",         gs[0, 0]),
        ("continuo",       gs[0, 1]),
        ("vicsek_all2all", gs[1, 0]),
        ("continuo_sr",    gs[1, 1]),
    ]

    subtitulos = {
        "vicsek":         "Vicsek Discreto — Radio local R = 2.0",
        "continuo":       "Materia Activa Continua — Radio local",
        "vicsek_all2all": "Vicsek All-to-All — Campo Medio",
        "continuo_sr":    "Continuo Sin Radio — Campo Medio",
    }

    for clave, pos in config:
        df = datasets.get(clave)
        ax = fig.add_subplot(pos)
        ax.set_title(subtitulos[clave], fontsize=10, pad=6)
        ax.set_xlabel("Ruido $\\eta$")
        ax.set_ylabel("$\\langle\\Phi\\rangle$")
        ax.set_xlim(0, 2.1)
        ax.set_ylim(-0.05, 1.1)
        ax.grid(True, alpha=0.3)

        if df is None:
            ax.text(0.5, 0.5, "Datos no disponibles\n(ejecutar fase3)", ha="center",
                    va="center", transform=ax.transAxes, color="#888", fontsize=10)
            continue

        color = COLORES[clave]
        etas  = df["eta"].values
        phis  = df["phi_promedio"].values
        stds  = df["phi_std"].values

        # Banda de incertidumbre
        ax.fill_between(etas, phis - stds, phis + stds,
                         alpha=0.20, color=color)
        # Línea principal
        ax.plot(etas, phis, "o-", color=color, lw=2, ms=6,
                label=ETIQUETAS[clave])

        # Marca del punto crítico
        if clave in eta_criticos:
            eta_c = eta_criticos[clave]["susceptibilidad"]
            ax.axvline(eta_c, color="#FFD700", lw=1.5, ls="--", alpha=0.85,
                       label=f"$\\eta_c$ = {eta_c:.3f}")

        ax.legend(fontsize=8, loc="upper left")

        # Anotación de régimen
        ax.text(0.05,  0.08, "Orden",  transform=ax.transAxes,
                color="#88FF88", fontsize=9, alpha=0.85)
        ax.text(0.72,  0.88, "Caos",   transform=ax.transAxes,
                color="#FF8888", fontsize=9, alpha=0.85)

    plt.savefig("fase4_transicion_fase.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("     ✓ fase4_transicion_fase.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURA 2: SUSCEPTIBILIDAD σ(Φ) vs η
# ─────────────────────────────────────────────────────────────────────────────

def graficar_susceptibilidad(datasets: dict, eta_criticos: dict) -> None:
    """
    Grafica la desviación estándar de Φ (proxy de susceptibilidad) vs η.
    El pico identifica el punto crítico η_c.
    """
    print("  → Generando Figura 2: Susceptibilidad σ(Φ) vs η...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        "Susceptibilidad $\\sigma(\\Phi)$ vs Ruido $\\eta$ — Localización del Punto Crítico $\\eta_c$",
        fontsize=14, fontweight="bold",
    )

    grupos = {
        "Modelos con Radio Local (Fase 3)":    ["vicsek", "continuo"],
        "Modelos de Campo Medio (Fase 3 v2)":  ["vicsek_all2all", "continuo_sr"],
    }

    for ax, (titulo, claves) in zip(axes, grupos.items()):
        ax.set_title(titulo, fontsize=11, pad=8)
        ax.set_xlabel("Ruido $\\eta$")
        ax.set_ylabel("$\\sigma(\\Phi)$  [desviación estándar]")
        ax.set_xlim(0, 2.1)
        ax.grid(True, alpha=0.3)

        for clave in claves:
            df = datasets.get(clave)
            if df is None:
                continue

            color = COLORES[clave]
            etas  = df["eta"].values
            stds  = df["phi_std"].values

            ax.plot(etas, stds, "o-", color=color, lw=2.2, ms=7,
                    label=ETIQUETAS[clave])

            # Marcar pico
            if clave in eta_criticos:
                eta_c = eta_criticos[clave]["susceptibilidad"]
                idx_max = np.argmax(stds)
                ax.annotate(
                    f"$\\eta_c = {eta_c:.3f}$",
                    xy=(etas[idx_max], stds[idx_max]),
                    xytext=(etas[idx_max] + 0.15, stds[idx_max] + 0.005),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.5),
                    color=color, fontsize=9,
                )

        ax.legend(fontsize=9, loc="upper right")

    plt.tight_layout()
    plt.savefig("fase4_susceptibilidad.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("     ✓ fase4_susceptibilidad.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURA 3: CORRELACIÓN ESPACIAL C(r) CON AJUSTES
# ─────────────────────────────────────────────────────────────────────────────

def graficar_correlacion(
    datasets_snap: dict,
    eta_criticos: dict,
    L: float = 20.0,
) -> dict:
    """
    Calcula y grafica C(r) en el punto crítico para cada modelo disponible.
    Ajusta decaimiento exponencial y ley de potencias.
    Retorna dict con los parámetros de ajuste para el resumen numérico.
    """
    print("  → Generando Figura 3: Correlación Espacial C(r)...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Función de Correlación Espacial $C(r)$ en el Punto Crítico $\\eta_c$",
        fontsize=14, fontweight="bold", y=0.99,
    )

    config_snap = [
        ("vicsek",         axes[0, 0], "vicsek",   False),
        ("continuo",       axes[0, 1], "continuo",  True),
        ("vicsek_all2all", axes[1, 0], "vicsek",   False),
        ("continuo_sr",    axes[1, 1], "continuo",  True),
    ]

    subtitulos_corr = {
        "vicsek":         "Vicsek Local — en $\\eta_c$",
        "continuo":       "Continuo Local — en $\\eta_c$",
        "vicsek_all2all": "Vicsek All-to-All — en $\\eta_c$",
        "continuo_sr":    "Continuo Sin Radio — en $\\eta_c$",
    }

    ajustes_todos = {}

    for clave, ax, tipo_snap, es_continuo in config_snap:
        snap_df = datasets_snap.get(clave)
        ax.set_title(subtitulos_corr[clave], fontsize=10, pad=6)
        ax.set_xlabel("Distancia $r$")
        ax.set_ylabel("$C(r)$")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="#555", lw=0.8, ls="--")

        if snap_df is None or clave not in eta_criticos:
            ax.text(0.5, 0.5, "Datos no disponibles", ha="center", va="center",
                    transform=ax.transAxes, color="#888", fontsize=10)
            continue

        eta_c = eta_criticos[clave]["susceptibilidad"]
        snap_eta = obtener_snapshot_eta(snap_df, eta_c)
        eta_real = float(snap_eta["eta"].iloc[0])

        pos = snap_eta[["x", "y"]].values
        dir_vec = (
            vectores_unitarios_continuo(snap_eta)
            if es_continuo
            else vectores_unitarios_vicsek(snap_eta)
        )

        r, C, C_err = calcular_correlacion_espacial(pos, dir_vec, L)

        color = COLORES[clave]
        ax.errorbar(r, C, yerr=C_err, fmt="o", color=color, ms=5,
                    elinewidth=1, capsize=3, label=f"Datos (η={eta_real:.3f})")

        ajustes = ajustar_correlacion(r, C)
        ajustes_todos[clave] = {"eta_c": eta_c, "ajustes": ajustes}

        r_fino = np.linspace(r[r > 0.3][0] if any(r > 0.3) else 0.5, r.max(), 200)

        # Curva exponencial
        if "exp" in ajustes:
            p = ajustes["exp"]
            C_exp = modelo_exponencial(r_fino, p["A"], p["xi"])
            ax.plot(r_fino, C_exp, "--", color="#FFD700", lw=1.8,
                    label=f"Exp:  ξ={p['xi']:.2f},  R²={p['R2']:.3f}")

        # Curva ley de potencias
        if "pow" in ajustes:
            p = ajustes["pow"]
            r_pos = r_fino[r_fino > 0]
            C_pow = modelo_ley_potencias(r_pos, p["A"], p["alpha"])
            ax.plot(r_pos, C_pow, "-.", color="#FF6BD6", lw=1.8,
                    label=f"Ley pot: α={p['alpha']:.2f},  R²={p['R2']:.3f}")

        # Interpretación del tipo de decaimiento
        mejor = _mejor_ajuste(ajustes)
        if mejor == "pow":
            tipo_txt = "⚡ Ley de Potencias (crítico)"
            col_txt  = "#FF6BD6"
        elif mejor == "exp":
            tipo_txt = "Exponencial (subcrítico)"
            col_txt  = "#FFD700"
        else:
            tipo_txt = "Sin ajuste"
            col_txt  = "#888"
        ax.text(0.97, 0.95, tipo_txt, transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, color=col_txt,
                bbox=dict(boxstyle="round,pad=0.3", fc="#1A1D2E", ec=col_txt, lw=1))

        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlim(left=0)

    plt.tight_layout()
    plt.savefig("fase4_correlacion_espacial.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("     ✓ fase4_correlacion_espacial.png")

    return ajustes_todos


def _mejor_ajuste(ajustes: dict) -> str:
    """Devuelve 'pow' o 'exp' según cuál tiene mayor R²."""
    r2_exp = ajustes.get("exp", {}).get("R2", -np.inf)
    r2_pow = ajustes.get("pow", {}).get("R2", -np.inf)
    if r2_pow > r2_exp and r2_pow > 0:
        return "pow"
    if r2_exp > 0:
        return "exp"
    return "ninguno"


# ─────────────────────────────────────────────────────────────────────────────
# FIGURA 4: COMPARATIVA CONJUNTA (panel resumen)
# ─────────────────────────────────────────────────────────────────────────────

def graficar_comparativa(datasets: dict, eta_criticos: dict) -> None:
    """
    Un único panel que superpone las 4 curvas Φ vs η para comparar
    directamente los modelos locales y de campo medio.
    """
    print("  → Generando Figura 4: Comparativa entre los 4 modelos...")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title(
        "Comparativa: Transición de Fase en los 4 Modelos",
        fontsize=14, fontweight="bold",
    )
    ax.set_xlabel("Ruido $\\eta$", fontsize=13)
    ax.set_ylabel("$\\langle\\Phi\\rangle$ — Parámetro de Orden Polar", fontsize=13)
    ax.set_xlim(0, 2.1)
    ax.set_ylim(-0.05, 1.12)
    ax.grid(True, alpha=0.3)

    estilos = {
        "vicsek":         ("o-",  2.0),
        "continuo":       ("s--", 2.0),
        "vicsek_all2all": ("^-",  2.0),
        "continuo_sr":    ("D--", 2.0),
    }

    for clave, (fmt, lw) in estilos.items():
        df = datasets.get(clave)
        if df is None:
            continue
        color = COLORES[clave]
        ax.plot(df["eta"], df["phi_promedio"], fmt,
                color=color, lw=lw, ms=7, label=ETIQUETAS[clave])
        ax.fill_between(
            df["eta"],
            df["phi_promedio"] - df["phi_std"],
            df["phi_promedio"] + df["phi_std"],
            alpha=0.10, color=color,
        )
        if clave in eta_criticos:
            eta_c = eta_criticos[clave]["susceptibilidad"]
            ax.axvline(eta_c, color=color, lw=1.0, ls=":", alpha=0.6)

    # Leyenda y anotaciones
    ax.legend(fontsize=10, loc="upper left", framealpha=0.9)
    ax.text(0.03, 0.12, "← Régimen Ordenado", transform=ax.transAxes,
            color="#88FF88", fontsize=10)
    ax.text(0.68, 0.88, "Régimen Desordenado →", transform=ax.transAxes,
            color="#FF8888", fontsize=10)

    plt.tight_layout()
    plt.savefig("fase4_comparativa.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("     ✓ fase4_comparativa.png")


# ─────────────────────────────────────────────────────────────────────────────
# RESUMEN NUMÉRICO
# ─────────────────────────────────────────────────────────────────────────────

def escribir_resumen(datasets: dict, eta_criticos: dict, ajustes_todos: dict) -> None:
    """Escribe un archivo de texto con todos los resultados numéricos clave."""
    print("  → Generando resumen numérico...")

    lineas = []
    sep = "=" * 70

    lineas.append(sep)
    lineas.append("  FASE 4 — RESUMEN DE RESULTADOS NUMÉRICOS")
    lineas.append(sep)

    for clave, etiq in ETIQUETAS.items():
        df = datasets.get(clave)
        lineas.append(f"\n{'─'*70}")
        lineas.append(f"  MODELO: {etiq}")
        lineas.append(f"{'─'*70}")

        if df is None:
            lineas.append("  [Sin datos — ejecutar la fase de recolección correspondiente]")
            continue

        # ── Estadísticas generales ──────────────────────────────────────────
        phi_max = df["phi_promedio"].max()
        phi_min = df["phi_promedio"].min()
        lineas.append(f"  Número de puntos eta     : {len(df)}")
        lineas.append(f"  Rango de eta             : [{df['eta'].min():.3f}, {df['eta'].max():.3f}]")
        lineas.append(f"  <Phi> máximo (bajo ruido): {phi_max:.5f}")
        lineas.append(f"  <Phi> mínimo (alto ruido): {phi_min:.5f}")

        # ── Punto crítico ───────────────────────────────────────────────────
        if clave in eta_criticos:
            ec = eta_criticos[clave]
            lineas.append(f"\n  Punto crítico (susceptibilidad máxima):")
            lineas.append(f"    η_c = {ec['susceptibilidad']:.4f}")
            lineas.append(f"    η_c = {ec['inflexion']:.4f}  (verificación por inflexión)")

        # ── Parámetros de ajuste de C(r) ────────────────────────────────────
        if clave in ajustes_todos:
            info = ajustes_todos[clave]
            lineas.append(f"\n  Correlación espacial (evaluada en η = {info['eta_c']:.4f}):")
            aj = info["ajustes"]
            if "exp" in aj:
                e = aj["exp"]
                lineas.append(f"    Ajuste exponencial  : ξ = {e['xi']:.4f}  |  R² = {e['R2']:.4f}")
            if "pow" in aj:
                p = aj["pow"]
                lineas.append(f"    Ajuste ley potencias: α = {p['alpha']:.4f}  |  R² = {p['R2']:.4f}")
            mejor = _mejor_ajuste(aj)
            if mejor == "pow":
                lineas.append("    ➜ Decaimiento tipo LEY DE POTENCIAS (señal de criticalidad)")
            elif mejor == "exp":
                lineas.append("    ➜ Decaimiento EXPONENCIAL (sistema subcrítico)")

    lineas.append(f"\n{sep}")
    lineas.append("  Archivos generados:")
    for f in ["fase4_transicion_fase.png", "fase4_susceptibilidad.png",
              "fase4_correlacion_espacial.png", "fase4_comparativa.png"]:
        lineas.append(f"    → {f}")
    lineas.append(sep + "\n")

    texto = "\n".join(lineas)
    with open("fase4_resumen.txt", "w", encoding="utf-8") as f:
        f.write(texto)
    print("     ✓ fase4_resumen.txt")
    print()
    print(texto)


# ─────────────────────────────────────────────────────────────────────────────
# PUNTO DE ENTRADA PRINCIPAL
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n" + "="*60)
    print("    FASE 4 -- ANALISIS VISUAL Y MATEMATICO")
    print("    (Modelos locales + Campo Medio)")
    print("="*60 + "\n")

    # ── 1. Cargar todos los resultados ───────────────────────────────────────
    print("[ 1/5 ] Cargando archivos CSV...")
    datasets = {
        "vicsek":         cargar_resultados("resultados_vicsek.csv"),
        "continuo":       cargar_resultados("resultados_continuo.csv"),
        "vicsek_all2all": cargar_resultados("resultados_vicsek_all2all.csv"),
        "continuo_sr":    cargar_resultados("resultados_continuo_sin_radio.csv"),
    }

    datasets_snap = {
        "vicsek":         cargar_snapshots("snapshots_vicsek.csv"),
        "continuo":       cargar_snapshots("snapshots_continuo.csv"),
        "vicsek_all2all": cargar_snapshots("snapshots_vicsek_all2all.csv"),
        "continuo_sr":    cargar_snapshots("snapshots_continuo_sin_radio.csv"),
    }

    # ── 2. Calcular puntos críticos ──────────────────────────────────────────
    print("\n[ 2/5 ] Calculando puntos críticos η_c...")
    eta_criticos = {}
    for clave, df in datasets.items():
        if df is None:
            continue
        ec_sus = eta_critico_por_susceptibilidad(df)
        ec_inf = eta_critico_por_inflexion(df)
        eta_criticos[clave] = {
            "susceptibilidad": ec_sus,
            "inflexion":       ec_inf,
        }
        print(f"  {ETIQUETAS[clave]}")
        print(f"    η_c (σ_max)   = {ec_sus:.4f}")
        print(f"    η_c (inflexión)= {ec_inf:.4f}")

    # ── 3. Gráficas principales ──────────────────────────────────────────────
    print("\n[ 3/5 ] Generando figuras de transición y susceptibilidad...")
    graficar_transicion_fase(datasets, eta_criticos)
    graficar_susceptibilidad(datasets, eta_criticos)
    graficar_comparativa(datasets, eta_criticos)

    # ── 4. Correlación espacial ──────────────────────────────────────────────
    print("\n[ 4/5 ] Calculando correlaciones espaciales C(r)...")
    print("        (esto puede tomar unos segundos por modelo...)")
    ajustes_todos = graficar_correlacion(datasets_snap, eta_criticos, L=20.0)

    # ── 5. Resumen numérico ──────────────────────────────────────────────────
    print("\n[ 5/5 ] Escribiendo resumen numérico...")
    escribir_resumen(datasets, eta_criticos, ajustes_todos)

    print("✓ Análisis completado. Archivos guardados en el directorio actual.")
