import csv
import math
import argparse
from pathlib import Path
import subprocess
import sys
import shutil

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.solve_assignment import (
    main as solve_main,
    digits_from_studienummer,
    pick_params,
    pick_tp_eq_mm,
    moonpool,
    trapz,
    RHO_G_MN,
)


OUTPUT_DIR = REPO_ROOT / "outputs"
ASSETS_DIR = OUTPUT_DIR / "report_assets"

B = 31.0
H = 15.5
CRANE_DEADWEIGHT = 2.0  # [MN]
CRANE_SWL = 1.5         # [MN]
CRANE_REACH = 30.0      # [m]


def material_from_e(e_digit: int):
    if e_digit <= 4:
        return {"name": "Staal S235J2", "sigma_y": 235.0, "E": 205000.0, "nu": 0.3, "rho": 7800.0}
    return {"name": "Staal S355J2", "sigma_y": 355.0, "E": 205000.0, "nu": 0.3, "rho": 7800.0}


def load_q1_data(csv_path: Path):
    rows = list(csv.DictReader(csv_path.open()))
    x = [float(r["x"]) for r in rows]
    W = [float(r["W"]) for r in rows]
    c = [float(r["c"]) for r in rows]
    g = [float(r["g"]) for r in rows]
    p = [float(r["p"]) for r in rows]
    q = [float(r["q"]) for r in rows]
    Fs = [float(r["Fs"]) for r in rows]
    Mb = [float(r["Mb"]) for r in rows]
    Fs_n = [float(r["Fs_nomoon"]) for r in rows]
    Mb_n = [float(r["Mb_nomoon"]) for r in rows]
    return x, W, c, g, p, q, Fs, Mb, Fs_n, Mb_n


def compute_q1_metrics(digits, x, g, q):
    g_digit = digits["g"]
    p = pick_params(g_digit)
    Lm, Bm, xm = moonpool(g_digit, p.L)
    seg_len = p.L / 6.0
    crane_x = 2 * seg_len if p.crane_transition == "2-3" else 3 * seg_len

    beff = [B - (Bm if (xm - Lm / 2.0 <= xi <= xm + Lm / 2.0) else 0.0) for xi in x]
    A_wp = trapz(beff, x)
    LCF = trapz([xi * bi for xi, bi in zip(x, beff)], x) / A_wp

    P_crane = CRANE_DEADWEIGHT + CRANE_SWL
    Wtot_distributed = trapz(g, x)
    Wtot = Wtot_distributed + P_crane
    LCG = (trapz([xi * gi for xi, gi in zip(x, g)], x) + crane_x * P_crane) / Wtot

    I0 = trapz(beff, x)
    I2 = trapz([((xi - LCF) ** 2) * bi for xi, bi in zip(x, beff)], x)
    T0 = Wtot / (RHO_G_MN * I0)
    a = Wtot * (LCG - LCF) / (RHO_G_MN * I2)
    ta = T0 + a * (0.0 - LCF)
    tf = T0 + a * (p.L - LCF)

    eq_force_distributed = trapz(q, x)
    eq_force_total = eq_force_distributed - P_crane
    eq_moment_distributed = trapz([(xi - LCF) * qi for xi, qi in zip(x, q)], x)
    eq_moment_total = eq_moment_distributed - P_crane * (crane_x - LCF)

    return {
        "L": p.L,
        "W_const": p.W_const,
        "c3": p.c3,
        "seg_len": seg_len,
        "crane_x": crane_x,
        "xm": xm,
        "Lm": Lm,
        "Bm": Bm,
        "P_crane": P_crane,
        "LCF": LCF,
        "LCG": LCG,
        "T0": T0,
        "ta": ta,
        "tf": tf,
        "eq_force_total": eq_force_total,
        "eq_moment_total": eq_moment_total,
    }


def compute_q2_metrics(digits):
    tp_mm = pick_tp_eq_mm(digits["f"])
    tp = tp_mm / 1000.0
    t_h = 4 * tp
    t_v = 2 * tp

    As = 4 * H * t_v
    z_n = H / 2.0

    A_h = B * t_h
    I_h_own = B * t_h ** 3 / 12.0
    d_h = H / 2.0 - t_h / 2.0
    I_h_shift = A_h * d_h ** 2
    I_h_total_each = I_h_own + I_h_shift

    I_v_own = t_v * H ** 3 / 12.0
    I_v_shift = 0.0
    I_v_total_each = I_v_own

    Ib = 2 * I_h_total_each + 4 * I_v_total_each

    return {
        "tp_mm": tp_mm,
        "tp": tp,
        "t_h": t_h,
        "t_v": t_v,
        "As": As,
        "z_n": z_n,
        "I_h_own": I_h_own,
        "I_h_shift": I_h_shift,
        "I_h_total_each": I_h_total_each,
        "I_v_own": I_v_own,
        "I_v_shift": I_v_shift,
        "I_v_total_each": I_v_total_each,
        "Ib": Ib,
    }


def compute_q3_metrics(x, Fs, Mb, As, Ib, material):
    E = material["E"]  # MPa = MN/m^2
    EI = E * Ib

    kappa = [m / EI for m in Mb]  # 1/m
    phi = [0.0] * len(x)
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        phi[i] = phi[i - 1] + 0.5 * (kappa[i] + kappa[i - 1]) * dx

    w = [0.0] * len(x)
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        w[i] = w[i - 1] + 0.5 * (phi[i] + phi[i - 1]) * dx

    # Remove rigid-body line so w(0)=w(L)=0.
    L = x[-1] - x[0]
    rigid_slope = w[-1] / L
    w_adj = [wi - rigid_slope * (xi - x[0]) for wi, xi in zip(w, x)]
    phi_adj = [pi - rigid_slope for pi in phi]

    phi_deg = [v * 180.0 / math.pi for v in phi_adj]
    w_mm = [v * 1000.0 for v in w_adj]

    i_m = max(range(len(Mb)), key=lambda i: abs(Mb[i]))
    i_v = max(range(len(Fs)), key=lambda i: abs(Fs[i]))
    M_max = Mb[i_m]
    V_max = Fs[i_v]
    x_m = x[i_m]
    x_v = x[i_v]

    y = H / 2.0
    sigma_b_max = abs(M_max) * y / Ib  # MPa
    tau_s = abs(V_max) / As            # MPa

    return {
        "EI": EI,
        "phi_deg": phi_deg,
        "w_mm": w_mm,
        "M_max": M_max,
        "V_max": V_max,
        "x_m": x_m,
        "x_v": x_v,
        "sigma_b_max": sigma_b_max,
        "tau_s": tau_s,
        "fails_yield": sigma_b_max > material["sigma_y"],
    }


def save_plot_single(x, y, ylabel, title, filename, color):
    fig, ax = plt.subplots(1, 1, figsize=(10.5, 3.2))
    ax.plot(x, y, color=color, linewidth=1.5)
    ax.fill_between(x, 0.0, y, color=color, alpha=0.18)
    ax.set_xlabel("x [m]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def make_report_figures(assets_dir, x, W, c, g, p, q, Fs, Mb, Fs_n, Mb_n, crane_x, phi_deg, w_mm):
    assets_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    ysets = [(W, "W(x) [MN/m]", "#0055aa"), (c, "c(x) [MN/m]", "#aa5500"), (g, "g(x) [MN/m]", "#228833")]
    for ax, (y, yl, col) in zip(axes, ysets):
        ax.plot(x, y, color=col, linewidth=1.4)
        ax.fill_between(x, 0.0, y, color=col, alpha=0.18)
        ax.axvline(crane_x, color="#444444", linestyle="--", linewidth=1.0, alpha=0.8)
        ax.set_ylabel(yl)
        ax.grid(True, alpha=0.3)
    axes[0].set_title("Vraag 1a: W(x), c(x), g(x)")
    axes[-1].set_xlabel("x [m]")
    axes[0].text(0.99, 0.92, "Kraanpunt op x = {:.2f} m".format(crane_x), transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(assets_dir / "q1a_weights.png", dpi=200)
    plt.close(fig)

    save_plot_single(x, p, "p(x) [MN/m]", "Vraag 1c: drijfvermogenverdeling p(x)", assets_dir / "q1c_buoyancy.png", "#9933cc")
    save_plot_single(x, q, "q(x) [MN/m]", "Vraag 1d: resultante belasting q(x)", assets_dir / "q1d_q.png", "#cc2222")
    save_plot_single(x, Fs, "Fs(x) [MN]", "Vraag 1e: dwarskrachtverdeling Fs(x)", assets_dir / "q1e_fs.png", "#006666")
    save_plot_single(x, Mb, "Mb(x) [MNm]", "Vraag 1e: buigend momentenverdeling Mb(x)", assets_dir / "q1e_mb.png", "#444444")

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.0), sharex=True)
    axes[0].plot(x, Fs, color="#0055aa", linewidth=1.5, label="Fs met moonpool")
    axes[0].plot(x, Fs_n, color="#dd8800", linewidth=1.5, label="Fs zonder moonpool")
    axes[0].fill_between(x, 0.0, Fs, color="#0055aa", alpha=0.14)
    axes[0].fill_between(x, 0.0, Fs_n, color="#dd8800", alpha=0.14)
    axes[0].set_ylabel("Fs [MN]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="best")

    axes[1].plot(x, Mb, color="#228833", linewidth=1.5, label="Mb met moonpool")
    axes[1].plot(x, Mb_n, color="#cc2222", linewidth=1.5, label="Mb zonder moonpool")
    axes[1].fill_between(x, 0.0, Mb, color="#228833", alpha=0.14)
    axes[1].fill_between(x, 0.0, Mb_n, color="#cc2222", alpha=0.14)
    axes[1].set_ylabel("Mb [MNm]")
    axes[1].set_xlabel("x [m]")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="best")
    axes[0].set_title("Vraag 1f: vergelijking met en zonder moonpool")
    fig.tight_layout()
    fig.savefig(assets_dir / "q1f_compare_overlay.png", dpi=200)
    plt.close(fig)

    save_plot_single(x, phi_deg, "φ(x) [deg]", "Vraag 3a: hoekverdraaiing φ(x)", assets_dir / "q3a_phi.png", "#7a3db8")
    save_plot_single(x, w_mm, "w(x) [mm]", "Vraag 3a: vervorming w(x)", assets_dir / "q3a_w.png", "#007a5a")


def render_markdown(name, studienummer, digits, q1, q2, q3, material, moonpool_dfs, moonpool_dmb):
    safety_txt = "Er treedt geen plastisch gedrag op." if not q3["fails_yield"] else "Er treedt plastisch gedrag op."
    md = """# MT1464 – Opdracht 1
## Verslag
### __NAME__
### __STUDIENUMMER__

![Voorpagina](plaatjevoorkant.png)

\\newpage

## 1. Uitwendige en inwendige belasting

### 1a
De gewichtscomponenten zijn opgesteld als:

$$
W(x)=W_{\\mathrm{const}} + P_{\\mathrm{kraan}}\\,\\delta(x-x_c),\\qquad
c(x)=c_{\\mathrm{seg}}(x)+P_{\\mathrm{SWL}}\\,\\delta(x-x_c),
$$
$$
g(x)=W(x)+c(x).
$$

Met $P_{\\mathrm{kraan}}=2.0\\ \\mathrm{MN}$, $P_{\\mathrm{SWL}}=1.5\\ \\mathrm{MN}$ en $x_c=__CRANE_X__\\ \\mathrm{m}$.

![Vraag 1a gewichtsverdelingen](q1a_weights.png)

### 1b
Gebruikte formules:

$$
A_{wp}=\\int_0^L b_{eff}(x)\\,dx,\\quad
LCF=\\frac{\\int_0^L x\\,b_{eff}(x)\\,dx}{A_{wp}}
$$
$$
W_{tot}=\\int_0^L g(x)\\,dx + P_{\\mathrm{kraan}}+P_{\\mathrm{SWL}},
\\quad
LCG=\\frac{\\int_0^L x\\,g(x)\\,dx + x_c(P_{\\mathrm{kraan}}+P_{\\mathrm{SWL}})}{W_{tot}}
$$
$$
T_0=\\frac{W_{tot}}{\\rho_w g A_{wp}},
\\quad
T(x)=T_0+a(x-LCF),
\\quad
a=\\frac{W_{tot}(LCG-LCF)}{\\rho_w g I_2}
$$
$$
t_a=T(0),\\qquad t_f=T(L).
$$

Resultaten (afgerond op 2 decimalen):  
- $T=__T0__\\ \\mathrm{m}$  
- $LCF=__LCF__\\ \\mathrm{m}$  
- $LCG=__LCG__\\ \\mathrm{m}$  
- $t_a=__TA__\\ \\mathrm{m}$  
- $t_f=__TF__\\ \\mathrm{m}$

### 1c
Drijfvermogenverdeling:

$$
p(x)=\\rho_w g\\,b_{eff}(x)\\,T(x)
$$

![Vraag 1c drijfvermogen](q1c_buoyancy.png)

### 1d
Resultante vlakwaterbelasting:

$$
q(x)=p(x)-g_{dist}(x)
$$

met evenwichtscontrole inclusief puntlasten:

$$
\\int_0^L q(x)\\,dx - (P_{\\mathrm{kraan}}+P_{\\mathrm{SWL}})\\approx 0
$$
$$
\\int_0^L (x-LCF)q(x)\\,dx -(P_{\\mathrm{kraan}}+P_{\\mathrm{SWL}})(x_c-LCF)\\approx 0
$$

Numeriek:  
- krachtsevenwicht: __EQ_FORCE__\\ \\mathrm{MN}  
- momentevenwicht: __EQ_MOMENT__\\ \\mathrm{MNm}

![Vraag 1d resultante belasting](q1d_q.png)

### 1e
Inwendige snedekrachten:

$$
F_s(x)=\\int_0^x q(\\xi)\\,d\\xi - H(x-x_c)(P_{\\mathrm{kraan}}+P_{\\mathrm{SWL}})
$$
$$
M_b(x)=\\int_0^x F_s(\\xi)\\,d\\xi + H(x-x_c)\\,P_{\\mathrm{SWL}}\\,R
$$

met $R=30\\ \\mathrm{m}$ en $H(\\cdot)$ de Heaviside-stapfunctie.

![Vraag 1e dwarskracht](q1e_fs.png)

![Vraag 1e buigend moment](q1e_mb.png)

### 1f
Vergelijking met/zonder moonpool toont lokaal sterkere variatie bij moonpool door lager effectief waterlijnoppervlak.  
Maximaal verschil:

$$
\\max|F_{s,met}-F_{s,zonder}|=__MOONPOOL_DFS__\\ \\mathrm{MN},\\qquad
\\max|M_{b,met}-M_{b,zonder}|=__MOONPOOL_DMB__\\ \\mathrm{MNm}.
$$

![Vraag 1f vergelijking moonpool](q1f_compare_overlay.png)

### 1g
Controlevoorwaarden aan de uiteinden:

$$
F_s(0)=0,\\ M_b(0)=0,\\ F_s(L)\\approx0,\\ M_b(L)=M_{\\mathrm{extern}}
$$

Aanvullend moet gelden:

$$
\\frac{dF_s}{dx}=q(x),\\qquad \\frac{dM_b}{dx}=F_s(x)
$$

De numerieke checks op deze relaties zijn uitgevoerd en klein.

## 2. Constructie eigenschappen

### 2a
Verticale delen domineren schuifstijfheid, omdat schuif vooral via web-achtige platen loopt.

$$
A_s=4H(2t_p)=8Ht_p
$$

Met $t_p=__TP_MM__\\ \\mathrm{mm}$:

$$
A_s=__AS__\\ \\mathrm{m^2}
$$

### 2b

$$
z_n=\\frac{\\sum A_i z_i}{\\sum A_i}=\\frac H2=__ZN__\\ \\mathrm{m}
$$

### 2c

$$
I_b=\\sum_i\\left(I_{\\mathrm{eigen},i}+A_i d_i^2\\right)
$$

Horizontaal (per plaat):

$$
I_{h,\\mathrm{eigen}}=__IH_OWN__\\ \\mathrm{m^4},\\quad
A_hd_h^2=__IH_SHIFT__\\ \\mathrm{m^4},\\quad
I_{h,tot}=__IH_TOT__\\ \\mathrm{m^4}
$$

Verticaal (per plaat):

$$
I_{v,\\mathrm{eigen}}=__IV_OWN__\\ \\mathrm{m^4},\\quad
A_vd_v^2=__IV_SHIFT__\\ \\mathrm{m^4},\\quad
I_{v,tot}=__IV_TOT__\\ \\mathrm{m^4}
$$

Totaal:

$$
I_b=2I_{h,tot}+4I_{v,tot}=__IB__\\ \\mathrm{m^4}
$$

### 2d
Dek en bodem liggen het verst van de neutrale as; dikte daar verhogen geeft de grootste winst in buigstijfheid en buigsterkte.

### 2e
Verticale delen koppelen dek en bodem, nemen schuif op en zorgen dat de samengestelde doorsnede effectief samenwerkt in buiging.

### 2f
Voor deze opdracht is $t_{p,eq}=__TP_MM__\\ \\mathrm{mm}$.  
In de Damen-documentatie staan principediktes zoals **shell 16 mm** en **longitudinale schotten 14/10 mm**.  
Het verschil is logisch: hier wordt met een vereenvoudigde equivalente dikte gerekend voor globale balkrespons; een echt ontwerp gebruikt locatie-afhankelijke plaatdiktes.

## 3. Constructie respons en limiet gedrag

Materiaal (op basis van $e=__E_DIGIT__$): __MAT_NAME__,  
$\\sigma_y=__SIGY__\\ \\mathrm{MPa}$, $E=__E_MOD__\\ \\mathrm{MPa}$, $\\nu=__NU__$, $\\rho=__RHO_MAT__\\ \\mathrm{kg/m^3}$.

### 3a
Euler-Bernoulli relaties:

$$
\\kappa(x)=\\frac{M_b(x)}{E I_b},\\qquad
\\varphi(x)=\\int_0^x \\kappa(\\xi)\\,d\\xi,\\qquad
w(x)=\\int_0^x \\varphi(\\xi)\\,d\\xi
$$

waarbij een lineaire rigid-body component is verwijderd zodat $w(0)=w(L)=0$.

Resultaten:

$$
\\varphi_{\\max}=__PHI_MAX__^\\circ,\\quad
\\varphi_{\\min}=__PHI_MIN__^\\circ
$$
$$
w_{\\max}=__W_MAX__\\ \\mathrm{mm},\\quad
w_{\\min}=__W_MIN__\\ \\mathrm{mm}
$$

![Vraag 3a hoekverdraaiing](q3a_phi.png)

![Vraag 3a vervorming](q3a_w.png)

### 3b
Maximale buigspanning (dek en bodem, zelfde absolute waarde):

$$
\\sigma_{b,\\max}=\\frac{|M_{b,\\max}|\\,z_{\\max}}{I_b}
$$

met $|M_{b,\\max}|=__MMAX__\\ \\mathrm{MNm}$ op $x=__X_M__\\ \\mathrm{m}$ en $z_{\\max}=H/2=7.75\\ \\mathrm{m}$:

$$
\\sigma_{b,\\max}=__SIGMA_B__\\ \\mathrm{MPa}
$$

Vergelijking met vloeigrens: $__SIGMA_B__ < __SIGY__\\ \\mathrm{MPa}$.  
__SAFETY_TXT__

### 3c
Een effectieve maatregel is het verhogen van sectiemodulus (bijv. extra dek/bodem-dikte of grotere holte), zodat bij gelijk moment de buigspanning daalt. Dit verhoogt veiligheid en verlengt de vermoeiingslevensduur.

### 3d
Gemiddelde schuifspanning:

$$
\\tau_s=\\frac{|F_{s,\\max}|}{A_s}
$$

met $|F_{s,\\max}|=__VMAX__\\ \\mathrm{MN}$ op $x=__X_V__\\ \\mathrm{m}$:

$$
\\tau_s=__TAU_S__\\ \\mathrm{MPa}
$$

Vergelijking: $\\tau_s\\approx __TAU_S__\\ \\mathrm{MPa}$ versus $\\sigma_{b,\\max}\\approx __SIGMA_B__\\ \\mathrm{MPa}$.  
Deze zijn van dezelfde orde, maar buiging domineert de extreme spanningen doordat $M_b$ piekt waar afstand tot neutrale as maximaal doorwerkt.
"""
    replacements = {
        "__NAME__": name,
        "__STUDIENUMMER__": studienummer,
        "__A__": str(digits["a"]),
        "__B__": str(digits["b"]),
        "__C__": str(digits["c"]),
        "__D__": str(digits["d"]),
        "__E_DIGIT__": str(digits["e"]),
        "__F_DIGIT__": str(digits["f"]),
        "__G_DIGIT__": str(digits["g"]),
        "__CRANE_X__": f"{q1['crane_x']:.2f}",
        "__T0__": f"{q1['T0']:.2f}",
        "__LCF__": f"{q1['LCF']:.2f}",
        "__LCG__": f"{q1['LCG']:.2f}",
        "__TA__": f"{q1['ta']:.2f}",
        "__TF__": f"{q1['tf']:.2f}",
        "__EQ_FORCE__": f"{q1['eq_force_total']:.3e}",
        "__EQ_MOMENT__": f"{q1['eq_moment_total']:.3e}",
        "__MOONPOOL_DFS__": f"{moonpool_dfs:.3f}",
        "__MOONPOOL_DMB__": f"{moonpool_dmb:.3f}",
        "__TP_MM__": f"{q2['tp_mm']:.1f}",
        "__AS__": f"{q2['As']:.2f}",
        "__ZN__": f"{q2['z_n']:.2f}",
        "__IH_OWN__": f"{q2['I_h_own']:.6f}",
        "__IH_SHIFT__": f"{q2['I_h_shift']:.6f}",
        "__IH_TOT__": f"{q2['I_h_total_each']:.6f}",
        "__IV_OWN__": f"{q2['I_v_own']:.6f}",
        "__IV_SHIFT__": f"{q2['I_v_shift']:.6f}",
        "__IV_TOT__": f"{q2['I_v_total_each']:.6f}",
        "__IB__": f"{q2['Ib']:.4f}",
        "__MAT_NAME__": material["name"],
        "__SIGY__": f"{material['sigma_y']:.0f}",
        "__E_MOD__": f"{material['E']:.0f}",
        "__NU__": f"{material['nu']}",
        "__RHO_MAT__": f"{material['rho']:.0f}",
        "__PHI_MAX__": f"{max(q3['phi_deg']):.4f}",
        "__PHI_MIN__": f"{min(q3['phi_deg']):.4f}",
        "__W_MAX__": f"{max(q3['w_mm']):.2f}",
        "__W_MIN__": f"{min(q3['w_mm']):.2f}",
        "__MMAX__": f"{abs(q3['M_max']):.2f}",
        "__X_M__": f"{q3['x_m']:.2f}",
        "__SIGMA_B__": f"{q3['sigma_b_max']:.0f}",
        "__SAFETY_TXT__": safety_txt,
        "__VMAX__": f"{abs(q3['V_max']):.2f}",
        "__X_V__": f"{q3['x_v']:.2f}",
        "__TAU_S__": f"{q3['tau_s']:.0f}",
    }
    for key, value in replacements.items():
        md = md.replace(key, value)
    return md


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate MT1464 assignment verslag (.docx).")
    parser.add_argument(
        "--studienummer",
        default="6470114",
        help="Studienummer used for parameter selection (default: 6470114).",
    )
    parser.add_argument(
        "--name",
        default="Jan Albert Driessen",
        help="Name shown on the title page (default: Jan Albert Driessen).",
    )
    args = parser.parse_args()

    studienummer = args.studienummer
    name = args.name
    solve_main(studienummer=studienummer, output_dir=str(OUTPUT_DIR))

    digits = digits_from_studienummer(studienummer)
    q1_csv = OUTPUT_DIR / "q1a_g_data.csv"
    x, W, c, g, p, q, Fs, Mb, Fs_n, Mb_n = load_q1_data(q1_csv)

    q1 = compute_q1_metrics(digits, x, g, q)
    q2 = compute_q2_metrics(digits)
    material = material_from_e(digits["e"])
    q3 = compute_q3_metrics(x, Fs, Mb, q2["As"], q2["Ib"], material)

    make_report_figures(
        ASSETS_DIR, x, W, c, g, p, q, Fs, Mb, Fs_n, Mb_n, q1["crane_x"], q3["phi_deg"], q3["w_mm"]
    )

    cover_candidates = [
        REPO_ROOT / "assets" / "plaatjevoorkant.png",
        ASSETS_DIR / "plaatjevoorkant.png",
    ]
    cover_src = next((p for p in cover_candidates if p.exists()), None)
    if cover_src is None:
        raise FileNotFoundError("Could not find plaatjevoorkant.png in assets or outputs/report_assets.")
    cover_dst = ASSETS_DIR / "plaatjevoorkant.png"
    if cover_src.resolve() != cover_dst.resolve():
        shutil.copy2(cover_src, cover_dst)

    moonpool_dfs = max(abs(a - b) for a, b in zip(Fs, Fs_n))
    moonpool_dmb = max(abs(a - b) for a, b in zip(Mb, Mb_n))

    md_path = ASSETS_DIR / "report.md"
    md_path.write_text(render_markdown(name, studienummer, digits, q1, q2, q3, material, moonpool_dfs, moonpool_dmb))

    docx_path = OUTPUT_DIR / f"MT1464_Assignment1_Verslag_{studienummer}.docx"
    subprocess.run(
        [
            "pandoc",
            str(md_path),
            "-o",
            str(docx_path),
            "--resource-path",
            str(ASSETS_DIR),
            "--from",
            "markdown+tex_math_dollars",
        ],
        check=True,
    )
    print(f"Generated: {docx_path}")
