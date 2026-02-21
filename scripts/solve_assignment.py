from dataclasses import dataclass
import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO = 1025.0
G = 9.81
RHO_G_MN = RHO * G / 1e6
B = 31.0
H = 15.5


@dataclass
class Params:
    g_digit: int
    L: float
    W_const: float
    c3: float
    crane_transition: str


def pick_params(g_digit: int) -> Params:
    if g_digit in [0, 1]:
        return Params(g_digit, 125.0, 0.55, 3.0, "2-3")
    if g_digit in [2, 3, 4]:
        return Params(g_digit, 135.0, 0.55, 2.5, "2-3")
    if g_digit in [5, 6, 7]:
        return Params(g_digit, 145.0, 0.75, 4.0, "3-4")
    return Params(g_digit, 155.0, 0.75, 3.5, "3-4")


def pick_tp_eq_mm(f_digit: int) -> float:
    return 7.5 if f_digit <= 4 else 8.5


def moonpool(g_digit: int, L: float):
    if g_digit in [0, 1]:
        xm = 72.92
    elif g_digit in [2, 3, 4]:
        xm = 78.75
    elif g_digit in [5, 6, 7]:
        xm = 108.75
    else:
        xm = 116.25
    Lm = L / 20.0
    Bm = B / 5.0
    return Lm, Bm, xm


def trapz(y, x):
    s = 0.0
    for i in range(1, len(x)):
        s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    return s


def cumulative_integral(y, x):
    out = [0.0] * len(x)
    for i in range(1, len(x)):
        out[i] = out[i - 1] + 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    return out


def draw_multiplot(filename: Path, x, series, labels, colors):
    n = len(series)
    fig, axes = plt.subplots(n, 1, figsize=(11, max(2.5 * n, 8)), sharex=True)
    if n == 1:
        axes = [axes]

    for k, (ax, y) in enumerate(zip(axes, series)):
        color = colors[k % len(colors)]
        ax.plot(x, y, color=color, linewidth=1.2)
        ax.fill_between(x, 0.0, y, color=color, alpha=0.18)
        ymin, ymax = min(y), max(y)
        if abs(ymax - ymin) < 1e-12:
            ymin -= 0.5
            ymax += 0.5
        pad = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_ylabel(labels[k], fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("x [m]")
    fig.tight_layout()
    fig.savefig(filename, dpi=180)
    plt.close(fig)


def draw_compare_overlay(filename: Path, x, Fs, Fs_n, Mb, Mb_n):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    panels = [
        (axes[0], [(Fs, "Fs met moonpool [MN]", "#0055aa"), (Fs_n, "Fs zonder moonpool [MN]", "#dd8800")], "Fs [MN]"),
        (axes[1], [(Mb, "Mb met moonpool [MNm]", "#228833"), (Mb_n, "Mb zonder moonpool [MNm]", "#cc2222")], "Mb [MNm]"),
    ]

    for ax, curves, ylabel in panels:
        y_all = []
        for y, label, color in curves:
            ax.plot(x, y, color=color, linewidth=1.5, label=label)
            ax.fill_between(x, 0.0, y, color=color, alpha=0.16)
            y_all.extend(y)
        ymin, ymax = min(y_all), max(y_all)
        if abs(ymax - ymin) < 1e-12:
            ymin -= 0.5
            ymax += 0.5
        pad = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)

    axes[-1].set_xlabel("x [m]")
    fig.tight_layout()
    fig.savefig(filename, dpi=180)
    plt.close(fig)




def apply_point_moment_jump(moment_curve, x, x0, delta_m):
    """Add a step (jump) in bending moment at x0 representing a concentrated moment."""
    out = []
    for xi, mi in zip(x, moment_curve):
        out.append(mi + (delta_m if xi >= x0 else 0.0))
    return out


def apply_point_force_jump(force_curve, x, x0, delta_f):
    """Add a step (jump) in shear force at x0 representing a concentrated vertical force."""
    out = []
    for xi, fi in zip(x, force_curve):
        out.append(fi + (delta_f if xi >= x0 else 0.0))
    return out


def digits_from_studienummer(studienummer: str):
    s = ''.join(ch for ch in str(studienummer) if ch.isdigit())
    if len(s) < 7:
        raise ValueError("Studienummer moet minstens 7 cijfers hebben.")
    vals = [int(ch) for ch in s[-7:]]
    return {'a': vals[0], 'b': vals[1], 'c': vals[2], 'd': vals[3], 'e': vals[4], 'f': vals[5], 'g': vals[6]}


def solve_q1(digits, output_dir: Path):
    g_digit = digits['g']
    p = pick_params(g_digit)
    Lm, Bm, xm = moonpool(g_digit, p.L)
    n = 3001
    x = [i * p.L / (n - 1) for i in range(n)]

    seg_len = p.L / 6
    c_vals = [2.0, 3.0, p.c3, 2.5, 2.0, 1.0]
    crane_x = 2 * seg_len if p.crane_transition == "2-3" else 3 * seg_len
    crane_deadweight = 2.0  # [MN] concentrated crane self-weight
    crane_swl = 1.5         # [MN] concentrated lifted load (SWL)
    crane_total_point_load = crane_deadweight + crane_swl

    W, c, gnet, beff = [], [], [], []
    for xi in x:
        s = min(int(xi / seg_len), 5)
        # Distributed baseline loads only; crane actions are handled as concentrated terms.
        ci = c_vals[s]
        wi = p.W_const
        gi = wi + ci
        in_moon = 1.0 if (xm - Lm / 2 <= xi <= xm + Lm / 2) else 0.0
        b = B - Bm * in_moon
        W.append(wi)
        c.append(ci)
        gnet.append(gi)
        beff.append(b)

    A_wp = trapz(beff, x)
    LCF = trapz([x[i] * beff[i] for i in range(n)], x) / A_wp
    Wtot_distributed = trapz(gnet, x)
    Wtot = Wtot_distributed + crane_total_point_load
    LCG = (trapz([x[i] * gnet[i] for i in range(n)], x) + crane_x * crane_total_point_load) / Wtot

    I0 = trapz(beff, x)
    I2 = trapz([((xi - LCF) ** 2) * beff[i] for i, xi in enumerate(x)], x)
    T0 = Wtot / (RHO_G_MN * I0)
    a = Wtot * (LCG - LCF) / (RHO_G_MN * I2)

    T = [T0 + a * (xi - LCF) for xi in x]
    ta = a * (0 - LCF)
    tf = a * (p.L - LCF)

    buoy = [RHO_G_MN * beff[i] * T[i] for i in range(n)]
    # q excludes Dirac point loads; concentrated loads are applied as jumps in Fs.
    q = [buoy[i] - gnet[i] for i in range(n)]
    Fs = cumulative_integral(q, x)
    Fs = apply_point_force_jump(Fs, x, crane_x, -crane_deadweight)
    Fs = apply_point_force_jump(Fs, x, crane_x, -crane_swl)
    Mb = cumulative_integral(Fs, x)

    # Crane-induced concentrated moment at crane location (requested jump in Mb).
    crane_reach = 30.0  # [m], from assignment text
    crane_moment = crane_swl * crane_reach  # [MNm]
    Mb = apply_point_moment_jump(Mb, x, crane_x, crane_moment)

    beff_n = [B] * n
    LCF_n = trapz([x[i] * beff_n[i] for i in range(n)], x) / trapz(beff_n, x)
    I2n = trapz([((xi - LCF_n) ** 2) * beff_n[i] for i, xi in enumerate(x)], x)
    T0n = Wtot / (RHO_G_MN * trapz(beff_n, x))
    an = Wtot * (LCG - LCF_n) / (RHO_G_MN * I2n)
    Tn = [T0n + an * (xi - LCF_n) for xi in x]
    buoy_n = [RHO_G_MN * beff_n[i] * Tn[i] for i in range(n)]
    qn = [buoy_n[i] - gnet[i] for i in range(n)]
    Fs_n = cumulative_integral(qn, x)
    Fs_n = apply_point_force_jump(Fs_n, x, crane_x, -crane_deadweight)
    Fs_n = apply_point_force_jump(Fs_n, x, crane_x, -crane_swl)
    Mb_n = cumulative_integral(Fs_n, x)
    Mb_n = apply_point_moment_jump(Mb_n, x, crane_x, crane_moment)

    W_signed = [-wi for wi in W]
    c_signed = [-ci for ci in c]
    g_signed = [-gi for gi in gnet]

    draw_multiplot(output_dir / 'q1a_g_plots.png', x,
                       [W_signed, c_signed, g_signed, buoy, q, Fs, Mb],
                       ['W(x) [MN/m]', 'c(x) [MN/m]', 'g(x) [MN/m]', 'p(x) [MN/m]', 'q(x) [MN/m]', 'Fs(x) [MN]', 'Mb(x) [MNm]'],
                       ['#0055aa', '#aa5500', '#228833', '#9933cc', '#cc2222', '#006666', '#444444'])

    draw_compare_overlay(output_dir / 'q1f_compare.png', x, Fs, Fs_n, Mb, Mb_n)

    with (output_dir / 'q1a_g_data.csv').open('w') as f:
        f.write('x,W,c,g,p,q,Fs,Mb,Fs_nomoon,Mb_nomoon\n')
        for i in range(n):
            f.write(f"{x[i]:.6f},{W_signed[i]:.6f},{c_signed[i]:.6f},{g_signed[i]:.6f},{buoy[i]:.6f},{q[i]:.6f},{Fs[i]:.6f},{Mb[i]:.6f},{Fs_n[i]:.6f},{Mb_n[i]:.6f}\n")

    eq_force_distributed = trapz(q, x)
    eq_force_total = eq_force_distributed - crane_total_point_load
    eq_moment_distributed = trapz([(x[i] - LCF) * q[i] for i in range(n)], x)
    eq_moment_total = eq_moment_distributed - crane_total_point_load * (crane_x - LCF)
    disc = [seg_len, 2 * seg_len, 3 * seg_len, 4 * seg_len, 5 * seg_len, crane_x, xm - Lm / 2, xm + Lm / 2]

    def smooth_idx(i):
        xi = x[i]
        return all(abs(xi - d) > 0.6 for d in disc)

    dFs_err = [abs(((Fs[i + 1] - Fs[i]) / (x[i + 1] - x[i])) - q[i]) for i in range(n - 1) if smooth_idx(i)]
    dMb_err = [abs(((Mb[i + 1] - Mb[i]) / (x[i + 1] - x[i])) - Fs[i]) for i in range(n - 1) if smooth_idx(i)]
    max_dFs = max(dFs_err) if dFs_err else 0.0
    max_dMb = max(dMb_err) if dMb_err else 0.0

    with (output_dir / 'check_q1a_f.md').open('w') as f:
        f.write('# Check vraag 1a t/m 1f\n\n')
        f.write(f'- g={g_digit}, L={p.L:.2f} m, W_const={p.W_const:.2f} MN/m, c3={p.c3:.2f} MN/m.\n')
        f.write(f'- T={T0:.2f} m, LCF={LCF:.2f} m, LCG={LCG:.2f} m, ta(trim)={ta:.2f} m, tf(trim)={tf:.2f} m.\n')
        f.write('- Notatie: t_a = voor; t_f = achter.\n')
        f.write(f'- Kraan puntlasten op x={crane_x:.2f} m: eigengewicht={crane_deadweight:.2f} MN, SWL={crane_swl:.2f} MN.\n')
        f.write(f'- Evenwicht (incl. puntlasten): ∫qdx-P={eq_force_total:.3e} MN; ∫(x-LCF)qdx-P(xc-LCF)={eq_moment_total:.3e} MNm.\n')
        f.write(f'- Randvoorwaarden: Fs(0)={Fs[0]:.3e}, Mb(0)={Mb[0]:.3e}, Fs(L)={Fs[-1]:.3e}, Mb(L)={Mb[-1]:.3e}.\n')
        f.write(f'- Afgeleide checks: max|dFs/dx-q|={max_dFs:.3e}, max|dMb/dx-Fs|={max_dMb:.3e}.\n')
        f.write(f'- Moonpool effect: max|Fs_met-Fs_zonder|={max(abs(Fs[i]-Fs_n[i]) for i in range(n)):.3f} MN; '
                f'max|Mb_met-Mb_zonder|={max(abs(Mb[i]-Mb_n[i]) for i in range(n)):.3f} MNm.\n')
        f.write(f'- Kraanmoment-jump in Mb op x={crane_x:.2f} m: +{crane_moment:.2f} MNm.\n')

    with (output_dir / 'answer_q1a_g.md').open('w') as f:
        f.write(f"# Uitwerking vraag 1a t/m 1g (studienummer {digits['a']}{digits['b']}{digits['c']}{digits['d']}{digits['e']}{digits['f']}{digits['g']})\n\n")
        f.write(f"- Afgeleide cijfers: a={digits['a']}, b={digits['b']}, c={digits['c']}, d={digits['d']}, e={digits['e']}, f={digits['f']}, g={digits['g']}.\n")
        f.write("- Voor vraag 1a t/m 1g wordt alleen cijfer g gebruikt.\n\n")
        f.write(f"- Kernwaarden: T={T0:.2f} m, LCF={LCF:.2f} m, LCG={LCG:.2f} m, ta(trim)={ta:.2f} m, tf(trim)={tf:.2f} m.\n")
        f.write("- Notatie: t_a = voor; t_f = achter.\n")
        f.write(f"- Kraan gemodelleerd als puntlasten op x={crane_x:.2f} m: 2.00 MN (eigengewicht in W) + 1.50 MN (SWL in c).\n")
        f.write(f"- Evenwichtcheck (incl. puntlasten): ∫qdx-P={eq_force_total:.3e} MN, ∫(x-LCF)qdx-P(xc-LCF)={eq_moment_total:.3e} MNm.\n")
        f.write(f"- Randvoorwaardencheck: Fs(L)={Fs[-1]:.3e} MN, Mb(L)={Mb[-1]:.3e} MNm.\n")
        f.write(f"- Toegepaste kraan-geïnduceerde moment-jump bij x={crane_x:.2f} m: +{crane_moment:.2f} MNm.\n")


def solve_q2(digits, output_dir: Path):
    f_digit = digits['f']
    tp_mm = pick_tp_eq_mm(f_digit)
    tp = tp_mm / 1000.0

    t_h = 4 * tp
    t_v = 2 * tp

    n_vertical = 5
    As = n_vertical * H * t_v
    z_n = H / 2.0

    A_h = B * t_h
    I_h_own = B * t_h ** 3 / 12
    d_h = H / 2 - t_h / 2
    I_h_shift = A_h * d_h ** 2
    I_h_total_each = I_h_own + I_h_shift

    I_v_own = t_v * H ** 3 / 12
    I_v_shift = 0.0
    I_v_total_each = I_v_own

    Ib = 2 * I_h_total_each + n_vertical * I_v_total_each

    with (output_dir / 'answer_q2a_f.md').open('w') as f:
        f.write('# Uitwerking vraag 2a t/m 2f\n\n')
        f.write(f"- Studienummer cijfers: a={digits['a']}, b={digits['b']}, c={digits['c']}, d={digits['d']}, e={digits['e']}, f={digits['f']}, g={digits['g']}.\n")
        f.write(f"- Voor vraag 2 wordt f gebruikt: f={f_digit} -> t_p,eq={tp_mm:.1f} mm ({tp:.4f} m).\n\n")

        f.write('## 2a) Schuifoppervlak As\n')
        f.write('- Verticale delen leveren de grootste bijdrage aan schuifstijfheid, omdat de schuifspanning/afschuifstroom in een slanke doorsnede vooral via web-achtige verticale platen loopt; deck en bodem dragen relatief minder bij in dit vereenvoudigde model.\n')
        f.write(f'- Formule: As(tp)={n_vertical}·H·(2tp)={2*n_vertical}Htp.\n')
        f.write(f'- As = {As:.2f} m².\n\n')

        f.write('## 2b) Hoogte neutrale as zn\n')
        f.write('- Formule: z_n = (Σ A_i z_i)/(Σ A_i). Door verticale symmetrie van de doorsnede volgt direct z_n=H/2.\n')
        f.write(f'- z_n = {z_n:.2f} m boven de basis.\n\n')

        f.write('## 2c) Oppervlaktetraagheidsmoment Ib\n')
        f.write('- Formule per deel: I_b = Σ(I_{eigen,i} + A_i d_i²).\n')
        f.write('- Horizontale delen (dek + bodem, elk):\n')
        f.write(f'  - I_eigen = b t³/12 = {I_h_own:.6f} m⁴\n')
        f.write(f'  - A d² = {I_h_shift:.6f} m⁴\n')
        f.write(f'  - I_totaal per plaat = {I_h_total_each:.6f} m⁴\n')
        f.write(f'- Verticale delen ({n_vertical} stuks, elk):\n')
        f.write(f'  - I_eigen = t h³/12 = {I_v_own:.6f} m⁴\n')
        f.write(f'  - A d² = {I_v_shift:.6f} m⁴\n')
        f.write(f'  - I_totaal per plaat = {I_v_total_each:.6f} m⁴\n')
        f.write(f'- Totaal I_b = 2·I_h + {n_vertical}·I_v = {Ib:.4f} m⁴.\n')
        f.write('- Grootste bijdrage: dek en bodem (grote afstand tot neutrale as -> grote A d²-term).\n\n')

        f.write('## 2d) Waarom vaak dikkere bodem/dek platen dan zijbeplating\n')
        f.write('- Bodem en dek liggen het verst van de neutrale as en dragen daardoor het meeste aan buigsterkte/stijfheid bij; diktevergroting daar is structureel het effectiefst.\n\n')

        f.write('## 2e) Functie verticale constructiedelen voor buigstijfheid\n')
        f.write('- Verticale delen koppelen dek en bodem, houden de doorsnede-vorm vast en leveren aanvullende I_b-bijdrage; vooral belangrijk voor schuifkracht-opname en voor het realiseren van samengestelde buigwerking van de hele doorsnede.\n\n')

        f.write('## 2f) Vergelijking met Damen Stan Pontoon tp\n')
        f.write('- In deze opdracht volgt t_p,eq = 7.5 mm uit tabel 4 (op basis van f).\n')
        f.write('- In praktijkdatasheets van pontons worden vaak verschillende plaatdiktes per locatie toegepast (bijv. dek/bodem zwaarder dan zijden), afhankelijk van sterkte, lokale belastingen, corrosiemarges, vermoeiing, en klasse-eisen.\n')
        f.write('- Dus vergelijkbaar in ordegrootte kan, maar exact gelijk hoeft niet omdat ontwerpdoelen en veiligheidsmarges verschillen.\n')


def resolve_output_dir(output_dir: str) -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    out = Path(output_dir)
    if not out.is_absolute():
        out = repo_root / out
    return out


def main(studienummer="6470114", output_dir="outputs"):
    out = resolve_output_dir(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    digits = digits_from_studienummer(studienummer)
    solve_q1(digits, out)
    solve_q2(digits, out)
    print(f'Generated outputs in: {out}')


if __name__ == '__main__':
    repo_root = Path(__file__).resolve().parents[1]
    default_output = "output" if (repo_root / "output").exists() else "outputs"
    parser = argparse.ArgumentParser(description='Los vraag 1a t/m 1g en 2a t/m 2f op voor een gegeven studienummer.')
    parser.add_argument('--studienummer', default='6470114', help='Studienummer (standaard: 6470114)')
    parser.add_argument('--output-dir', default=default_output, help='Directory voor gegenereerde bestanden (relatief t.o.v. repo root)')
    args = parser.parse_args()
    main(args.studienummer, args.output_dir)
