from dataclasses import dataclass
import argparse
from pathlib import Path

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


def delta_rect(x, x0, width):
    return (1.0 / width) if abs(x - x0) <= width / 2 else 0.0


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


def draw_multiplot_svg(filename: Path, x, series, labels, colors):
    w, h = 1100, 1700
    margin_l, margin_r = 70, 20
    margin_t, margin_b = 30, 30
    n = len(series)
    panel_h = (h - margin_t - margin_b) / n

    def xmap(xv):
        return margin_l + (xv - x[0]) / (x[-1] - x[0]) * (w - margin_l - margin_r)

    with filename.open("w") as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">\n')
        f.write('<style>text{font-family:Arial;font-size:12px}.title{font-weight:bold}</style>\n')
        for k, y in enumerate(series):
            y0 = margin_t + k * panel_h
            ymin, ymax = min(y), max(y)
            if ymax - ymin < 1e-9:
                ymax = ymin + 1.0

            def ymap(yv):
                return y0 + panel_h - 20 - (yv - ymin) / (ymax - ymin) * (panel_h - 40)

            f.write(f'<rect x="{margin_l}" y="{y0+10}" width="{w-margin_l-margin_r}" height="{panel_h-20}" fill="none" stroke="#ccc"/>\n')
            pts = " ".join(f"{xmap(xi):.2f},{ymap(yi):.2f}" for xi, yi in zip(x, y))
            f.write(f'<polyline points="{pts}" fill="none" stroke="{colors[k%len(colors)]}" stroke-width="1.2"/>\n')
            f.write(f'<text x="{margin_l+5}" y="{y0+24}" class="title">{labels[k]}</text>\n')
            f.write(f'<text x="10" y="{y0+24}">min {ymin:.2f}</text>\n')
            f.write(f'<text x="10" y="{y0+40}">max {ymax:.2f}</text>\n')
        f.write('</svg>\n')


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
    spread = 0.2

    W, c, gnet, beff = [], [], [], []
    for xi in x:
        s = min(int(xi / seg_len), 5)
        ci = c_vals[s] + 1.5 * delta_rect(xi, crane_x, spread)
        wi = p.W_const + 2.0 * delta_rect(xi, crane_x, spread)
        gi = wi + ci
        in_moon = 1.0 if (xm - Lm / 2 <= xi <= xm + Lm / 2) else 0.0
        b = B - Bm * in_moon
        W.append(wi)
        c.append(ci)
        gnet.append(gi)
        beff.append(b)

    A_wp = trapz(beff, x)
    LCF = trapz([x[i] * beff[i] for i in range(n)], x) / A_wp
    Wtot = trapz(gnet, x)
    LCG = trapz([x[i] * gnet[i] for i in range(n)], x) / Wtot

    I0 = trapz(beff, x)
    I2 = trapz([((xi - LCF) ** 2) * beff[i] for i, xi in enumerate(x)], x)
    T0 = Wtot / (RHO_G_MN * I0)
    a = Wtot * (LCG - LCF) / (RHO_G_MN * I2)

    T = [T0 + a * (xi - LCF) for xi in x]
    ta = T0 + a * (0 - LCF)
    tf = T0 + a * (p.L - LCF)

    buoy = [RHO_G_MN * beff[i] * T[i] for i in range(n)]
    q = [buoy[i] - gnet[i] for i in range(n)]
    Fs = cumulative_integral(q, x)
    Mb = cumulative_integral(Fs, x)

    beff_n = [B] * n
    LCF_n = trapz([x[i] * beff_n[i] for i in range(n)], x) / trapz(beff_n, x)
    I2n = trapz([((xi - LCF_n) ** 2) * beff_n[i] for i, xi in enumerate(x)], x)
    T0n = Wtot / (RHO_G_MN * trapz(beff_n, x))
    an = Wtot * (LCG - LCF_n) / (RHO_G_MN * I2n)
    Tn = [T0n + an * (xi - LCF_n) for xi in x]
    buoy_n = [RHO_G_MN * beff_n[i] * Tn[i] for i in range(n)]
    qn = [buoy_n[i] - gnet[i] for i in range(n)]
    Fs_n = cumulative_integral(qn, x)
    Mb_n = cumulative_integral(Fs_n, x)

    draw_multiplot_svg(output_dir / 'q1a_g_plots.svg', x,
                       [W, c, gnet, buoy, q, Fs, Mb],
                       ['W(x) [MN/m]', 'c(x) [MN/m]', 'g(x) [MN/m]', 'p(x) [MN/m]', 'q(x) [MN/m]', 'Fs(x) [MN]', 'Mb(x) [MNm]'],
                       ['#0055aa', '#aa5500', '#228833', '#9933cc', '#cc2222', '#006666', '#444444'])

    draw_multiplot_svg(output_dir / 'q1f_compare.svg', x,
                       [Fs, Fs_n, Mb, Mb_n],
                       ['Fs met moonpool [MN]', 'Fs zonder moonpool [MN]', 'Mb met moonpool [MNm]', 'Mb zonder moonpool [MNm]'],
                       ['#0055aa', '#dd8800', '#228833', '#cc2222'])

    with (output_dir / 'q1a_g_data.csv').open('w') as f:
        f.write('x,W,c,g,p,q,Fs,Mb,Fs_nomoon,Mb_nomoon\n')
        for i in range(n):
            f.write(f"{x[i]:.6f},{W[i]:.6f},{c[i]:.6f},{gnet[i]:.6f},{buoy[i]:.6f},{q[i]:.6f},{Fs[i]:.6f},{Mb[i]:.6f},{Fs_n[i]:.6f},{Mb_n[i]:.6f}\n")

    eq_force = trapz(q, x)
    eq_moment = trapz([(x[i] - LCF) * q[i] for i in range(n)], x)
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
        f.write(f'- T={T0:.2f} m, LCF={LCF:.2f} m, LCG={LCG:.2f} m, ta={ta:.2f} m, tf={tf:.2f} m.\n')
        f.write(f'- Evenwicht: ∫qdx={eq_force:.3e} MN; ∫(x-LCF)qdx={eq_moment:.3e} MNm.\n')
        f.write(f'- Randvoorwaarden: Fs(0)={Fs[0]:.3e}, Mb(0)={Mb[0]:.3e}, Fs(L)={Fs[-1]:.3e}, Mb(L)={Mb[-1]:.3e}.\n')
        f.write(f'- Afgeleide checks: max|dFs/dx-q|={max_dFs:.3e}, max|dMb/dx-Fs|={max_dMb:.3e}.\n')
        f.write(f'- Moonpool effect: max|Fs_met-Fs_zonder|={max(abs(Fs[i]-Fs_n[i]) for i in range(n)):.3f} MN; '
                f'max|Mb_met-Mb_zonder|={max(abs(Mb[i]-Mb_n[i]) for i in range(n)):.3f} MNm.\n')

    with (output_dir / 'answer_q1a_g.md').open('w') as f:
        f.write(f"# Uitwerking vraag 1a t/m 1g (studienummer {digits['a']}{digits['b']}{digits['c']}{digits['d']}{digits['e']}{digits['f']}{digits['g']})\n\n")
        f.write(f"- Afgeleide cijfers: a={digits['a']}, b={digits['b']}, c={digits['c']}, d={digits['d']}, e={digits['e']}, f={digits['f']}, g={digits['g']}.\n")
        f.write("- Voor vraag 1a t/m 1g wordt alleen cijfer g gebruikt.\n\n")
        f.write(f"- Kernwaarden: T={T0:.2f} m, LCF={LCF:.2f} m, LCG={LCG:.2f} m, ta={ta:.2f} m, tf={tf:.2f} m.\n")
        f.write(f"- Evenwichtcheck: ∫qdx={eq_force:.3e} MN, ∫(x-LCF)qdx={eq_moment:.3e} MNm.\n")
        f.write(f"- Randvoorwaardencheck: Fs(L)={Fs[-1]:.3e} MN, Mb(L)={Mb[-1]:.3e} MNm.\n")


def solve_q2(digits, output_dir: Path):
    f_digit = digits['f']
    tp_mm = pick_tp_eq_mm(f_digit)
    tp = tp_mm / 1000.0

    t_h = 4 * tp
    t_v = 2 * tp

    As = 4 * H * t_v
    z_n = H / 2.0

    A_h = B * t_h
    I_h_own = B * t_h ** 3 / 12
    d_h = H / 2 - t_h / 2
    I_h_shift = A_h * d_h ** 2
    I_h_total_each = I_h_own + I_h_shift

    I_v_own = t_v * H ** 3 / 12
    I_v_shift = 0.0
    I_v_total_each = I_v_own

    Ib = 2 * I_h_total_each + 4 * I_v_total_each

    with (output_dir / 'answer_q2a_f.md').open('w') as f:
        f.write('# Uitwerking vraag 2a t/m 2f\n\n')
        f.write(f"- Studienummer cijfers: a={digits['a']}, b={digits['b']}, c={digits['c']}, d={digits['d']}, e={digits['e']}, f={digits['f']}, g={digits['g']}.\n")
        f.write(f"- Voor vraag 2 wordt f gebruikt: f={f_digit} -> t_p,eq={tp_mm:.1f} mm ({tp:.4f} m).\n\n")

        f.write('## 2a) Schuifoppervlak As\n')
        f.write('- Verticale delen leveren de grootste bijdrage aan schuifstijfheid, omdat de schuifspanning/afschuifstroom in een slanke doorsnede vooral via web-achtige verticale platen loopt; deck en bodem dragen relatief minder bij in dit vereenvoudigde model.\n')
        f.write('- Formule: As(tp)=4·H·(2tp)=8Htp.\n')
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
        f.write('- Verticale delen (4 stuks, elk):\n')
        f.write(f'  - I_eigen = t h³/12 = {I_v_own:.6f} m⁴\n')
        f.write(f'  - A d² = {I_v_shift:.6f} m⁴\n')
        f.write(f'  - I_totaal per plaat = {I_v_total_each:.6f} m⁴\n')
        f.write(f'- Totaal I_b = 2·I_h + 4·I_v = {Ib:.4f} m⁴.\n')
        f.write('- Grootste bijdrage: dek en bodem (grote afstand tot neutrale as -> grote A d²-term).\n\n')

        f.write('## 2d) Waarom vaak dikkere bodem/dek platen dan zijbeplating\n')
        f.write('- Bodem en dek liggen het verst van de neutrale as en dragen daardoor het meeste aan buigsterkte/stijfheid bij; diktevergroting daar is structureel het effectiefst.\n\n')

        f.write('## 2e) Functie verticale constructiedelen voor buigstijfheid\n')
        f.write('- Verticale delen koppelen dek en bodem, houden de doorsnede-vorm vast en leveren aanvullende I_b-bijdrage; vooral belangrijk voor schuifkracht-opname en voor het realiseren van samengestelde buigwerking van de hele doorsnede.\n\n')

        f.write('## 2f) Vergelijking met Damen Stan Pontoon tp\n')
        f.write('- In deze opdracht volgt t_p,eq = 7.5 mm uit tabel 4 (op basis van f).\n')
        f.write('- In praktijkdatasheets van pontons worden vaak verschillende plaatdiktes per locatie toegepast (bijv. dek/bodem zwaarder dan zijden), afhankelijk van sterkte, lokale belastingen, corrosiemarges, vermoeiing, en klasse-eisen.\n')
        f.write('- Dus vergelijkbaar in ordegrootte kan, maar exact gelijk hoeft niet omdat ontwerpdoelen en veiligheidsmarges verschillen.\n')


def main(studienummer="6470114", output_dir="outputs"):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    digits = digits_from_studienummer(studienummer)
    solve_q1(digits, out)
    solve_q2(digits, out)
    print(f'Generated outputs in: {out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Los vraag 1a t/m 1g en 2a t/m 2f op voor een gegeven studienummer.')
    parser.add_argument('--studienummer', default='6470114', help='Studienummer (standaard: 6470114)')
    parser.add_argument('--output-dir', default='outputs', help='Directory voor gegenereerde bestanden')
    args = parser.parse_args()
    main(args.studienummer, args.output_dir)
