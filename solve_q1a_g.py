from dataclasses import dataclass

RHO = 1025.0
G = 9.81
RHO_G_MN = RHO * G / 1e6
B = 31.0

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


def draw_multiplot_svg(filename, x, series, labels, colors):
    w, h = 1100, 1700
    margin_l, margin_r = 70, 20
    margin_t, margin_b = 30, 30
    n = len(series)
    panel_h = (h - margin_t - margin_b) / n

    def xmap(xv):
        return margin_l + (xv - x[0]) / (x[-1] - x[0]) * (w - margin_l - margin_r)

    with open(filename, "w") as f:
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


def main(g_digit=4):
    p = pick_params(g_digit)
    Lm, Bm, xm = moonpool(g_digit, p.L)
    n = 3001
    x = [i * p.L / (n - 1) for i in range(n)]

    seg_len = p.L / 6
    c_vals = [2.0, 3.0, p.c3, 2.5, 2.0, 1.0]

    crane_x = 2 * seg_len if p.crane_transition == "2-3" else 3 * seg_len
    spread = 0.2

    W, c, gnet = [], [], []
    beff = []
    for xi in x:
        s = int(xi / seg_len)
        if s > 5:
            s = 5
        ci = c_vals[s] + 1.5 * delta_rect(xi, crane_x, spread)
        wi = p.W_const + 2.0 * delta_rect(xi, crane_x, spread)
        gi = wi + ci
        in_moon = 1.0 if (xm - Lm / 2 <= xi <= xm + Lm / 2) else 0.0
        b = B - Bm * in_moon
        W.append(wi); c.append(ci); gnet.append(gi); beff.append(b)

    A_wp = trapz(beff, x)
    xbeff = [x[i] * beff[i] for i in range(n)]
    LCF = trapz(xbeff, x) / A_wp

    Wtot = trapz(gnet, x)
    xg = [x[i] * gnet[i] for i in range(n)]
    LCG = trapz(xg, x) / Wtot

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
    A_wp_n = trapz(beff_n, x)
    LCF_n = trapz([x[i] * beff_n[i] for i in range(n)], x) / A_wp_n
    I0n = trapz(beff_n, x)
    I2n = trapz([((xi - LCF_n) ** 2) * beff_n[i] for i, xi in enumerate(x)], x)
    T0n = Wtot / (RHO_G_MN * I0n)
    an = Wtot * (LCG - LCF_n) / (RHO_G_MN * I2n)
    Tn = [T0n + an * (xi - LCF_n) for xi in x]
    buoy_n = [RHO_G_MN * beff_n[i] * Tn[i] for i in range(n)]
    qn = [buoy_n[i] - gnet[i] for i in range(n)]
    Fs_n = cumulative_integral(qn, x)
    Mb_n = cumulative_integral(Fs_n, x)

    with open('q1a_g_data.csv', 'w') as f:
        f.write('x,W,c,g,p,q,Fs,Mb,Fs_nomoon,Mb_nomoon\n')
        for i in range(n):
            f.write(f"{x[i]:.6f},{W[i]:.6f},{c[i]:.6f},{gnet[i]:.6f},{buoy[i]:.6f},{q[i]:.6f},{Fs[i]:.6f},{Mb[i]:.6f},{Fs_n[i]:.6f},{Mb_n[i]:.6f}\n")

    draw_multiplot_svg(
        'q1a_g_plots.svg',
        x,
        [W, c, gnet, buoy, q, Fs, Mb],
        ['W(x) [MN/m]', 'c(x) [MN/m]', 'g(x) [MN/m]', 'p(x) [MN/m]', 'q(x) [MN/m]', 'Fs(x) [MN]', 'Mb(x) [MNm]'],
        ['#0055aa', '#aa5500', '#228833', '#9933cc', '#cc2222', '#006666', '#444444']
    )

    draw_multiplot_svg(
        'q1f_compare.svg',
        x,
        [Fs, Fs_n, Mb, Mb_n],
        ['Fs met moonpool [MN]', 'Fs zonder moonpool [MN]', 'Mb met moonpool [MNm]', 'Mb zonder moonpool [MNm]'],
        ['#0055aa', '#dd8800', '#228833', '#cc2222']
    )

    eq_force = trapz(q, x)
    eq_moment = trapz([(x[i] - LCF) * q[i] for i in range(n)], x)

    with open('answer_q1a_g.md', 'w') as f:
        f.write(f"# Uitwerking vraag 1a t/m 1g (aanname cijfer g = {g_digit})\n\n")
        f.write("## 1a\n")
        f.write("- W(x), c(x), g(x)=W(x)+c(x) zijn berekend per segment met expliciete puntlastbijdrage van de kraan (eigen gewicht in W, SWL in c).\n")
        f.write("- Data staat in `q1a_g_data.csv`; figuren in `q1a_g_plots.svg`.\n\n")
        f.write("## 1b\n")
        f.write(f"- Gemiddelde draft T = {T0:.2f} m.\n")
        f.write(f"- LCF = {LCF:.2f} m.\n")
        f.write(f"- LCG = {LCG:.2f} m.\n")
        f.write(f"- ta (achter) = {ta:.2f} m, tf (voor) = {tf:.2f} m.\n\n")
        f.write("## 1c\n")
        f.write("- Buoyancy verdeling: p(x)=rho*g*b_eff(x)*T(x), met b_eff=B-Bm in moonpoolgebied, anders B.\n\n")
        f.write("## 1d\n")
        f.write("- Resultante langsscheepse vlakwater belasting: q(x)=p(x)-g(x).\n")
        f.write(f"- Evenwichtcheck: integraal q dx = {eq_force:.3e} MN; integraal (x-LCF)q dx = {eq_moment:.3e} MNm (beide ~0).\n\n")
        f.write("## 1e\n")
        f.write("- Dwarskracht: Fs(x)=integral(q dx), buigend moment: Mb(x)=integral(Fs dx), met Fs(0)=0 en Mb(0)=0.\n")
        f.write(f"- Eindecontrole numeriek: Fs(L)={Fs[-1]:.3e} MN, Mb(L)={Mb[-1]:.3e} MNm.\n\n")
        f.write("## 1f\n")
        f.write("- Vergelijking met/zonder moonpool staat in `q1f_compare.svg`.\n")
        f.write("- Consequentie: lokaal minder drijfvermogen bij moonpool veroorzaakt extra neerwaartse netto belasting, zichtbaar als lokale verandering in Fs en sterkere kromming in Mb rond x=xm.\n\n")
        f.write("## 1g\n")
        f.write("- Controle via randvoorwaarden voor een vrij drijvende romp: Fs(0)=0, Mb(0)=0, Fs(L)=0, Mb(L)=0.\n")
        f.write("- Extra controle: numeriek afleiden moet opleveren dFs/dx=q en dMb/dx=Fs.\n")

    print('Generated: answer_q1a_g.md, q1a_g_data.csv, q1a_g_plots.svg, q1f_compare.svg')


if __name__ == '__main__':
    main(4)
