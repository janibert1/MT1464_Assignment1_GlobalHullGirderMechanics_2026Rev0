# Uitwerking vraag 1a t/m 1g (studienummer 6470114)

- Afgeleide cijfers: a=6, b=4, c=7, d=0, e=1, f=1, g=4.
- Voor vraag 1a t/m 1g wordt in de opgave alleen cijfer g gebruikt, dus g = 4.

## 1a
- W(x), c(x), g(x)=W(x)+c(x) zijn berekend per segment met expliciete puntlastbijdrage van de kraan (eigen gewicht in W, SWL in c).
- Data staat in `q1a_g_data.csv`; figuren in `q1a_g_plots.svg`.

## 1b
- Gemiddelde draft T = 8.90 m.
- LCF = 67.39 m.
- LCG = 61.79 m.
- ta (achter) = 11.09 m, tf (voor) = 6.70 m.

## 1c
- Buoyancy verdeling: p(x)=rho*g*b_eff(x)*T(x), met b_eff=B-Bm in moonpoolgebied, anders B.

## 1d
- Resultante langsscheepse vlakwater belasting: q(x)=p(x)-g(x).
- Evenwichtcheck: integraal q dx = 1.040e-11 MN; integraal (x-LCF)q dx = -5.164e-10 MNm (beide ~0).

## 1e
- Dwarskracht: Fs(x)=integral(q dx), buigend moment: Mb(x)=integral(Fs dx), met Fs(0)=0 en Mb(0)=0.
- Eindecontrole numeriek: Fs(L)=1.040e-11 MN, Mb(L)=-1.854e-04 MNm.

## 1f
- Vergelijking met/zonder moonpool staat in `q1f_compare.svg`.
- Consequentie: lokaal minder drijfvermogen bij moonpool veroorzaakt extra neerwaartse netto belasting, zichtbaar als lokale verandering in Fs en sterkere kromming in Mb rond x=xm.

## 1g
- Controle via randvoorwaarden voor een vrij drijvende romp: Fs(0)=0, Mb(0)=0, Fs(L)=0, Mb(L)=0.
- Extra controle: numeriek afleiden moet opleveren dFs/dx=q en dMb/dx=Fs.
