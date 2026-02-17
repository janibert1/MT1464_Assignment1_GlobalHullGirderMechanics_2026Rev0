# Uitwerking vraag 3a t/m 3d

- Materiaal (op basis van e=1): staal S235J2, E=205000 MPa, sigma_y=235 MPa.
- Doorsnede: Ib=129.9037 m⁴, As=0.93 m², zn=7.75 m.

## 3a) Hoekverdraaiing en vervorming
- Gebruikte relaties: kappa=M/(EI), theta=∫kappa dx, w=∫theta dx met referentie BCs theta(0)=0 en w(0)=0.
- max|theta| = 0.1001 deg.
- max|w| = 117.62 mm.
- Plotbestand: `q3a_plots.svg`.

## 3b) Buigspanning bij grootste buigend moment
- Grootste |Mb| = 614.19 MNm op x=68.81 m.
- sigma_dek = 37 MPa, sigma_bodem = -37 MPa.
- Vergelijking met vloeigrens: sigma_max=37 MPa vs sigma_y=235 MPa -> ELASTISCH.

## 3c) Maatregel voor hogere veiligheid/duurzaamheid
- Verhoog sectiestijfheid/sterkte waar het het meest effectief is (bijv. lokaal extra dek/bodem dikte of langsscheepse verstijvers), zodat spanningen dalen en levensduur toeneemt.

## 3d) Gemiddelde schuifspanning
- Formule: tau_avg = |Fs|/As.
- Grootste |Fs| = 17.79 MN -> tau_avg = 19 MPa.
- Vergelijking met buigspanning: tau_avg/sigma_max = 0.522; schuifspanning is zelfde orde dan buigspanning binnen dit balkmodel.
