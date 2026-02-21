# Uitwerking vraag 2a t/m 2f

- Studienummer cijfers: a=6, b=4, c=7, d=0, e=1, f=1, g=4.
- Voor vraag 2 wordt f gebruikt: f=1 -> t_p,eq=7.5 mm (0.0075 m).

## 2a) Schuifoppervlak As
- Verticale delen leveren de grootste bijdrage aan schuifstijfheid, omdat de schuifspanning/afschuifstroom in een slanke doorsnede vooral via web-achtige verticale platen loopt; deck en bodem dragen relatief minder bij in dit vereenvoudigde model.
- Formule: As(tp)=5·H·(2tp)=10Htp.
- As = 1.16 m².

## 2b) Hoogte neutrale as zn
- Formule: z_n = (Σ A_i z_i)/(Σ A_i). Door verticale symmetrie van de doorsnede volgt direct z_n=H/2.
- z_n = 7.75 m boven de basis.

## 2c) Oppervlaktetraagheidsmoment Ib
- Formule per deel: I_b = Σ(I_{eigen,i} + A_i d_i²).
- Horizontale delen (dek + bodem, elk):
  - I_eigen = b t³/12 = 0.000070 m⁴
  - A d² = 55.642109 m⁴
  - I_totaal per plaat = 55.642179 m⁴
- Verticale delen (5 stuks, elk):
  - I_eigen = t h³/12 = 4.654844 m⁴
  - A d² = 0.000000 m⁴
  - I_totaal per plaat = 4.654844 m⁴
- Totaal I_b = 2·I_h + 5·I_v = 134.5586 m⁴.
- Grootste bijdrage: dek en bodem (grote afstand tot neutrale as -> grote A d²-term).

## 2d) Waarom vaak dikkere bodem/dek platen dan zijbeplating
- Bodem en dek liggen het verst van de neutrale as en dragen daardoor het meeste aan buigsterkte/stijfheid bij; diktevergroting daar is structureel het effectiefst.

## 2e) Functie verticale constructiedelen voor buigstijfheid
- Verticale delen koppelen dek en bodem, houden de doorsnede-vorm vast en leveren aanvullende I_b-bijdrage; vooral belangrijk voor schuifkracht-opname en voor het realiseren van samengestelde buigwerking van de hele doorsnede.

## 2f) Vergelijking met Damen Stan Pontoon tp
- In deze opdracht volgt t_p,eq = 7.5 mm uit tabel 4 (op basis van f).
- In praktijkdatasheets van pontons worden vaak verschillende plaatdiktes per locatie toegepast (bijv. dek/bodem zwaarder dan zijden), afhankelijk van sterkte, lokale belastingen, corrosiemarges, vermoeiing, en klasse-eisen.
- Dus vergelijkbaar in ordegrootte kan, maar exact gelijk hoeft niet omdat ontwerpdoelen en veiligheidsmarges verschillen.
