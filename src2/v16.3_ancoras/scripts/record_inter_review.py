"""Record the agent's visual review of all 213 maps, not a learned classifier.

Positions refer to the complete group index in indice_galerias.json. The maps
were inspected individually in contact sheets with common axes and colour
limits. Labels do not affect selection, calculations, or the Pearson matrix.
"""
import sys
sys.dont_write_bytecode=True
import json
from collections import Counter
from common import ROOT, sha, write_json
from pair_identity import pair_key
from run_full_inter import OUT
from run_full_intra import now

# Approximate wavelength bands observed in the contact sheets. These are
# annotations, not numerical thresholds or new spectral measurements.
HIGH = {
    "plug × slug": {2:"350–550 e perto de 1.000",8:"350–650 e perto de 1.000"},
    "plug × churn": {6:"350–600",7:"450–700",10:"400–650",16:"350–650",20:"350–650",22:"600–1.000, com contribuicao perto de 350"},
    "plug × wavy annular": {4:"400–650",7:"450–700",8:"500–1.000",16:"450–650",17:"350–600",19:"600–1.000",23:"600–1.000",26:"250–450"},
    "slug × churn": {1:"200–300 e 500–700",4:"450–1.000",11:"550–1.000",21:"350–550",22:"250–400 e 600–1.000",28:"550–1.000"},
    "slug × wavy annular": {3:"500–700",7:"550–850",9:"500–700",10:"650–1.000, com contribuicao perto de 300",16:"500–1.000",21:"350–700",23:"300–400 e 650–900",24:"300–450 e 600–1.000",25:"200–300 e 650–1.000"},
    "churn × wavy annular": {
        1:"600–1.000",3:"550–800",7:"400–1.000",10:"650–1.000, com contribuicoes perto de 250–450",
        12:"300–450",14:"300–450 e 650–900",15:"300–450 e 600–1.000",17:"250–400 e perto de 1.000",
        20:"200–300 e 400–650",21:"350–450, com contribuicoes perto de 200 e 700",22:"350–500 e 750–1.000",
        24:"650–1.000, com contribuicao perto de 350",26:"300–450 e perto de 1.000",28:"400–700",
        29:"650–950",30:"550–700",34:"550–850",36:"300–450 e 500–750",38:"500–900",40:"350–650",
        43:"400–550 e 700–1.000",44:"500–750",47:"200–300 e 550–800",48:"200–300 e 550–850",
        49:"350–550",52:"450–700",53:"350–700",54:"650–1.000, com contribuicao perto de 250",
        57:"450–650",62:"450–650",65:"350–500 e 750–950",66:"300–450 e 550–750",67:"180–550 e 750–1.000",
        68:"300–450 e 650–950",69:"300–450 e 600–1.000",71:"300–500",74:"250–550 e 650–900",
        75:"350–450 e 650–850",76:"300–550",77:"300–400",80:"350–550",83:"200–300 e 550–750",
        85:"350–550",88:"450–650",89:"450–700"}
}


def main():
    rows=json.loads((OUT/"metricas_wavelets.json").read_text())
    index=json.loads((OUT/"indice_galerias.json").read_text())
    lookup={r["pair_id"]:r for r in index}
    assert len(rows)==len(index)==213
    labels={}
    for row in rows:
        key=pair_key(row)
        evidence=lookup[key]
        group,position=evidence["grupo"],evidence["posicao"]
        assert group==row["padrao"] and group in HIGH
        bands=HIGH[group].get(position)
        reason=(f"Faixa(s) de coerencia elevada junto ao alinhamento, aproximadamente em lambda {bands} px. "
                "A semelhança localizada entre padroes diferentes desfavorece a separacao."
                if bands else "Coerencia predominantemente baixa ou moderada junto ao alinhamento; "
                "regioes altas restritas, nas maiores escalas ou afastadas do lag zero. Favoravel neste sentido local.")
        labels[key]=dict(run_x=row["run_x"],run_y=row["run_y"],padrao=group,
            perfil="comparavel",comprimento=2048,rotulo="desfavoravel" if bands else "favoravel",
            motivo_visual=reason,natureza="Orientacao visual do agente apos inspecao de todos os mapas; "
            "favoravel significa WCO predominantemente baixa entre padroes diferentes. "
            "Nao e um teste, limiar estatistico, taxa de acerto ou ausencia de toda coerencia localizada.",
            mapa_sha256=row["mapa_sha256"],evidencia=evidence["galeria"],
            evidencia_sha256=sha(ROOT/evidence["galeria"]),posicao_na_galeria=position,data_utc=now())
    for group,positions in HIGH.items():
        assert set(positions)<={r["posicao"] for r in index if r["grupo"]==group}
    write_json(OUT/"rotulos_visuais.json",labels)
    print(dict(Counter(v["rotulo"] for v in labels.values())))


if __name__=="__main__":
    main()
