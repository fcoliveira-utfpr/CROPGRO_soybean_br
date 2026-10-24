"""Extrai GMR, tipo de crescimento e ciclo (dias) por cultivar dos catálogos de cultivares de soja da Embrapa.

O ciclo vem como "Ciclo médio: N dias" ou numa tabela por faixa de altitude ("até 500 m   114 – 120").
Cada catálogo é associado à(s) macrorregião(ões)/estado(s) que cobre (tabela CATALOGOS).
Saída: dados/literatura/catalogos_embrapa_ciclo.csv (uma linha por cultivar × catálogo × faixa de altitude).
"""
import re
import subprocess
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_PDF = RAIZ / "dados" / "literatura" / "catalogos_embrapa"
SAIDA = RAIZ / "dados" / "literatura" / "catalogos_embrapa_ciclo.csv"

# arquivo → (ano da edição, região de indicação)
CATALOGOS = {
    "FolhetoSoja.pdf": (2010, "MA/TO (Imperatriz)"),
    "FolhetoSoja201011MGfinalcapa.pdf": (2010, "MG"),
    "CAT.22011.pdf": (2011, "MT"),
    "SojaBA201112VE.pdf": (2011, "BA"),
    "catalogoSojaGO2012FINAL.pdf": (2012, "GO"),
    "CATALOGOSojaCentrosul2013.pdf": (2013, "Centro-Sul (MR 1–3)"),
    "CatalogoSoja2013MTRO.pdf": (2013, "MT/RO"),
    "cataologosoja.5.2015.pdf": (2015, "Centro-Sul (MR 1–3)"),
    "Cultivaressojanortenordeste2016OL1.pdf": (2016, "Norte/Nordeste (MR 4–5)"),
    "Catalogo4LivretoSojaFM2017OL.pdf": (2017, "Centro-Sul (MR 1–3)"),
    "catalogoSojaGO2018OL.pdf": (2018, "GO"),
    "Catalogo4Soja2019OL.pdf": (2019, "Brasil"),
    "Cultivaressojaroraima.pdf": (2019, "RR"),
}

RE_GMR = re.compile(r"maturidade\s+relativa\s*:?\s*(\d{1,2}[.,]\d)", re.I)
RE_NOME = re.compile(r"\b(BRS\s?[0-9]{3,4}[A-Z]{0,4}(?:\s?IPRO|\s?RR|\s?IPRO\s?STS|\s?STS)?|BRS\s+[A-Z][a-zé]+(?:\s[A-Z][a-zé]+)?)")
RE_CICLO_MEDIO = re.compile(r"ciclo\s+m[ée]dio\s*:?\s*(?:de\s*)?(\d{2,3})(?:\s*[-–a]\s*(\d{2,3}))?\s*dias", re.I)
RE_ALTITUDE = re.compile(r"(at[ée]\s*\d{3}\s*m|\d{3}\s*a\s*\d{3,4}\s*m|acima\s+de\s+\d{3,4}\s*m)\s+(\d{2,3})\s*[-–]\s*(\d{2,3})", re.I)
RE_HABITO = re.compile(r"(?:tipo|h[áa]bito)\s+de\s+crescimento\s*:?\s*(semi[- ]?determinado|indeterminado|determinado)", re.I)


def texto(pdf):
    return subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True).stdout


def extrair(pdf, janela_antes=25, janela_depois=45):
    linhas = texto(pdf).splitlines()
    registros = []
    for i, l in enumerate(linhas):
        m = RE_GMR.search(l)
        if not m:
            continue
        bloco_antes = linhas[max(0, i - janela_antes):i + 1]
        bloco = linhas[max(0, i - 5):i + janela_depois]
        nomes = [n for b in reversed(bloco_antes) for n in RE_NOME.findall(b)]
        nome = nomes[0].strip() if nomes else None
        hab = next((h.group(1).lower() for b in bloco for h in [RE_HABITO.search(b)] if h), None)
        base = dict(arquivo=pdf.name, cultivar=nome, gmr=float(m.group(1).replace(",", ".")), habito=hab, linha=i)
        alt = [(a.group(1), int(a.group(2)), int(a.group(3))) for b in bloco for a in RE_ALTITUDE.finditer(b)]
        cm = next((c for b in bloco for c in [RE_CICLO_MEDIO.search(b)] if c), None)
        if alt:
            vistos = set()
            for faixa, a, b in alt:
                chave = re.sub(r"\s+", " ", faixa.lower())
                if chave in vistos:
                    continue
                vistos.add(chave)
                registros.append({**base, "altitude": chave, "ciclo_min": a, "ciclo_max": b})
        elif cm:
            a = int(cm.group(1)); b = int(cm.group(2)) if cm.group(2) else a
            registros.append({**base, "altitude": None, "ciclo_min": a, "ciclo_max": b})
        else:
            registros.append({**base, "altitude": None, "ciclo_min": None, "ciclo_max": None})
    return registros


def compilar():
    regs = []
    for arq, (ano, regiao) in CATALOGOS.items():
        pdf = DIR_PDF / arq
        if pdf.exists():
            regs += [{**r, "ano_catalogo": ano, "regiao": regiao} for r in extrair(pdf)]
    d = pd.DataFrame(regs)
    d["ciclo_medio"] = (d.ciclo_min + d.ciclo_max) / 2
    return d


if __name__ == "__main__":
    d = compilar()
    d.to_csv(SAIDA, index=False)
    print(d.groupby("arquivo").agg(registros=("gmr", "size"), com_ciclo=("ciclo_medio", "count"), com_nome=("cultivar", "count")))
