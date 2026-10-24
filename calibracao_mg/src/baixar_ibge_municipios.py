"""Baixa a produção de soja por município (IBGE, PAM, tabela 1612 do SIDRA), 1980–2024.

Variáveis: 109 área plantada (ha), 216 área colhida (ha), 214 quantidade produzida (t), 112 rendimento médio (kg/ha).
Uma consulta por ano (limite de ~50 mil valores da API). Ano do IBGE = ano da colheita (2024 = safra 2023/24).
Saída: dados/brutos/ibge/sidra_1612_soja_municipios_1980_2024.csv
Uso: .venv/bin/python src/baixar_ibge_municipios.py
"""
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
SAIDA = RAIZ / "dados" / "brutos" / "ibge" / "sidra_1612_soja_municipios_1980_2024.csv"
URL = "https://apisidra.ibge.gov.br/values/t/1612/n6/all/v/109,216,214,112/p/{ano}/c81/2713"
NOMES = {"109": "area_plantada_ha", "216": "area_colhida_ha", "214": "producao_t", "112": "rendimento_kg_ha"}


def ano(a, tentativas=4):
    for k in range(tentativas):
        try:
            js = json.load(urllib.request.urlopen(URL.format(ano=a), timeout=300))
            break
        except Exception:
            if k == tentativas - 1:
                raise
            time.sleep(20 * (k + 1))
    d = pd.DataFrame(js[1:]).rename(columns={"D1C": "cod_municipio", "D1N": "municipio", "D2C": "var", "D3N": "ano", "V": "valor"})
    d["valor"] = pd.to_numeric(d.valor, errors="coerce")
    d = d.pivot_table(index=["cod_municipio", "municipio", "ano"], columns="var", values="valor").reset_index()
    d.columns.name = None
    return d.rename(columns=NOMES)


if __name__ == "__main__":
    partes = []
    for a in range(1980, 2025):
        p = ano(a)
        partes.append(p)
        print(a, len(p), "municípios com dado", flush=True)
    d = pd.concat(partes, ignore_index=True)
    d.to_csv(SAIDA, index=False)
    print(d.shape, "→", SAIDA.relative_to(RAIZ))
