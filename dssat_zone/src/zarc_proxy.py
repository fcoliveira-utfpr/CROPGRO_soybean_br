"""Índice de satisfação de necessidade de água (ISNA) no estilo do ZARC, como base de comparação do zoneamento.

Reproduz a lógica do ZARC da soja (balanço hídrico sequencial diário, ciclo fixo por grupo, ISNA = ETr/ETm na fase
crítica), não a implementação oficial (SARRA/BIPZON), cujos parâmetros exatos não são públicos:
  - CAD = AD do ZARC (mm/cm, equação oficial da IN SPA/MAPA nº 1/2022 com a textura de 0–40 cm) × ZR (50 cm,
    SUPOSIÇÃO) — aqui, e só aqui, a CAD com profundidade fixa de raízes entra, porque é assim que o ZARC trabalha;
  - ETm = Kc × ETo (ETo do Xavier, FAO-56); Kc por fração do ciclo (inicial 0,4 até 15%, sobe a 1,15 em 40%, fica até
    80%, cai a 0,5 no fim);
  - ETr = ETm × Ks, Ks = 1 se o armazenamento ≥ (1 − p)·CAD, senão linear (p = 0,5, FAO-56);
  - ciclo fixo por grupo: GMR 5,0–5,5 → 105 d; 6,0–6,5 → 120 d; 7,0–7,5 → 135 d;
  - fase crítica: 40–80% do ciclo (floração e enchimento de grãos);
  - armazenamento inicial: CAD cheia na semeadura (o ZARC inicia o balanço antes; SUPOSIÇÃO).
"""
import numpy as np
import pandas as pd

ZR = 0.5
P_FAO = 0.5
CICLO = {5.0: 105, 5.5: 105, 6.0: 120, 6.5: 120, 7.0: 135, 7.5: 135}
FASE_CRITICA = (0.40, 0.80)


# Classes de AD do ZARC (IN SPA/MAPA nº 1/2022, redação da IN nº 2/2022), em mm/cm: limite inferior de AD1…AD6
LIMITES_AD = [0.34, 0.46, 0.61, 0.80, 1.06, 1.40]


def ad_zarc(areia, silte, argila):
    """Água disponível (mm/cm) pela equação oficial do ZARC (textura em % na camada de 0–40 cm)."""
    at, sil, arg = (np.asarray(v, float) for v in (areia, silte, argila))
    x = (-0.02128887 * at - 0.01005814 * sil - 0.01901894 * arg + 0.0001171219 * at * sil
         + 0.0002073924 * at * arg + 0.00006118707 * sil * arg - 0.000006373789 * at * sil * arg)
    return (1 + 0.3591 * x) ** 2.78474 * 10


def classe_ad(ad):
    """'AD1'…'AD6' (abaixo de 0,34 mm/cm: 'AD0', fora das classes)."""
    return np.array([f"AD{int(np.searchsorted(LIMITES_AD, v, side='right'))}" for v in np.atleast_1d(ad)])


def ad_perfil(p):
    """AD do ZARC de um perfil do MapBiomas (médias das camadas de 0–40 cm)."""
    cam = ["000_010", "010_020", "020_030", "030_040"]
    return float(ad_zarc(np.mean([p[f"sand_{c}"] for c in cam]), np.mean([p[f"silt_{c}"] for c in cam]),
                         np.mean([p[f"clay_{c}"] for c in cam])))


def kc(frac):
    return np.interp(frac, [0, 0.15, 0.40, 0.80, 1.0], [0.4, 0.4, 1.15, 1.15, 0.5])


def isna(chuva, eto, cad, ciclo):
    """ISNA da fase crítica para uma semeadura (arrays diários a partir da semeadura, comprimento ≥ ciclo)."""
    arm = cad
    etr_c = etm_c = 0.0
    for d in range(ciclo):
        f = d / ciclo
        etm = kc(f) * eto[d]
        ks = 1.0 if arm >= (1 - P_FAO) * cad else arm / ((1 - P_FAO) * cad)
        etr = min(etm * ks, arm + chuva[d])
        arm = min(cad, max(0.0, arm + chuva[d] - etr))
        if FASE_CRITICA[0] <= f < FASE_CRITICA[1]:
            etr_c += etr
            etm_c += etm
    return etr_c / etm_c if etm_c > 0 else np.nan


def isna_vetor(chuva, eto, cad, ciclo):
    """Mesmo balanço de isna(), vetorizado: chuva e eto com forma (n_semeaduras, ciclo)."""
    arm = np.full(chuva.shape[0], float(cad))
    etr_c = np.zeros(chuva.shape[0])
    etm_c = np.zeros(chuva.shape[0])
    for d in range(ciclo):
        f = d / ciclo
        etm = kc(f) * eto[:, d]
        ks = np.where(arm >= (1 - P_FAO) * cad, 1.0, arm / ((1 - P_FAO) * cad))
        etr = np.minimum(etm * ks, arm + chuva[:, d])
        arm = np.clip(arm + chuva[:, d] - etr, 0.0, cad)
        if FASE_CRITICA[0] <= f < FASE_CRITICA[1]:
            etr_c += etr
            etm_c += etm
    return np.where(etm_c > 0, etr_c / np.where(etm_c > 0, etm_c, 1), np.nan)


def isna_celula(clima, ad, gmrs, semeaduras, anos):
    """DataFrame (gmr, semeadura "mm-dd", ano, indice_isna) para uma célula e um solo (AD do ZARC em mm/cm)."""
    cad = ad * ZR * 100
    chuva, eto = clima.pr.clip(lower=0).to_numpy(), clima.ET.clip(lower=0).to_numpy()
    idx0 = clima.index[0]
    linhas = []
    for g in gmrs:
        n = CICLO[g]
        chaves, ini = [], []
        for a in anos:
            for m, d in semeaduras:
                i = (pd.Timestamp(a, m, d) - idx0).days
                if 0 <= i and i + n <= len(chuva):
                    chaves.append((g, f"{m:02d}-{d:02d}", a))
                    ini.append(i)
        j = np.asarray(ini)[:, None] + np.arange(n)[None]
        v = isna_vetor(chuva[j], eto[j], cad, n)
        linhas += [(*k, x) for k, x in zip(chaves, v)]
    return pd.DataFrame(linhas, columns=["gmr", "semeadura", "ano", "indice_isna"])
