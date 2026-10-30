"""Caminhos e constantes do zoneamento da soja no Paraná (dssat_zone). Reaproveita os módulos de ../calibracao_mg/src."""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

RAIZ = Path(__file__).resolve().parents[1]
CALIB = RAIZ.parent / "calibracao_mg"
sys.path.insert(0, str(CALIB / "src"))

DIR_RUN = RAIZ / ".run"
DIR_RES = RAIZ / "resultados"
DIR_DADOS = RAIZ / "dados"

# Cultivares brasileiras (fenologia v13 com 9.000 passos + crescimento da BRS 284), as mesmas dos dois artigos
CUL_BR = CALIB / "produtos" / "dados" / "SBGRO048_BR.CUL"
ECO_BR = CALIB / "produtos" / "dados" / "SBGRO048_BR.ECO"
GMR = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]


def cultivar_id(gmr):
    return f"BRG{int(round(gmr * 10)):03d}"


# Semeadura: 12 decêndios de 1/set a 31/dez (dia 1, 11 e 21 de cada mês)
DECENDIOS = [(m, d) for m in (9, 10, 11, 12) for d in (1, 11, 21)]

ANO_INI, NYERS = 1980, 45            # safras 1980/81 a 2024/25
INICIO = (3, 15)                     # SUPOSIÇÃO: início do balanço hídrico em 15/mar, como na replicação (bs2019)
