"""Escreve coeficientes numa linha de cultivar (.CUL) ou ecótipo (.ECO) do DSSAT, respeitando as colunas do cabeçalho."""
import re

FORMATOS = {"CSDL": "{:5.2f}", "PPSEN": "{:5.3f}", "EM-FL": "{:5.1f}", "FL-SD": "{:5.1f}", "SD-PM": "{:5.2f}",
            "R1PPO": "{:5.3f}"}


def substituir(texto, codigo, valores):
    """Troca, na linha que começa por `codigo`, os campos do cabeçalho @VAR#/@ECO# indicados em `valores`."""
    linhas = texto.splitlines()
    cab = None
    for i, l in enumerate(linhas):
        if l.startswith(("@VAR#", "@ECO#")):
            cab = l
        elif cab and l.startswith(codigo + " "):
            fins = {m.group(): m.end() for m in re.finditer(r"\S+", cab)}
            novo = l.ljust(len(cab))
            for nome, v in valores.items():
                fim = fins[nome]
                larg = len(FORMATOS[nome].format(v))
                novo = novo[:fim - larg] + FORMATOS[nome].format(v) + novo[fim:]
            linhas[i] = novo.rstrip() if len(l) < len(cab) else novo
            return "\n".join(linhas) + "\n"
    raise KeyError(codigo)
