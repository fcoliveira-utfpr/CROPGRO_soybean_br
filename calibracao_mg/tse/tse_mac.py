"""Roda o TSE v4.3 original (Memic et al., 2021; github.com/DSSAT/TSE) no macOS, sem alterar o algoritmo.

O TSE é um programa Python 3.7 + PyQt5 empacotado para Windows (TSE.exe, PyInstaller). Aqui o bytecode original
(TSE_v4.3.pyc, extraído do TSE.exe) é executado pelo Python 3.7, com duas camadas finas por fora do código do TSE:

1. Tradução de caminhos. O TSE supõe o DSSAT em C:\\DSSAT48 e monta todos os caminhos com "\\". As funções do
   sistema que ele usa (open, os.*, shutil.*, glob.glob, subprocess.*) são envolvidas para traduzir C:\\DSSAT48\\...
   em RAIZ/DSSAT48/... (e de volta, nos retornos de os.getcwd e glob.glob). A chamada ao DSCSM048.EXE vira o
   executável do DSSAT 4.8 para macOS (o mesmo do DSSATTools e do dssat_zone), com o arquivo batch traduzido.
2. Interface sem tela. Qt em modo offscreen; as caixas de diálogo (QMessageBox) são registradas no log em vez de
   esperar um clique. O roteiro (roteiro_*) dispara os mesmos sinais que os cliques do usuário, na ordem do manual.

Uso: tools/py37/bin/python3.7 tse_mac.py <roteiro.json>
"""
import builtins
import glob
import importlib.machinery
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

AQUI = Path(__file__).resolve().parent
RAIZ = Path(os.environ.get("TSE_RAIZ", AQUI / "raiz"))  # contém DSSAT48/ (instalação simulada; fora do git)
PYC = AQUI / "pyc"                                    # TSE_v4.3.pyc, calibrator_Interface.pyc, CoefficientCombinations.pyc
BIN = Path(os.environ.get("DSSAT_BIN", "/Users/coted-sh/Python/DSSAT/calibracao_mg/.venv/lib/python3.11/"
                                       "site-packages/DSSATTools/bin/dscsm048"))
CURTO = os.environ.get("TSE_CURTO", "/tmp/tse48")   # o DSSAT lê o DSSATPRO em colunas fixas: caminho longo
DSSAT_HOME = CURTO + os.sep                          # empurra o nome do modelo para fora do campo
LOG = RAIZ / "log_tse.txt"

_open, _exists, _isfile, _isdir = builtins.open, os.path.exists, os.path.isfile, os.path.isdir
_mkdir, _makedirs, _remove, _chdir, _listdir, _getcwd, _rename = (os.mkdir, os.makedirs, os.remove, os.chdir,
                                                                   os.listdir, os.getcwd, os.rename)
_copy, _copyfile, _copy2, _move, _rmtree = shutil.copy, shutil.copyfile, shutil.copy2, shutil.move, shutil.rmtree
_glob, _glob1 = glob.glob, glob.glob1
_call, _popen, _run = subprocess.call, subprocess.Popen, subprocess.run


def log(*a):
    with _open(LOG, "a") as f:
        f.write(" ".join(str(x) for x in a) + "\n")


def mac(p):
    """C:\\DSSAT48\\x, /DSSAT48/x ou C:\\DSSAT48/x (os.path.join misturado) → RAIZ/DSSAT48/x."""
    if not isinstance(p, str):
        return p
    q = p.replace("\\", "/")
    if q[:2].upper() == "C:":
        q = q[2:]
    q = q.replace("//", "/")
    if q.startswith("/DSSAT48"):
        return str(RAIZ) + q
    return q if "\\" in p else p


def win(p):
    s = str(RAIZ)
    return "C:" + p[len(s):].replace("/", "\\") if isinstance(p, str) and p.startswith(s) else p


def _abrir(f, *a, **k):
    return _open(mac(f), *a, **k)


def _dssat(cmd):
    """DSCSM048.EXE <modelo> B <batch> → executável do macOS, com o batch reescrito em caminhos do Mac."""
    partes = cmd.split()
    modelo, modo, batch = partes[1], partes[2], partes[3]
    cwd = _getcwd()
    texto = _open(batch, encoding="utf-8", errors="ignore").read()
    linhas = []
    for linha in texto.splitlines():                  # caminho do FileX nas colunas 1–94; o DSSAT corta caminhos
        if "\\" in linha or linha.strip().upper().startswith("C:"):        # longos, e ele roda na pasta do FileX
            caminho, resto = linha[:94].strip(), linha[94:]
            alvo = Path(mac(caminho))
            linha = "{:<94}".format(alvo.name if alvo.parent == Path(cwd) else str(alvo)) + resto
        linhas.append(linha)
    with _open(batch, "w") as f:
        f.write("\n".join(linhas) + "\n")
    _copy(str(RAIZ / "DSSAT48" / "DSSATPRO.L48"), cwd)
    t0 = time.time()
    p = _popen([str(BIN), modelo, modo, batch], cwd=cwd, env={"DSSAT_HOME": DSSAT_HOME},     # Popen original: o
               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)                              # run() usaria o substituto
    p.communicate()
    log(f"DSSAT {modelo} {modo} {batch} em {cwd}: {time.time() - t0:.2f} s, código {p.returncode}")
    return p.returncode


def _chamar(cmd, *a, **k):
    if isinstance(cmd, str) and "DSCSM048" in cmd.upper():
        return _dssat(cmd)
    log("subprocess ignorado:", cmd)                  # GBuild e afins (Windows)
    return 0


class _Popen:
    def __init__(self, cmd, *a, **k):
        self.returncode = _chamar(cmd)

    def wait(self, *a, **k):
        return self.returncode


def instalar_traducao():
    builtins.open = _abrir
    os.path.exists = lambda p: _exists(mac(p))
    os.path.isfile = lambda p: _isfile(mac(p))
    os.path.isdir = lambda p: _isdir(mac(p))
    os.mkdir = lambda p, *a, **k: _mkdir(mac(p), *a, **k)
    os.makedirs = lambda p, *a, **k: _makedirs(mac(p), *a, **k)
    os.remove = lambda p: _remove(mac(p))
    os.chdir = lambda p: _chdir(mac(p))
    os.listdir = lambda p=".": _listdir(mac(p))
    os.getcwd = lambda: win(_getcwd())
    os.rename = lambda a, b: _rename(mac(a), mac(b))
    shutil.copy = lambda a, b, **k: _copy(mac(a), mac(b), **k)
    shutil.copyfile = lambda a, b, **k: _copyfile(mac(a), mac(b), **k)
    shutil.copy2 = lambda a, b, **k: _copy2(mac(a), mac(b), **k)
    shutil.move = lambda a, b, **k: _move(mac(a), mac(b), **k)
    shutil.rmtree = lambda p, *a, **k: _rmtree(mac(p), *a, **k)
    glob.glob = lambda pat, **k: [win(x) for x in _glob(mac(pat), **k)]
    glob.glob1 = lambda d, pat: _glob1(mac(d), pat)
    subprocess.call = _chamar
    subprocess.Popen = _Popen
    subprocess.run = lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, _chamar(cmd))


def silenciar_dialogos(QtW):
    def registrar(nome):
        def f(*a, **k):
            log(f"QMessageBox.{nome}:", " | ".join(str(x) for x in a[1:3]))
            return QtW.QMessageBox.Ok
        return f
    for nome in ("warning", "information", "critical", "question", "about"):
        setattr(QtW.QMessageBox, nome, staticmethod(registrar(nome)))
    QtW.QMessageBox.exec = lambda self: log("QMessageBox.exec:", self.text()) or QtW.QMessageBox.Ok
    QtW.QMessageBox.exec_ = QtW.QMessageBox.exec
    QtW.QDialog.exec = lambda self: log("QDialog.exec") or 1
    QtW.QDialog.exec_ = QtW.QDialog.exec


def carregar_tse():
    """Importa o bytecode original como módulo (não como __main__: o roteiro cria a janela)."""
    sys.path.insert(0, str(PYC))
    os.getcwd = lambda: "C:\\DSSAT48\\Tools\\TSE"     # o TSE deduz a instalação do diretório em que roda
    loader = importlib.machinery.SourcelessFileLoader("tse_v43", str(PYC / "TSE_v4.3.pyc"))
    spec = importlib.util.spec_from_loader("tse_v43", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    os.getcwd = lambda: win(_getcwd())
    mod.pathToCropForLater = "C:\\DSSAT48\\Tools\\TSE"
    return mod


if __name__ == "__main__":
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    _open(LOG, "w").close()
    if os.path.islink(CURTO):
        _remove(CURTO)
    os.symlink(str(RAIZ / "DSSAT48"), CURTO)
    instalar_traducao()
    from PyQt5 import QtWidgets as QtW
    silenciar_dialogos(QtW)
    app = QtW.QApplication(sys.argv)
    tse = carregar_tse()
    janela = tse.Window()
    import roteiro
    roteiro.executar(janela, tse, json.load(_open(sys.argv[1])), log)
