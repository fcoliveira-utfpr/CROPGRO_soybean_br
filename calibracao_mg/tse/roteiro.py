"""Roteiro do TSE: dispara os sinais da interface na ordem que o usuário seguiria (manual do TSE v4.1, seção III).

Configuração (JSON):
  cultura, modelo        nome na lista de culturas (p. ex. "Soybean") e módulo ("CRGRO")
  cultivar               código ou nome na lista de cultivares (p. ex. "BRAGG")
  tratamentos            linhas da lista de tratamentos, como "UFGA7801.SBX" e número do TRT
  modo                   "crescimento" (variáveis do FileT) ou "fenologia" (eventos do FileA)
  variaveis              códigos do FileT (modo crescimento) ou eventos (ADAP, PD1P, PDFP, MDAP)
  coeficientes           {nome: [min, max, n_comb]}
  metodo                 "faixa" (Run Range setup) ou "reducao" (Run Range Reduction)
  ate                    para depurar: "init", "cultivar", "setup" ou "coef" (para antes de rodar)
"""
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication


def processar():
    for _ in range(3):
        QApplication.processEvents()


def clicar_item(lista, texto, log, contem=True):
    for i in range(lista.count()):
        it = lista.item(i)
        if (texto in it.text()) if contem else (it.text().strip() == texto):
            lista.setCurrentRow(i)
            it.setSelected(True)
            lista.clicked.emit(lista.indexFromItem(it))
            processar()
            return it
    raise RuntimeError(f"'{texto}' não está na lista: {[lista.item(i).text() for i in range(lista.count())][:30]}")


def estado(ui, log, titulo):
    log(f"── {titulo}")
    log("textBrowser_2:\n" + ui.textBrowser_2.toPlainText()[-3000:])


def executar(w, tse, c, log):
    ui = w.ui
    ui.lineEdit_58.setText("C:\\DSSAT48")
    # 1. cultura e modelo
    clicar_item(ui.listWidget_2, c["cultura"], log, contem=False)
    i = ui.comboBox_6.findText(c["modelo"])
    ui.comboBox_6.setCurrentIndex(i)
    ui.comboBox_6.activated.emit(i)
    processar()
    log("modelo:", w.modelNew, w.cultivarFileNew, w.cropModelDCSMNew, w.cropDirNew)
    if c.get("arquivo") == "ECO":                    # caixa "Ecotype" da seção 1 (coeficientes do .ECO)
        ui.checkBox_12.setChecked(False)
        ui.checkBox_13.setChecked(True)
        ui.checkBox_13.clicked.emit(True)
        processar()
    log("cultivares:", [ui.listWidget_7.item(k).text() for k in range(ui.listWidget_7.count())][:15])
    if c.get("ate") == "init":
        return estado(ui, log, "init")
    # 2. cultivar, FileX e tratamentos
    clicar_item(ui.listWidget_7, c["cultivar"], log)
    log("FileX:", [ui.listWidget_8.item(k).text() for k in range(ui.listWidget_8.count())])
    ui.listWidget_9.clear()
    arquivos = {x for x, _ in c["tratamentos"]}       # seleção múltipla dos FileX (Ctrl+clique no manual)
    ui.listWidget_8.clearSelection()
    ultimo = None
    for k in range(ui.listWidget_8.count()):
        it = ui.listWidget_8.item(k)
        if it.text()[:12] in arquivos:
            it.setSelected(True)
            ultimo = it
    ui.listWidget_8.clicked.emit(ui.listWidget_8.indexFromItem(ultimo))
    processar()
    log("TRT:", ui.listWidget_9.count(), "tratamentos listados")
    for k in range(ui.listWidget_9.count()):
        it = ui.listWidget_9.item(k)
        if any(it.text().startswith(x) and int(it.text()[14:16]) == trt for x, trt in c["tratamentos"]):
            it.setSelected(True)
            ui.listWidget_9.clicked.emit(ui.listWidget_9.indexFromItem(it))
    processar()
    if c.get("ate") == "cultivar":
        return estado(ui, log, "cultivar")
    ui.pushButton_11.click()                         # 2.4 Check the DSSAT file setup
    processar()
    if c.get("ate") == "setup":
        return estado(ui, log, "setup")
    # 3. variáveis-alvo (crescimento: FileT; fenologia: eventos do FileA)
    if c["modo"] == "crescimento":
        ui.pushButton_15.click()
        processar()
        ui.listWidget_5.clearSelection()
        for v in c["variaveis"]:
            for k in range(ui.listWidget_5.count()):
                it = ui.listWidget_5.item(k)
                if v in it.text().split():            # "Grain wt kg/ha  GWAD ..."
                    it.setSelected(True)
        processar()
        ui.pushButton_9.click()                      # Initialize variable/s
    else:
        ui.pushButton_14.click()
        processar()
        caixas = {"ADAP": ui.checkBox_8, "PD1P": ui.checkBox_4, "PDFP": ui.checkBox_10, "MDAP": ui.checkBox_9}  # 1º–4º evento
        if c["variaveis"] == "todos":
            # "All selected": marca os eventos com rótulo (os observados em todos os tratamentos). As caixas de
            # cada evento são exclusivas — clicar numa desmarca as outras —, então só esta caixa usa vários eventos.
            ui.checkBox_7.setChecked(True)
            ui.checkBox_7.clicked.emit(True)
        else:
            for v in c["variaveis"]:
                caixas[v].setChecked(True)
                caixas[v].clicked.emit(True)
        processar()
        log("eventos marcados:", {k: cx.isChecked() for k, cx in caixas.items()}, "todos:", ui.checkBox_7.isChecked())
    processar()
    log("alvos:", [getattr(ui, f"lineEdit_{k}").text() for k in (33, 34, 35, 36)],
        "fenologia:", [getattr(ui, f"lineEdit_{k}").text() for k in (54, 57, 59, 60)])
    # 4. coeficientes: seleção, mínimo, máximo e número de combinações
    ui.listWidget.clearSelection()
    for nome in c["coeficientes"]:
        for k in range(ui.listWidget.count()):
            it = ui.listWidget.item(k)
            if it.text().split()[0] == nome:
                it.setSelected(True)
    processar()
    ui.pushButton_5.click()                          # Initialize
    processar()
    rotulos = (46, 38, 45, 40, 39, 30, 31)
    campos = ((3, 4, 24), (6, 7, 25), (9, 10, 26), (12, 13, 27), (15, 16, 28), (18, 19, 29), (21, 22, 32))
    for rot, (mn, mx, nc) in zip(rotulos, campos):
        nome = getattr(ui, f"lineEdit_{rot}").text().strip()
        if nome in c["coeficientes"]:
            a, b, n = c["coeficientes"][nome]
            getattr(ui, f"lineEdit_{mn}").setText(str(a))
            getattr(ui, f"lineEdit_{mx}").setText(str(b))
            getattr(ui, f"lineEdit_{nc}").setText(str(n))
            processar()
            log(f"coeficiente {nome}: min {a} max {b} n {n} → passo",
                getattr(ui, f"lineEdit_{campos[rotulos.index(rot)][0] + 2}").text())
    if c.get("ate") == "coef":
        return estado(ui, log, "coef")
    ui.pushButton_4.click()                          # 4.4 Check input setup
    processar()
    estado(ui, log, "input setup")
    # 5. execução
    if c.get("metodo", "faixa") == "reducao":
        ui.pushButton_24.click()                     # Run Range Reduction
    else:
        ui.pushButton_2.click()                      # Run Range setup
    processar()
    estado(ui, log, "resultado")
    if c.get("multi_trt", True) and ui.pushButton_7.isEnabled():
        ui.pushButton_7.click()                      # Multi TRT comb (calculate_mean; manual, apêndice A.2)
        processar()
        estado(ui, log, "multi-TRT")
    else:
        log("Multi TRT comb indisponível (pushButton_7 desabilitado)")
