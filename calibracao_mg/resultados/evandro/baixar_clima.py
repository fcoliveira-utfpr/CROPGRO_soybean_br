"""Clima do Xavier (pixel 0,1°) 1980–2025 nas células dos experimentos de Evandro Silva."""
import sys
sys.path.insert(0, "src")
import pandas as pd
import gee
gee.iniciar()
cel = pd.DataFrame({"lon": [-47.65, -42.75, -55.75, -53.55, -53.75], "lat": [-22.75, -5.05, -12.75, -28.65, -28.85]})
print("falhas:", gee.baixar_series(cel, 1980, 2025, no=False, n_threads=5), flush=True)
