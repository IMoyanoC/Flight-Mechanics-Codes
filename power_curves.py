import pandas as pd
import matplotlib.pyplot as plt
import re

# ============================================================
# READ EXCEL FILE
# ============================================================

df = pd.read_excel("potencias.xlsx")

# ============================================================
# IDENTIFY VELOCITY AND ALTITUDES
# ============================================================

V = df["V [m/s]"]

# Detect altitudes automatically from columns like:
# Preq_0m, Pdisp_0m, Preq_500m, Pdisp_500m, etc.
altitudes = []

for col in df.columns:
    match = re.match(r"Preq_(\d+)m", str(col))

    if match:
        altitudes.append(int(match.group(1)))

altitudes = sorted(altitudes)

# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(10, 7))

for h in altitudes:

    Preq = df[f"Preq_{h}m"]
    Pdisp = df[f"Pdisp_{h}m"]

    # Required power: solid line
    linea, = plt.plot(
        V,
        Preq,
        linewidth=2,
        label=f"{h} m"
    )

    # Available power: dashed line, same color
    plt.plot(
        V,
        Pdisp,
        linestyle="--",
        linewidth=2,
        color=linea.get_color()
    )

# ============================================================
# FORMAT
# ============================================================

plt.xlabel("Velocidad [m/s]")
plt.ylabel("Potencia [W]")
plt.title("Potencia requerida y disponible vs Velocidad")

plt.grid(True)
plt.legend(title="Altitud")

plt.tight_layout()
plt.show()
