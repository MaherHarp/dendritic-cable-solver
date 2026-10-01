"""Render figures from results/results.json (run `npm run exp` first)."""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "results", "results.json")))
OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID = "#2b2b2a", "#6f6e69", "#e4e3dd"
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2, "lines.markersize": 6, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)


# 1. Solver scaling.
e1 = R["e1"]
fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.loglog(e1["hinesSizes"], e1["hinesMs"], "o-", color=BLUE, label="Hines tree solve, O(n)")
ax.loglog(e1["denseSizes"], e1["denseMs"], "s-", color=ORANGE, label="dense LU, O(n³)")
ax.set_xlabel("nodes n")
ax.set_ylabel("wall-clock per solve (ms)")
ax.set_title("Tree solve vs dense LU (Node, single thread)", loc="left")
ax.legend()
save(fig, "solver_scaling.png")

# 2. Spatial convergence.
conv = R["e2"]["conv"]
fig, ax = plt.subplots(figsize=(6.2, 4.2))
colors = {"cable": BLUE, "rall": ORANGE, "random": AQUA}
for key, d in conv.items():
    geom, mode = key.split("/")
    ls, mk = ("-", "o") if mode == "junction" else ("--", "^")
    ax.loglog(d["nsegs"], d["errs"], ls, marker=mk, color=colors[geom], label=f"{geom}, {mode} branch points")
ns = conv["random/junction"]["nsegs"]
ax.loglog(ns, [0.3 / n**2 for n in ns], ":", color=MUTED, lw=1)
ax.loglog(ns, [0.05 / n for n in ns], ":", color=MUTED, lw=1)
ax.annotate("slope −2", (ns[-2], 0.3 / ns[-2] ** 2), color=MUTED, xytext=(4, -10), textcoords="offset points")
ax.annotate("slope −1", (ns[-2], 0.05 / ns[-2]), color=MUTED, xytext=(4, 4), textcoords="offset points")
ax.set_xlabel("compartments per section")
ax.set_ylabel("relative error in input resistance")
ax.set_title("Junction nodes restore second-order accuracy at branch points", loc="left")
ax.legend(fontsize=8)
save(fig, "spatial_convergence.png")

# 3. Temporal convergence.
res = R["e3"]["res"]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=False)
for ax, label in zip(axes, ["passive", "HH"]):
    d = res[label]
    ax.loglog(d["dts"], d["be"], "s-", color=ORANGE, label="backward Euler")
    ax.loglog(d["dts"], d["cn"], "o-", color=BLUE, label="staggered Crank–Nicolson")
    ax.set_xlabel("dt (ms)")
    ax.set_title(f"{label} membrane", loc="left")
axes[0].set_ylabel("max |V − V_ref| (mV)")
axes[0].legend()
save(fig, "temporal_convergence.png")

# 4. CN ringing.
ring = R["e3"]["ring"]
plateau = R["e3"]["plateau"]
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ref = ring["ref"]
ax.plot([i * ref["dt"] for i in range(len(ref["v"]))], [v / plateau for v in ref["v"]], color=MUTED, lw=1.5, label="reference (λ·dt < 1)")
for key, color, mk in [("be@0.5", ORANGE, "s"), ("cn@0.5", BLUE, "o")]:
    d = ring[key]
    t = [i * d["dt"] for i in range(len(d["v"]))]
    ax.plot(t, [v / plateau for v in d["v"]], marker=mk, color=color, label=key.replace("@", ", dt = ").upper().replace("DT", "dt") + " ms")
ax.set_xlim(0, 6)
ax.set_xlabel("time (ms)")
ax.set_ylabel("root voltage / plateau")
ax.set_title("CN rings on stiff modes at large dt; BE is monotone", loc="left")
ax.legend()
save(fig, "cn_ringing.png")

# 5. Rall equivalence transient.
rall = R["e2"]["rall"]
tr = rall["traces"]
dt = rall["sampleDt"]
fig, ax = plt.subplots(figsize=(6.2, 3.8))
for key, color, ls, lab in [("cylinder", MUTED, "-", "equivalent cylinder"), ("tree", BLUE, "--", "3/2-rule tree"),
                            ("broken", ORANGE, "-", "daughters ×1.3 diameter")]:
    ax.plot([i * dt for i in range(len(tr[key]))], tr[key], ls, color=color, label=lab)
ax.set_xlabel("time (ms)")
ax.set_ylabel("root voltage (mV above rest)")
ax.set_title("Rall's equivalent cylinder holds through the transient", loc="left")
ax.legend()
save(fig, "rall_transient.png")

# 6. Conduction velocity.
e5 = R["e5"]
fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.loglog(e5["diams"], e5["vs"], "o", color=BLUE, label=f"simulated (dx = {e5['dx']} µm)")
k = e5["vsScaled"][0] / e5["diams"][0] ** 0.5
ax.loglog(e5["diams"], [k * d**0.5 for d in e5["diams"]], "-", color=MUTED, lw=1, label="∝ √d")
ax.loglog([476], [18.8], "D", color=ORANGE, markersize=8, label="Hodgkin & Huxley 1952 (computed)")
ax.set_xlabel("axon diameter (µm)")
ax.set_ylabel("conduction velocity (m/s)")
ax.set_title("Unmyelinated HH axon, 18.5 °C", loc="left")
ax.legend()
save(fig, "conduction_velocity.png")

print("figures written to", OUT)
