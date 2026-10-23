# 06_model_diagram.py - JAGS model diagram -of Kruschke (2015) / Baath's distribution diagrams

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SKY, BLUE, ORANGE, INK, GREY = "#8fcbea", "#2a78d6", "#eb6834", "#14232b", "#5b6b73"
H = 8.0
fig, ax = plt.subplots(figsize=(8.0, 9.3))
ax.set_xlim(0, 100); ax.set_ylim(-14.5, 110); ax.axis("off")

def glyph(x, yb, kind):
    if kind == "normal":            # symmetric bell
        xs = np.linspace(-3.2, 3.2, 200); ys = np.exp(-xs ** 2 / 2); xp = x + xs / 3.2 * 5.5
    else:                           # gamma: skewed, steep rise then long right tail
        xs = np.linspace(0, 6, 300); ys = xs ** 0.5 * np.exp(-xs); ys = ys / ys.max(); xp = x - 5.5 + xs / 6 * 11
    ax.plot(xp, yb + H * ys, color=SKY, lw=2.6, solid_capstyle="round")
    ax.plot([x - 5.8, x + 5.8], [yb, yb], color="black", lw=2.0)

def node(x, yb, kind, mean=None, prec=None, shaperate=None):
    glyph(x, yb, kind)
    if kind == "normal":
        ax.text(x, yb + 2.9, mean, ha="center", va="center", fontsize=11.5, color=INK)
        ax.text(x + 6.3, yb + 6.0, prec, ha="center", va="center", fontsize=11.5, color=INK)
        ax.text(x, yb + 0.9, "normal", ha="center", va="center", fontsize=9, color=INK)
    else:
        ax.text(x + 2.4, yb + 4.8, shaperate, ha="center", va="center", fontsize=11, color=INK)
        ax.text(x + 0.6, yb + 0.9, "gamma", ha="center", va="center", fontsize=9, color=INK)

def arrow(p0, p1, sym, dx=1.4, dy=0.0):
    col = BLUE if sym == "~" else ORANGE
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="-|>", lw=1.5, color=INK, shrinkA=0, shrinkB=0, mutation_scale=15))
    mx, my = (p0[0] + p1[0]) / 2 + dx, (p0[1] + p1[1]) / 2 + dy
    ax.text(mx, my, sym, ha="left", va="center", fontsize=20 if sym == "~" else 17, fontweight="bold", color=col,
            bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none"))

def sym(x, y, txt, size=14):
    ax.text(x, y, txt, ha="center", va="center", fontsize=size, color=INK)

def tag(x, y, txt):
    ax.text(x, y, txt, ha="left", va="center", fontsize=9.5, style="italic", color=GREY)

#  legend 
ax.text(3, 108, r"$\sim$", fontsize=18, fontweight="bold", color=BLUE, va="center")
ax.text(7, 108, "is distributed as (stochastic)", fontsize=11, va="center", color=INK)
ax.text(54, 108, "=", fontsize=17, fontweight="bold", color=ORANGE, va="center")
ax.text(57.5, 108, "is defined as (deterministic)", fontsize=11, va="center", color=INK)

#  row 1: priors 
R1 = 93; R2 = 62; TOP2 = R2 + H + 0.4
prior_normal = [11.5, 46.0, 78.0, 92.0]
for x in prior_normal:
    node(x, R1, "normal", "0", "1.0E-4")
node(26.0, R1, "gamma", shaperate="0.01, 0.01")
node(60.0, R1, "gamma", shaperate="0.01, 0.01")
tag(1.5, 82, "k = 1 ... 6")
arrow((11.5, R1 - 1.6), (17.0, TOP2), "~", dx=0.8)       # prior -> alpha_k mean
arrow((26.0, R1 - 1.6), (22.8, TOP2), "~", dx=0.8)       # gamma -> tau_a
arrow((46.0, R1 - 1.6), (51.0, TOP2), "~", dx=0.8)       # prior -> beta_k mean
arrow((60.0, R1 - 1.6), (56.8, TOP2), "~", dx=0.8)       # gamma -> tau_b
arrow((78.0, R1 - 1.6), (78.0, 74.8), "~"); arrow((92.0, R1 - 1.6), (92.0, 74.8), "~")
sym(78, 71.8, r"$q_k$"); sym(92, 71.8, r"$g_k$")

#  row 2: country-level distributions (partial pooling) 
node(18, R2, "normal", r"$\alpha_k$", r"$\tau_a$")
node(52, R2, "normal", r"$\beta_k$", r"$\tau_b$")
for x, s in [(18, r"$a_i$"), (52, r"$b_i$")]:
    arrow((x, R2 - 1.6), (x, 49.2), "~")
    sym(x, 46.2, s, 15)
tag(27, 61, "i = 1 ... 236")

# row 3: deterministic mean 
ax.add_patch(plt.Rectangle((5, 22), 92, 9.5, fc="#edf2f3", ec=INK, lw=1.4))
ax.text(51, 26.8, r"$\mu_n \;=\; a_i \;+\; b_i\, t_n \;+\; q_k\, t_n^{2} \;+\; g_k\, \mathrm{COVID}_n$", ha="center", va="center", fontsize=14, color=INK)
arrow((18, 43.4), (18, 31.9), "=")
arrow((52, 43.4), (52, 31.9), "=")
arrow((78, 69.2), (78, 31.9), "=", dy=-10)
arrow((92, 69.2), (92, 31.9), "=", dy=-10)

# row 4: likelihood 
R4 = 3
node(46, R4, "normal", r"$\mu_n$", r"$\tau$")
arrow((46, 21.6), (46, R4 + H + 0.6), "=", dx=1.4, dy=-0.5)
node(84, R4 + 2, "gamma", shaperate="0.01, 0.01")
arrow((78.2, R4 + 6.5), (52.6, R4 + 5.6), "~", dx=-1.0, dy=2.4)
arrow((46, R4 - 1.6), (46, -7.2), "~")
sym(46, -9.4, r"$y_n$", 17)
tag(54, -2.5, "n = 1 ... 17,464 country-years")

ax.text(1, -13.6, "Parameters beside each distribution: (mean, precision) for normal, (shape, rate) for gamma. Curves are schematic.",
        fontsize=8.8, color=GREY, va="bottom")
plt.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.005)
plt.savefig("fig/fig_jags_diagram.png", dpi=300, facecolor="white")
print("saved")
