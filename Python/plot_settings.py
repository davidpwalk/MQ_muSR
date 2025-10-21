import seaborn as sns
import matplotlib.pyplot as plt

def set_demonlab_style(tex=True):
    # Seaborn base
    sns.set_theme(context="notebook", style="white")

    # Activate TeX and fourier font if requested
    if tex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": r"\usepackage{fourier}",
        })

    # Matplotlib rcParams to mimic your Plotly template
    plt.rcParams.update({
        # Figure
        "figure.figsize": (6, 4),
        "figure.dpi": 100,
        "figure.facecolor": "white",   # outer background
        "savefig.facecolor": "white",  # saved figure background

        # Axes background
        "axes.facecolor": "white",  # inner plotting area

        # Fonts
        "font.size": 16,
        "text.color": "black",

        # Axes
        "axes.edgecolor": "black",
        "axes.linewidth": 1.25,
        "axes.labelcolor": "black",
        "axes.titlesize": 16,
        "axes.titlecolor": "black",

        # # Ticks
        # "xtick.direction": "out",
        # "ytick.direction": "out",
        # "xtick.color": "black",
        # "ytick.color": "black",
        # "xtick.major.width": 1.25,
        # "ytick.major.width": 1.25,
        # "xtick.minor.visible": True,
        # "ytick.minor.visible": True,
        # "xtick.minor.width": 1.0,
        # "ytick.minor.width": 1.0,
        # "xtick.major.size": 6,
        # "ytick.major.size": 6,
        # "xtick.minor.size": 3,
        # "ytick.minor.size": 3,

        # Grid
        "axes.grid": False,

        # Legend
        "legend.frameon": True,
        "legend.facecolor": (1, 1, 1, 0.55),
        "legend.edgecolor": "black",
        "legend.framealpha": None,
        "legend.fontsize": 16,
        "legend.borderaxespad": 1,
    })

    # Use G10 qualitative colors from Plotly
    g10_colors = [
        "#3366CC", "#DC3912", "#FF9900", "#109618", "#990099",
        "#0099C6", "#DD4477", "#66AA00", "#B82E2E", "#316395"
    ]
    sns.set_palette(g10_colors)

# Define a custom plotly template
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

DemonLab_template = go.layout.Template()

DemonLab_template.layout = dict(
    template="simple_white",
    width=600,
    height=400,
    font=dict(family="URW Fourier", size=24, color="black"),
    margin=dict(l=90, r=20, t=20, b=70),
    colorway=px.colors.qualitative.G10,
    xaxis=dict(
        showline=True,
        linecolor="black",
        linewidth=1.25,
        mirror=True,
        minor=dict(ticks="outside"),
        ticks="outside",
        tickwidth=1.25,
        tickcolor="black",
        showgrid=False,
        zeroline=False,
        tickformat="~",
        ticklen=6
    ),
    yaxis=dict(
        showline=True,
        linecolor="black",
        linewidth=1.25,
        mirror=True,
        minor=dict(ticks="outside"),
        ticks="outside",
        tickwidth=1.25,
        tickcolor="black",
        showgrid=False,
        zeroline=False,
        tickformat="~",
        ticklen=6
    ),
    legend=dict(
        x=0.99, y=0.985, xanchor="right", yanchor="top",
        bgcolor="rgba(255, 255, 255, 0.55)",
        bordercolor="black",
        borderwidth=1,
        font=dict(color="black", size=20),
        tracegroupgap=10,
    ),
)

# Set this template as the default
pio.templates["DemonLab"] = DemonLab_template
