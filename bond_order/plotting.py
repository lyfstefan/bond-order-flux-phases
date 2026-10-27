"""Plots return figures; no global rcParams, TeX requirement or show() calls."""
import numpy as np
import matplotlib.pyplot as plt


def plot_bands(model,currents):
    bands=model.bands(currents)
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,band,label in zip(axes,bands,["Lower band","Upper band"]):
        im=ax.pcolormesh(model.kx/np.pi,model.ky/np.pi,band,shading="auto",cmap="viridis")
        ax.set(xlabel="kx / pi",ylabel="ky / pi",title=label)
        fig.colorbar(im,ax=ax,label="Energy (same units as t and V)")
    fig.suptitle(f"Four-current model: {model.mesh} x {model.mesh} mesh, t={model.t:g}, t'={model.tp:g}, V={model.interaction:g}")
    return fig
