import numpy as np
import pylab as pl


def pltimg(ax, xarr, yarr, zarr, xlabel, ylabel, zlabel, cmap,
           CB_levels, CB_ticklabels, flag_contour):
    min_val, max_val = np.amin(zarr), np.amax(zarr)
    ax.set_xlabel(xlabel, labelpad=-2)
    ax.set_ylabel(ylabel)
    im = ax.imshow(zarr.transpose(),
                   interpolation='bicubic', origin='lower',
                   cmap=cmap, aspect='auto', alpha=0.7,
                   extent=(min(xarr), max(xarr),
                           min(yarr), max(yarr)))
    im.set_clim(vmin=min_val, vmax=max_val)

    CB = pl.colorbar(im, ax=ax, ticks=CB_levels)
    CB.ax.set_yticklabels(CB_ticklabels)
    CB.ax.set_ylabel(zlabel, labelpad=3)
    CB.ax.minorticks_off()

    if flag_contour:
        X, Y = np.meshgrid(xarr, yarr)
        CS = ax.contour(X, Y, zarr.transpose(),
                        CB_levels, linestyles='solid',
                        colors='k', linewidths=2, alpha=0.5)
        fmt = {}
        for l, s in zip(CS.levels, CB_ticklabels):
            fmt[l] = s
        pl.clabel(CS, CS.levels, inline=True, fmt=fmt,
                  fontsize=30, colors=None)