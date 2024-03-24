# import woma as woma1
import woma
from gadget import Snapshot
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import numpy as np
import h5py
from IPython.display import Video
from IPython.display import HTML
import swiftsimio as sw
from unyt import m, cm, s, g, kg
import unyt
import sys
import seaborn as sns

# from plot_solution import load_snapshot, plot_snapshot
import scipy
from copy import copy, deepcopy
import subprocess
import os

# import scipy.interpolate
from scipy import interpolate
from IPython.display import Video
import time
import pandas as pd
import subprocess
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

plt.rcParams["patch.force_edgecolor"] = True
woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "SS08_water"])
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import colorcet as cc
from planetboundmass import Bound, Snap
import cmcrameri.cm as cmc
from matplotlib.pyplot import get_cmap
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from slurmpy import Slurm
import subprocess
from pathlib import Path
from eostable import extEOStable
import colormaps as local_cmaps


def plot_spherical_profiles(
    planet=None,
    useU=False,
    Till=False,
    planet1loc=None,
    planet2loc=None,
    loadfromfile=False,
):

    fig, ax = plt.subplots(1, 5, figsize=(20, 4.5))

    if not loadfromfile:

        # planet.calculate_entropies(useU=useU)

        colour = np.empty(len(planet.A1_r), dtype=object)
        for id_c, c in Di_id_colour.items():
            colour[planet.A1_mat_id == id_c] = c

        ax[0].scatter(planet.A1_r / R_earth, planet.A1_rho, s=1, c=colour)
        ax[0].set_xlabel(r"Radius $[R_\oplus]$", fontsize=14)
        ax[0].set_ylabel(r"Density [kg m$^{-3}$]", fontsize=14)
        # ax[0, 0].set_yscale("log")
        # ax[0, 0].set_xlim(0, None)

        ax[1].scatter(
            planet.A1_r / R_earth, np.log10(1e-9 * planet.A1_P), c=colour, s=0.1
        )
        ax[1].set_xlabel(r"Radius $[R_\oplus]$", fontsize=14)
        ax[1].set_ylabel(r"Pressure [log$_{10}$(GPa)]", fontsize=14)
        # ax[1].set_yscale("log")
        # ax[0, 1].set_xlim(0, None)

        ax[2].scatter(planet.A1_r / R_earth, planet.A1_u, s=1, c=colour)
        ax[2].set_xlabel(r"Radius $[R_\oplus]$", fontsize=14)
        # ax[1, 0].set_ylabel(r"Internal energy, $S$ [J kg$^{-1} K{^-1}$]")
        ax[2].set_ylabel(r"Internal energy [J kg$^{-1} $]", fontsize=14)
        # ax[1, 0].set_xlim(0, None)
        # ax[1, 0].set_ylim(0, None)

        ax[3].scatter(planet.A1_r / R_earth, planet.A1_T, s=0.1, c=colour)
        ax[3].set_xlabel(r"Radius $[R_\oplus]$", fontsize=14)
        ax[3].set_ylabel(r"Temperature [K]", fontsize=14)
        # ax[1, 1].set_xlim(0, None)
        # ax[1, 1].set_ylim(0, None)

        if not Till:
            planet.calculate_entropies()

            ax[4].scatter(planet.A1_r / R_earth, planet.A1_s, s=1, c=colour)
            ax[4].set_xlabel(r"Radius $[R_\oplus]$", fontsize=14)
            ax[4].set_ylabel(r"Entropy [J kg$^{-1} K{^-1}$]", fontsize=14)
            # ax[1, 0].set_xlim(0, None)
            # ax[1, 0].set_ylim(0, None)
    else:
        planet1 = woma.Planet(load_file=planet1loc)
        planet2 = woma.Planet(load_file=planet2loc)

        ax[0, 0].scatter(planet1.A1_r / R_earth, planet1.A1_rho, s=1)
        ax[0, 0].scatter(planet2.A1_r / R_earth, planet2.A1_rho, s=1)
        ax[0, 0].set_xlabel(r"Radius, $[R_\oplus]$")
        ax[0, 0].set_ylabel(r"Density, [kg m$^{-3}$]")

        ax[0, 1].scatter(planet1.A1_r / R_earth, planet1.A1_P, s=0.1)
        ax[0, 1].scatter(planet2.A1_r / R_earth, planet2.A1_P, s=0.1)
        ax[0, 1].set_xlabel(r"Radius, $[R_\oplus]$")
        ax[0, 1].set_ylabel(r"Pressure, [Pa]")
        ax[0, 1].set_yscale("log")

        ax[0, 2].scatter(planet1.A1_r / R_earth, planet1.A1_u, s=1)
        ax[0, 2].scatter(planet2.A1_r / R_earth, planet2.A1_u, s=1)
        ax[0, 2].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
        ax[0, 2].set_ylabel(r"Internal energy, $S$ [J kg$^{-1} $]")

        ax[0, 3].scatter(planet1.A1_r / R_earth, planet1.A1_T, s=0.1)
        ax[0, 3].scatter(planet2.A1_r / R_earth, planet2.A1_T, s=0.1)
        ax[0, 3].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
        ax[0, 3].set_ylabel(r"Temperature, $T$ [K]")

        if not Till:
            planet1.calculate_entropies()
            planet2.calculate_entropies()

            ax[0, 4].scatter(planet1.A1_r / R_earth, planet1.A1_s, s=1)
            ax[0, 4].scatter(planet2.A1_r / R_earth, planet2.A1_s, s=1)
            ax[0, 4].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
            ax[0, 4].set_ylabel(r"Entropy, $S$ [J kg$^{-1} K{^-1}$]")

    for axs in ax:
        axs.tick_params(axis="x", labelsize=11)  # Set font size for x-axis tick labels
        axs.tick_params(axis="y", labelsize=11)  # Set font size for y-axis tick labels

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.35)
    fig.savefig(
        "/user/home/qb20321/hammer/paper2_image/profile_m1d0.png",
        dpi=600,
        format="png",
        transparent=True,
        bbox_inches="tight",
    )
    # plt.show()

    outtime = 500


txtsize = 15
font_size = 20
idoff = 200000000
bodyoff = 100000000

params = {
    "axes.labelsize": font_size,
    "font.size": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "font.family": "serif",
}
matplotlib.rcParams.update(params)

# Material IDs ( = type_id * type_factor + unit_id )
type_factor = 100
type_ANEOS = 4
type_Til = 1
type_HHe = 2
type_SESAME = 3
type_idg = 0
id_body = 200000000
# Name and ID
Di_mat_id = {
    "ANEOS_iron": type_ANEOS * type_factor + 1,
    "ANEOS_iron_2": type_ANEOS * type_factor + 1 + id_body,
    "ANEOS_alloy": type_ANEOS * type_factor + 2,
    "ANEOS_alloy_2": type_ANEOS * type_factor + 2 + id_body,
    "ANEOS_forsterite": type_ANEOS * type_factor,
    "ANEOS_forsterite_2": type_ANEOS * type_factor + id_body,
    "HM80_HHE": type_HHe * type_factor,
    "HM80_HHE_2": type_HHe * type_factor + id_body,
    "HM80_ice": type_HHe * type_factor + 1,
    "HM80_ice_2": type_HHe * type_factor + 1 + id_body,
    "SS08_water": type_SESAME * type_factor + 3,
    "SS08_water_2": type_SESAME * type_factor + 3 + id_body,
    "AQUA": type_SESAME * type_factor + 4,
    "AQUA_2": type_SESAME * type_factor + 4 + id_body,
    "CMS19_HHe": type_SESAME * type_factor + 7,
    "idg_HHe": type_idg * type_factor,
    "Til_iron": type_Til * type_factor,
    "Til_iron_2": type_Til * type_factor + id_body,
    "Til_granite": type_Til * type_factor + 1,
    "Til_granite_2": type_Til * type_factor + 1 + id_body,
    "Til_water": type_Til * type_factor + 2,
    "Til_water_2": type_Til * type_factor + 2 + id_body,
    "Til_basalt": type_Til * type_factor + 3,
    "Til_basalt2": type_Til * type_factor + 3 + id_body,
}
# Colour
Di_mat_colour = {
    "ANEOS_iron": "tomato",  # np.array([237, 70, 14])
    "ANEOS_alloy": "tomato",
    "ANEOS_forsterite": "royalblue",  # np.array([14, 237, 59])
    "ANEOS_iron_2": "sandybrown",  # np.array([169, 237, 14])
    "ANEOS_alloy_2": "sandybrown",
    "ANEOS_forsterite_2": "pink",  # np.array([54, 14, 237])
    "SS08_water": "skyblue",
    "SS08_water_2": "skyblue",
    "AQUA": "skyblue",
    "AQUA_2": "skyblue",
    "HM80_HHE": "aliceblue",
    "HM80_HHE_2": "aliceblue",
    "HM80_ice": "skyblue",
    "HM80_ice_2": "skyblue",
    "CMS19_HHe": "lavenderblush",
    "idg_HHe": "lavenderblush",
    "Til_iron": "tomato",
    "Til_iron_2": "sandybrown",
    "Til_granite": "mediumseagreen",
    "Til_granite_2": "pink",
    "Til_water": "skyblue",
    "Til_water_2": "skyblue",
    "Til_basalt": "mediumseagreen",
    "Til_basalt2": "mediumseagreen",
}
# marker size
mks = 0.1
Di_mat_size = {
    "ANEOS_iron": mks,
    "ANEOS_alloy": mks,
    "ANEOS_forsterite": mks,
    "ANEOS_alloy_2": mks,
    "ANEOS_iron_2": mks,
    "ANEOS_forsterite_2": mks,
    "HM80_HHE": mks,
    "HM80_HHE_2": mks,
    "HM80_ice": mks,
    "HM80_ice_2": mks,
    "SS08_water": mks,
    "SS08_water_2": mks,
    "AQUA": mks,
    "AQUA_2": mks,
    "CMS19_HHe": mks,
    "idg_HHe": mks,
    "Til_iron": mks,
    "Til_iron_2": mks,
    "Til_granite": mks,
    "Til_granite_2": mks,
    "Til_water": mks,
    "Til_water_2": mks,
    "Til_basalt": mks,
    "Til_basalt2": mks,
}

Di_id_colour = {Di_mat_id[mat]: colour for mat, colour in Di_mat_colour.items()}

Di_id_size = {Di_mat_id[mat]: size * 3.5 for mat, size in Di_mat_size.items()}

# density limits
rhomin = 5e-5
rhomax = 15.0

# entropy limits
cmin = 1.5
cmax = 12.0

# number of cells for grid
Ng = 601j
Ngz = 21j

# region to grid/plot
xmax = 3
xmin = -xmax
ymin = -xmax
ymax = xmax
zmin = -3
zmax = 3

zcut = 6.0


def sw_plot(
    loc,
    snapshot_id=None,
    ax=None,
    npt=1e9,
    ax_lim=3.0,
    offcenter=None,
    outtime=None,
    midplane=False,
    rhocontour=False,
    quarter=False,
    belowZ=True,
    figsz=10,
    plotxy=True,
    selpid=None,
    EOS=None,
    mantleplot=True,
    plotwhole=False,
    pureplot=False,
    psize=5,
    PADX=10,
    PADY=10,
    tcx="w",
    tcy="w",
    SPH=True,
    load_region=None,
):
    """Select and load the particles to plot."""
    # Snapshot to load
    #     if snapshot_id is not None:
    #         snapshot = "snapshot_%04d.hdf5" % snapshot_id
    #         if load_region is not None:
    #             # mask = sw.mask(loc+snapshot)
    #             # boxhalf = mask.metadata.boxsize[0]/2
    # #            load_region *= R_earth # load_region is a 3x2 list with unit in R_earth (mks)
    # #            load_region[:,0] = boxhalf-load_region[:,0]
    #             load_region += boxhalf
    #             mask.constrain_spatial(load_region)
    #             data = sw.load(loc+snapshot,mask=mask)
    #         else:
    #             data = sw.load(loc+snapshot)
    #     else:
    #         if load_region is not None:
    #             mask = sw.mask(loc)
    # #            boxhalf = mask.metadata.boxsize[0]/2
    #             #load_region *= R_earth # load_region is a 3x2 list with unit in R_earth (mks)
    # #            load_region += boxhalf
    # #            print(load_region)
    #             mask.constrain_spatial(load_region)
    #             data = sw.load(loc,mask=mask)
    #         else:
    #             data = sw.load(loc)

    # Only load data with the axis limits and below z=0

    # box_mid = 0.5 * mask.metadata.boxsize[0].to(unyt.Rearth)
    # box_mid = 0.5 * mask.metadata.boxsize[0].to(unyt.m)

    # Load
    if snapshot_id is not None:
        data = sw.load(loc + snapshot)
    else:
        data = sw.load(loc)

    if SPH:
        box_mid = 0.5 * data.metadata.boxsize[0]
        pos = data.gas.coordinates - box_mid
        # box_mid = 0.5 * mask.metadata.boxsize[0].to(unyt.m)
        # pos = data.gas.coordinates.to(unyt.m) - box_mid
        id = data.gas.particle_ids
        mat_id = data.gas.material_ids.value
        m = data.gas.masses
    else:
        box_mid = 0.5 * data.metadata.boxsize[0].to(unyt.m)
        data.dark_matter.coordinates.convert_to_mks()
        pos = data.dark_matter.coordinates - box_mid
        pos = np.array(pos)
        id = data.dark_matter.particle_ids
        mat_id = np.ones(len(pos)) * 400
        m = data.dark_matter.masses

    if rhocontour:
        s = Snapshot()
        zi = np.linspace(zmin, zmax, int(Ngz.imag))
        X, Y, Z = np.mgrid[xmin:xmax:(Ng), ymin:ymax:(Ng), zmin:zmax:(Ngz)]
        XX, YY = np.mgrid[xmin:xmax:(Ng), ymin:ymax:(Ng)]

        cmap = plt.get_cmap("plasma").copy()
        cmap.set_under("w")

        data.gas.velocities.convert_to_cgs()
        vel = np.array(data.gas.velocities)

        data.gas.densities.convert_to_cgs()
        rho_cgs = np.array(data.gas.densities)

        s.pot = data.gas.potentials
        s.vel = np.array(vel)
        s.vx = np.array(vel[:, 0])
        s.vy = np.array(vel[:, 1])
        s.vz = np.array(vel[:, 2])
        s.x = np.array(pos[:, 0])
        s.y = np.array(pos[:, 1])
        s.z = np.array(pos[:, 2])
        s.rho = np.array(rho_cgs)

        modz = np.abs(s.z)

        vcut = 2 * s.vel.max()

        rhoi = scipy.interpolate.griddata(
            (
                s.x[(s.vx < vcut) * (modz < zcut)],
                s.y[(s.vx < vcut) * (modz < zcut)],
                s.z[(s.vx < vcut) * (modz < zcut)],
            ),
            s.rho[(s.vx < vcut) * (modz < zcut)],
            (X, Y, Z),
            method="linear",
            fill_value=1.0e-18,
        )

        coz = (s.z[s.pot == s.pot.min()])[0]  # s.z[modz<zcut]
        nn = (npy.nonzero(zi == (zi[zi <= coz])[-1])[0])[0]

        RHO_sh = rhoi[:, :, nn].T
        step = 0.8
        m = np.amax(RHO_sh)
        levels = np.arange(1.0, m, step)

    # pos_centerM = np.sum(pos * m[:,np.newaxis], axis=0) / np.sum(m)
    # pos -= pos_centerM
    sel_x = pos[:, 0] > 0
    pos_centerM = np.sum(pos[sel_x] * m[sel_x, np.newaxis], axis=0) / np.sum(m[sel_x])
    pos -= pos_centerM

    if midplane:
        sel = np.where(
            np.logical_and(pos[:, 2] < 0.05 * R_earth, pos[:, 2] > -0.05 * R_earth)
        )[0]
        mat_id = mat_id[sel]
        m = m[sel]
        id = id[sel]
        pos = pos[sel]

    # Restrict to z < 0
    elif belowZ:
        if plotxy:
            sel = np.where(
                np.logical_and(pos[:, 2] < 0.1 * R_earth, pos[:, 2] > -600 * R_earth)
            )[0]
        else:
            sel = np.where(
                np.logical_and(pos[:, 1] < 0.1 * R_earth, pos[:, 1] > -600 * R_earth)
            )[0]
        mat_id = mat_id[sel]
        id = id[sel]
        m = m[sel]
        pos = pos[sel]

    # Sort in z order so higher particles are plotted on top
    if plotxy:
        sort = np.argsort(pos[:, 2])
    else:
        sort = np.argsort(pos[:, 1])
    mat_id = mat_id[sort]
    id = id[sort]
    m = m[sort]
    pos = pos[sort]

    if npt > 0:
        mat_id[npt <= id] += id_body
    else:
        mat_id[np.logical_and(id <= 2 * idoff, id > idoff + bodyoff)] += id_body
        mat_id[np.logical_and(id <= idoff, id > bodyoff)] += id_body

    if EOS == "iron":
        core_sel = np.logical_or(
            mat_id == Di_mat_id["ANEOS_iron_2"], mat_id == Di_mat_id["ANEOS_iron"]
        )
        mantel_sel = np.logical_or(
            mat_id == Di_mat_id["ANEOS_forsterite_2"],
            mat_id == Di_mat_id["ANEOS_forsterite"],
        )
    elif EOS == "alloy":
        core_sel = np.logical_or(
            mat_id == Di_mat_id["ANEOS_alloy_2"], mat_id == Di_mat_id["ANEOS_alloy"]
        )
        mantel_sel = np.logical_or(
            mat_id == Di_mat_id["ANEOS_forsterite_2"],
            mat_id == Di_mat_id["ANEOS_forsterite"],
        )
    elif EOS == "Till":
        core_sel = np.logical_or(
            mat_id == Di_mat_id["Til_iron_2"], mat_id == Di_mat_id["Til_iron"]
        )
        mantle_sel = np.logical_or(
            mat_id == Di_mat_id["Til_granite_2"], mat_id == Di_mat_id["Til_granite"]
        )
    else:
        raise TypeError("Wrong EOS name, please check!")
    # impactor_mantle = mat_id==Di_mat_id["ANEOS_granite_2"]
    # target_mantle   = mat_id==Di_mat_id["ANEOS_granite"]

    # Plot the particles, coloured by their material.
    #    fig = plt.figure(figsize=(figsz, figsz))
    fig = plt.figure(figsize=(8, 8))
    fig.patch.set_facecolor("#111111")
    if ax is not None:
        ax = ax
    else:
        ax = plt.gca()

    ax.set_aspect("equal", anchor="C")
    ax.set_facecolor("k")

    # colour = np.zeros((len(pos),3),dtype=object)
    colour = np.empty(len(pos), dtype=object)
    for id_c, c in Di_id_colour.items():
        colour[mat_id == id_c] = c
    #    colour/=255.0

    size = np.empty(len(pos), dtype=float)
    for id_s, s in Di_id_size.items():
        size[mat_id == id_s] = s

    # print(np.unique(mat_id))
    # print(set(colour))
    # return pos, size, colour, impactor_core
    if plotwhole:
        ax.scatter(
            pos[:, 0],
            pos[:, 1],
            s=psize * size,
            c=colour,
            edgecolors="none",
            marker=".",
            alpha=1,
        )
        if selpid is not None:
            ax.scatter(
                pos[np.in1d(id, selpid), 0], pos[np.in1d(id, selpid), 1], s=1, c="lime"
            )
    else:
        if plotxy:
            if mantleplot:
                ax.scatter(
                    pos[mantel_sel, 0],
                    pos[mantel_sel, 1],
                    s=psize,
                    c=colour[mantel_sel],
                    edgecolors="none",
                    marker=".",
                    alpha=1,
                )
            ax.scatter(
                pos[core_sel, 0],
                pos[core_sel, 1],
                s=psize,
                c=colour[core_sel],
                edgecolors="none",
                marker=".",
                alpha=1,
            )

            if selpid is not None:
                ax.scatter(
                    pos[np.in1d(id, selpid), 0],
                    pos[np.in1d(id, selpid), 1],
                    s=1,
                    c="lime",
                )
        else:
            if mantleplot:
                ax.scatter(
                    pos[mantel_sel, 0],
                    pos[mantel_sel, 2],
                    s=4.5,
                    c=colour[mantel_sel],
                    edgecolors="none",
                    marker=".",
                    alpha=1,
                )
            ax.scatter(
                pos[core_sel, 0],
                pos[core_sel, 2],
                s=4.5,
                c=colour[core_sel],
                edgecolors="none",
                marker=".",
                alpha=1,
            )
            if selpid is not None:
                axx = ax.scatter(
                    pos[np.in1d(id, selpid), 0],
                    pos[np.in1d(id, selpid), 2],
                    s=1,
                    c="lime",
                )

    # ax.scatter(pos[target_mantle,0],pos[target_mantle,1],s=1,c='mediumseagreen',edgecolors="none")
    # ax.scatter(pos[impactor_core,0],pos[impactor_core,1],s=1,c='sandybrown',edgecolors="none")
    #     ax.scatter(pos[target_core,0],pos[target_core,1],s=1,c='tomato',edgecolors="none")
    if rhocontour:
        con = ax.contour(
            XX,
            YY,
            RHO_sh,
            levels,
            origin="lower",
            extent=[-ax_lim, ax_lim, -ax_lim, ax_lim],
            alphas=0,
            cmap="YlGn",
        )

    # ax.set_title('0.897 $M_\oplus$ Protoearth with 1.2% $M_{total}$ atmosphere ')
    if offcenter is not None:
        ax.set_xlim(-(ax_lim + offcenter), ax_lim - offcenter)
        ax.set_yticks(ax.get_xticks())
        ax.set_ylim(-ax_lim, ax_lim)

    ax.yaxis.label.set_color("w")
    ax.xaxis.label.set_color("w")
    ax.tick_params(axis="x", colors=tcx, labelsize=14, direction="in", pad=PADX)
    ax.tick_params(axis="y", colors=tcy, labelsize=14, direction="in", pad=PADY)

    if plotxy:
        ax.set_ylabel(r"y Position ($R_\oplus$)", fontsize=22, c="w")
    else:
        ax.set_ylabel(r"z Position ($R_\oplus$)", fontsize=22, c="w")

    ax.set_xlabel(r"x Position ($R_\oplus$)", fontsize=22, c="w")

    #    ax.tick_params(axis="y",direction="in")

    if quarter == 1:
        ax.set_xlim(0, ax_lim)
        ax.set_ylim(0, ax_lim)
    elif quarter == 2:
        ax.set_xlim(0, ax_lim)
        ax.set_ylim(-ax_lim, 0)
    elif quarter == 3:
        ax.set_xlim(-ax_lim, 0)
        ax.set_ylim(-ax_lim, 0)
    elif quarter == 4:
        ax.set_xlim(-ax_lim, 0)
        ax.set_ylim(0, ax_lim)
    else:
        if quarter != 0:
            raise TypeError("Wrong quarter numer, choose from 1, 2, 3, 4")

    if outtime is not None:
        # ax.text(0.97,0.97,r'v = 9.17km/s, b = 0.7            {0:6.2f} h'.format(snapshot_id*outtime/3600),
        #              transform=ax.transAxes,ha='right',va='top',fontsize=txtsize,color='w')
        ax.text(
            0.97,
            0.97,
            r"{0:6.2f} h".format(outtime / 3600),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=18,
            color="w",
            fontweight="bold",
        )

    if pureplot:
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.set_xlim(ax_lim[0], ax_lim[1])
        ax.set_ylim(ax_lim[2], ax_lim[3])

        scalebar = ScaleBar(
            1,
            "m",
            fixed_value=0.5 * Bound.R_earth,
            scale_formatter=lambda value, unit: f"{value/Bound.R_earth}"
            + r"$R_\oplus$",
            color="w",
            frameon=False,
            location="lower left",
            font_properties={"size": 14},
        )
        ax.add_artist(scalebar)

        # ax.xaxis.set_ticks([-0.5*ax_lim,0.0,0.5*ax_lim])
        # ax.yaxis.set_ticks([-0.5*ax_lim,0.0,0.5*ax_lim])

        # ax.set_xticks([0,0.25,0.5,0.75,1.])
        # ax.set_yticks([0,0.25,0.5,0.75,1.])

        # ax.get_xaxis().set_ticks([])
        # ax.get_yaxis().set_ticks([])
        ax.axis("off")

        fig.savefig(
            "/user/home/qb20321/hammer/denseplanets/figures/select_particle.png",
            format="png",
            dpi=300,
            facecolor="k",
            bbox_inches="tight",
        )

    # plt.tight_layout()
    # plt.show()
    # # save = loc + "plot/snaplot_%04d.png" % snapshot_id
    # # plt.savefig(save)
    # plt.cla()
    # plt.clf()
    # plt.close()
    return ax


def print_impact_info(loc_tar, loc_imp, Xves=1.0, verbose=1):

    pos_tar, vel_tar, h_tar, m_tar, rho_tar, p_tar, u_tar, matid_tar, R_tar = (
        loadsw_to_woma(loc_tar)
    )
    pos_imp, vel_imp, h_imp, m_imp, rho_imp, p_imp, u_imp, matid_imp, R_imp = (
        loadsw_to_woma(loc_imp)
    )

    M_t = np.sum(m_tar)
    M_i = np.sum(m_imp)
    R_t = R_tar
    R_i = R_imp
    gammma = M_i / M_t
    # Mutual escape speed
    v_esc = np.sqrt(2 * G * (M_t + M_i) / (R_t + R_i))
    times = Xves
    if verbose:
        print("Total mass    = %.5f" % ((M_t + M_i) / M_earth))
        print(
            "Target mass   = %.5f, number of partiles %d, Radius=%.5f"
            % (M_t / M_earth, len(pos_tar), R_tar / R_earth)
        )
        print(
            "Impactor mass = %.5f, number of partiles %d, Radius=%.5f"
            % (M_i / M_earth, len(pos_imp), R_imp / R_earth)
        )
        print("Mass ratio    = %.4f" % (gammma))
        print("Mutual escape velocity = %f m/s" % v_esc)
        print("%.3f times the mutual escape velocity = %f m/s" % (times, v_esc * times))

    return v_esc, gamma
