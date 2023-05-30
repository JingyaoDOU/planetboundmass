def gadget_plot(
    loc,
    ax=None,
    plotxy=True,
    ax_lim=2.0,
    figsz=6,
    offcenter=0.0,
    sel_id=None,
    plotmantle=True,
):
    s = Snapshot()
    s.load(loc, center=True)

    if plotxy:
        sel = np.where(np.logical_and(s.z / Rearth < 0.1, s.z / Rearth > -6))[0]
    else:
        sel = np.where(np.logical_and(s.y / Rearth < 0.1, s.y / Rearth > -6))[0]
    id = s.id[sel]
    x = s.x[sel] / Rearth
    y = s.y[sel] / Rearth
    z = s.z[sel] / Rearth

    if plotxy:
        sort = np.argsort(z)
    else:
        sort = np.argsort(y)
    id = id[sort]
    x = x[sort]
    z = z[sort]
    y = y[sort]

    xm1 = x[np.logical_and(id > IDOFF, id <= IDOFF + BODYOFF)]
    ym1 = y[np.logical_and(id > IDOFF, id <= IDOFF + BODYOFF)]
    zm1 = z[np.logical_and(id > IDOFF, id <= IDOFF + BODYOFF)]

    xm2 = x[np.logical_and(id <= 2 * IDOFF, id > IDOFF + BODYOFF)]
    ym2 = y[np.logical_and(id <= 2 * IDOFF, id > IDOFF + BODYOFF)]
    zm2 = z[np.logical_and(id <= 2 * IDOFF, id > IDOFF + BODYOFF)]

    xc1 = x[np.logical_and(id <= IDOFF, id <= BODYOFF)]
    yc1 = y[np.logical_and(id <= IDOFF, id <= BODYOFF)]
    zc1 = z[np.logical_and(id <= IDOFF, id <= BODYOFF)]

    xc2 = x[np.logical_and(id <= IDOFF, id > BODYOFF)]
    yc2 = y[np.logical_and(id <= IDOFF, id > BODYOFF)]
    zc2 = z[np.logical_and(id <= IDOFF, id > BODYOFF)]

    fig = plt.figure(figsize=(figsz, figsz))
    if ax is not None:
        ax = ax
    else:
        ax = plt.gca()

    ax.set_aspect("equal")

    s = 3.5
    if plotxy:
        if plotmantle:
            ax.scatter(
                xm1,
                ym1,
                s=s,
                c="mediumseagreen",
                edgecolors="none",
                marker=".",
                alpha=1,
            )
            ax.scatter(xm2, ym2, s=s, c="pink", edgecolors="none", marker=".", alpha=1)
        ax.scatter(xc1, yc1, s=s, c="tomato", edgecolors="none", marker=".", alpha=1)
        ax.scatter(
            xc2, yc2, s=s, c="sandybrown", edgecolors="none", marker=".", alpha=1
        )
        if sel_id is not None:
            ax.scatter(
                x[np.in1d(id, sel_id)],
                y[np.in1d(id, sel_id)],
                s=15,
                c="lime",
                edgecolors="none",
                marker=".",
                alpha=1,
            )
    else:
        if plotmantle:
            ax.scatter(
                xm1,
                zm1,
                s=s,
                c="mediumseagreen",
                edgecolors="none",
                marker=".",
                alpha=1,
            )
            ax.scatter(xm2, zm2, s=s, c="pink", edgecolors="none", marker=".", alpha=1)
        ax.scatter(xc1, zc1, s=s, c="tomato", edgecolors="none", marker=".", alpha=1)
        ax.scatter(
            xc2, zc2, s=s, c="sandybrown", edgecolors="none", marker=".", alpha=1
        )
        if sel_id is not None:
            ax.scatter(
                x[np.in1d(id, sel_id)],
                z[np.in1d(id, sel_id)],
                s=15,
                c="lime",
                edgecolors="none",
                marker=".",
                alpha=1,
            )

    ax.set_xlim(-(ax_lim + offcenter), ax_lim - offcenter)
    ax.set_yticks(ax.get_xticks())
    ax.set_ylim(-ax_lim, ax_lim)
    if plotxy:
        ax.set_ylabel(r"y Position ($R_\oplus$)")
    else:
        ax.set_ylabel(r"z Position ($R_\oplus$)")

    ax.set_xlabel(r"x Position ($R_\oplus$)")

    ax.yaxis.label.set_color("w")
    ax.xaxis.label.set_color("w")
    ax.tick_params(axis="x", colors="w", labelsize=14)
    ax.tick_params(axis="y", colors="w", labelsize=14)

    return ax


pid[pid_tar_core_sel] = np.arange(1, np.sum(pid_tar_core_sel) + 1)
pid[pid_tar_mantle_sel] = np.arange(1, np.sum(pid_tar_mantle_sel) + 1) + IDOFF

pid[pid_imp_core_sel] = np.arange(1, np.sum(pid_imp_core_sel) + 1) + BODYOFF
pid[pid_imp_mantle_sel] = np.arange(1, np.sum(pid_imp_mantle_sel) + 1) + IDOFF + BODYOFF
