import numpy as np
import pandas as pd
import swiftsimio as sw
import unyt
import os

G_cgs = 6.67408e-8  # in cgs
Mearth_cgs = 5.97240e27  # in cgs


def loadsw_to_woma(snapshot, unit="mks", if_R_atmos=False, if_core_radius=False):
    """load swift hdf5 snapshot date and calculate some necessary variables

    Args:
        snapshot (str, manditary): the path of the snapshot
        unit (str, optional): are we load mks unit or cgs unit. Defaults to "mks".
        R_atmos (bool, optional): Do we count thickness of the atmosphere when
            calculating the radius of the planet. Defaults to False.
        if_core_radius (bool, optional): Whether to return the core radius. Defaults to False.


    """
    # Load
    data = sw.load(snapshot)
    if unit == "mks":
        box_mid = 0.5 * data.metadata.boxsize[0].to(unyt.m)
        data.gas.coordinates.convert_to_mks()
        pos = np.array(data.gas.coordinates - box_mid)
        data.gas.velocities.convert_to_mks()
        vel = np.array(data.gas.velocities)
        data.gas.smoothing_lengths.convert_to_mks()
        h = np.array(data.gas.smoothing_lengths)
        data.gas.masses.convert_to_mks()
        m = np.array(data.gas.masses)
        data.gas.densities.convert_to_mks()
        rho = np.array(data.gas.densities)
        data.gas.pressures.convert_to_mks()
        p = np.array(data.gas.pressures)
        data.gas.internal_energies.convert_to_mks()
        u = np.array(data.gas.internal_energies)
        matid = np.array(data.gas.material_ids)
        # pid     = np.array(data.gas.particle_ids)

    elif unit == "cgs":
        box_mid = 0.5 * data.metadata.boxsize[0].to(unyt.cm)
        data.gas.coordinates.convert_to_cgs()
        pos = np.array(data.gas.coordinates - box_mid)
        data.gas.velocities.convert_to_cgs()
        vel = np.array(data.gas.velocities)
        data.gas.smoothing_lengths.convert_to_cgs()
        h = np.array(data.gas.smoothing_lengths)
        data.gas.masses.convert_to_cgs()
        m = np.array(data.gas.masses)
        data.gas.densities.convert_to_cgs()
        rho = np.array(data.gas.densities)
        data.gas.pressures.convert_to_cgs()
        p = np.array(data.gas.pressures)
        data.gas.internal_energies.convert_to_cgs()
        u = np.array(data.gas.internal_energies)
        matid = np.array(data.gas.material_ids)
        # pid     = np.array(data.gas.particle_ids)
    else:
        raise TypeError("Wrong unit selection, please check!!")

    pos_centerM = np.sum(pos * m[:, np.newaxis], axis=0) / np.sum(m)
    vel_centerM = np.sum(vel * m[:, np.newaxis], axis=0) / np.sum(m)

    pos -= pos_centerM
    vel -= vel_centerM

    core_key_list = np.array([401, 402])
    atmos_key_list = np.array([0, 1, 2, 200, 305, 306, 307])
    uniq_mat = np.unique(matid)
    atmos_id = np.intersect1d(atmos_key_list, uniq_mat)
    core_id = np.intersect1d(core_key_list, uniq_mat)[0]

    # if atmos_id or core_id is empty, set to np.nan
    if len(atmos_id) == 0:
        atmos_id = np.nan

    if not if_R_atmos:
        pos_to_use = pos[matid != atmos_id]
        matid = matid[matid != atmos_id]
    else:
        pos_to_use = np.squeeze(pos)[0]

    pos_to_use = np.squeeze(pos_to_use)
    matid = np.squeeze(matid)

    xy = np.hypot(pos_to_use[:, 0], pos_to_use[:, 1])
    r = np.hypot(xy, pos_to_use[:, 2])
    if if_core_radius:
        R_core = np.mean(np.sort(r[matid == core_id])[-100:])
    r = np.sort(r)
    R = np.mean(r[-200:])

    if if_core_radius:
        return pos, vel, h, m, rho, p, u, matid, R, R_core
    else:
        return pos, vel, h, m, rho, p, u, matid, R


def edacm(
    R=0.0,  # radius of the target in cgs
    r=0.0,  # radius of the impactor in cgs
    b=None,
    v=None,  # Impact velocity in cgs
    M_tar=0.0,  # in csg
    M_tot=0.0,  # in csg
    loc_tar=None,
    loc_imp=None,
    if_R_atmos=False,
):
    if (loc_tar is not None) and (loc_imp is not None):
        _, _, _, m_tar, _, _, _, _, R_tar = loadsw_to_woma(
            loc_tar, unit="cgs", if_R_atmos=if_R_atmos
        )
        _, _, _, m_imp, _, _, _, _, R_imp = loadsw_to_woma(
            loc_imp, unit="cgs", if_R_atmos=if_R_atmos
        )

        M_tar = np.sum(m_tar)
        M_imp = np.sum(m_imp)
        M_tot = M_tar + M_imp
        R = R_tar
        r = R_imp
    if (b == None) or (v == None):
        raise ValueError("Please provide a valid velocity and impact parameters")

    c_star = 1.9  # +-0.3
    rho1 = 1.0  # g * cm-3
    mu_bar = 0.36  # +-0.01
    M_imp = M_tot - M_tar
    gamma = M_imp / M_tar
    mu = gamma * M_tar / (1 + gamma)
    # angle correction
    l = (R + r) * (1 - b)
    alpha = (3 * r * l**2 - l**3) / (4 * r**3)

    mu_alpha = (alpha * M_imp * M_tar) / (alpha * M_imp + M_tar)
    R_c1 = np.power(3 * (1 + gamma) * M_tar / (np.pi * 4), 1 / 3)  # in cm
    Q_RD_star_gamma1 = c_star * (4 * np.pi / 5) * rho1 * G_cgs * R_c1**2
    Q_RD_star = Q_RD_star_gamma1 * ((1 + gamma) ** 2 / (4 * gamma)) ** (
        2 / (3 * mu_bar) - 1
    )
    Q_RD_star_prime = Q_RD_star * (mu / mu_alpha) ** (2 - (3 * mu_bar) / 2)

    Q_R = 0.5 * mu * v**2 / M_tot

    Q_R_norm = Q_R / Q_RD_star_prime

    return (
        Q_R,
        Q_R_norm,
        Q_RD_star_prime,
        M_tot / Mearth_cgs,
        M_tar / Mearth_cgs,
        M_imp / Mearth_cgs,
    )


def load_melt_vapor_curve(mat_id):
    this_dir, _ = os.path.split(__file__)
    if mat_id == 400:
        dataloc_meltCurve = os.path.join(this_dir, "data/s19_forsterite_meltCurve.csv")
        dataloc_vapourCurve = os.path.join(
            this_dir, "data/s19_forsterite_vaporCurve.csv"
        )
    elif mat_id == 401:
        dataloc_meltCurve = os.path.join(this_dir, "data/s20_iron_meltCurve.csv")
        dataloc_vapourCurve = os.path.join(this_dir, "data/s20_iron_vaporCurve.csv")
    elif mat_id == 402:
        dataloc_meltCurve = os.path.join(this_dir, "data/s20_alloy_meltCurve.csv")
        dataloc_vapourCurve = os.path.join(this_dir, "data/s20_alloy_vaporCurve.csv")
    else:
        raise ValueError(
            "Currently only have SESAME iron (401), SESAME alloy (402) and SESAME forsterite (400) vapour curve data"
        )
    data_meltCurve = pd.read_csv(dataloc_meltCurve, index_col=False)
    data_vaporCurve = pd.read_csv(dataloc_vapourCurve, index_col=False)

    return data_meltCurve, data_vaporCurve


class VapourFrc:
    iron_critical_P = 0.658993  # iron critical point Pressue in Gpa
    forsterite_critical_P = 0.159304  # forsterite critical point Pressue in Gpa

    iron_critical_S = 3.78634  # iron critical point
    forsterite_critical_S = 6.37921  # forsterite critical point entropy in KJ/kg/K

    def __init__(self, mat_id, entropy, pressure):
        """initialization variables needed to calculated the vapour fraction.

        Args:
            mat_id (int): material id
            entropy (float): in mks J/kg/K
            pressure (float): in mks pa

        Raises:
            ValueError: _description_
        """
        if mat_id not in [400, 401, 402]:
            raise ValueError(
                "Currently only have iron and forsterite vapour curve loaded"
            )
        self.mat_id = mat_id
        entropy *= 1e-3  # switch to kJ/kg/K
        pressure *= 1e-9  # switch to Gpa
        self.vapour_frac = np.zeros(len(entropy))
        self.super_frac = np.zeros(len(entropy))
        if mat_id == 400:
            self.super_sel = np.logical_and(
                pressure > VapourFrc.forsterite_critical_P,
                entropy > VapourFrc.forsterite_critical_S,
            )
            self.sel = pressure < VapourFrc.forsterite_critical_P
            self.entropy = entropy[self.sel]
            self.pressure = pressure[self.sel]
        elif mat_id == 401:
            self.super_sel = np.logical_and(
                pressure > VapourFrc.iron_critical_P,
                entropy > VapourFrc.iron_critical_S,
            )
            self.sel = pressure < VapourFrc.iron_critical_P
            self.entropy = entropy[self.sel]
            self.pressure = pressure[self.sel]
        else:
            raise ValueError(
                "Currently only have iron and forsterite vapour curve loaded"
            )
        self.super_frac[self.super_sel] = (
            1  # set supercritical vapour fraction factor to 2
        )
        self.meltCurve, self.vaporCurve = load_melt_vapor_curve(mat_id)

    def lever(self, left_point, right_point, s):
        """given the liquid side vapor curve entropy and the vapour side vapor curve entropy,
        calculate the fraction of vapour using lever-rule.
        """

        v_frac = np.zeros(len(self.entropy))
        v_frac = (s - left_point) / (right_point - left_point)
        v_frac[v_frac < 0] = 0
        v_frac[v_frac > 1] = 1

        return v_frac

    def vapour_fraction(self):
        liquid_side_entropy = np.interp(
            self.pressure,
            self.vaporCurve["P_vapor_curve_liquid"],
            self.vaporCurve["S_vapor_curve_liquid"],
        )
        vapour_side_entropy = np.interp(
            self.pressure,
            self.vaporCurve["P_vapor_curve_gas"],
            self.vaporCurve["S_vapor_curve_gas"],
        )

        vapour_frac = self.lever(
            liquid_side_entropy * 1e3, vapour_side_entropy * 1e3, self.entropy
        )
        self.vapour_frac[self.sel] = vapour_frac
        return self.vapour_frac

    def super_critical(self):
        return self.super_sel


class PhaseFinder:
    """Find the phase of the material based on the entropy and pressure
    phase 1: solid
    phase 2: liquid+solid
    phase 3: liquid
    phase 4: liquid+vapour
    phase 5: vapour
    phase 6: supercritical

    Phase number is stored as an integer array.
    """

    iron_critical_P = 0.658993  # iron critical point Pressue in Gpa
    forsterite_critical_P = 0.159304  # forsterite critical point Pressue in Gpa

    iron_critical_S = 3.78634  # iron critical point
    forsterite_critical_S = 6.37921  # forsterite critical point entropy in KJ/kg/K

    def __init__(self, mat_id, entropy, pressure):
        if mat_id not in [400, 401, 402]:
            raise ValueError(
                "Currently only have iron and forsterite vapour curve loaded"
            )
        self.mat_id = mat_id
        self.entropy *= 1e-3  # switch to kJ/kg/K
        self.pressure *= 1e-9  # switch to Gpa
        self.phase = np.zeros(len(entropy))
        if mat_id == 400:
            self.critical_P = 0.159304
            self.critical_S = 6.37921
            self.phase[
                np.logical_and(
                    pressure > self.critical_P,
                    entropy > self.critical_S,
                )
            ] = 6

        elif mat_id == 401:
            critical_P = 0.658993
            critical_S = 3.78634

            self.phase[
                np.logical_and(
                    pressure > self.critical_P,
                    entropy > self.critical_S,
                )
            ] = 6

        else:
            raise ValueError(
                "Currently only have iron and forsterite vapour curve loaded"
            )

        self.meltCurve, self.vaporCurve = load_melt_vapor_curve(mat_id)

    # def phase_finder(self):
    #     # deal with material with pressure higher than critical point but not supercritical

    #     meltC_liquid_side_entropy = np.interp(
    #         self.pressure[(self.pressure>self.critical_P)&(self.entropy<self.critical_S)],
    #         self.vaporCurve["P_vapor_curve_liquid"],
    #         self.vaporCurve["S_vapor_curve_liquid"],
    #     )
    #     meltC_vapour_side_entropy = np.interp(
    #         self.pressure[],
    #         self.vaporCurve["P_vapor_curve_gas"],
    #         self.vaporCurve["S_vapor_curve_gas"],
    #     )

    #     if self.en


def main():
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d

    PVsl = load_PSvc_data(mat_id=401)

    interp_func = interp1d(
        np.flip(PVsl[2]), np.flip(PVsl[0]), kind="linear", fill_value="extrapolate"
    )
    s_array = np.linspace(np.min(PVsl[2]), np.max(PVsl[2]), 100)
    p_array = interp_func(s_array)

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(1e3 * s_array, np.log10(p_array), c="b")
    ax.scatter(1e3 * PVsl[2], np.log10(PVsl[0]), c="r", s=1)
    ax.scatter(1e3 * PVsl[2][0], np.log10(PVsl[0][0]), c="lime", s=20)

    plt.show()
    # pass


if __name__ == "__main__":
    main()
