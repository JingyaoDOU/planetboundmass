import numpy as np
import swiftsimio as sw
import unyt
import os

G_cgs = 6.67408e-8  # in cgs
Mearth_cgs = 5.97240e27  # in cgs


def loadsw_to_woma(snapshot, unit="mks", if_R_atmos=False):
    """load swift hdf5 snapshot date and calculate some necessary variables

    Args:
        snapshot (str, manditary): the path of the snapshot
        unit (str, optional): are we load mks unit or cgs unit. Defaults to "mks".
        R_atmos (bool, optional): Do we count thickness of the atmosphere when
            calculating the radius of the planet. Defaults to False.


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

    atmos_key_list = np.array([0, 1, 2, 200, 305, 306, 307])
    uniq_mat = np.unique(matid)
    atmos_id = np.intersect1d(atmos_key_list, uniq_mat)

    if not if_R_atmos:
        pos_to_use = pos[matid != atmos_id]
    else:
        pos_to_use = pos

    xy = np.hypot(pos_to_use[:, 0], pos_to_use[:, 1])
    r = np.hypot(xy, pos_to_use[:, 2])
    r = np.sort(r)
    R = np.mean(r[-200:])

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


def load_PSvc_data(mat_id):
    this_dir, _ = os.path.split(__file__)
    if mat_id == 400:
        dataloc = os.path.join(this_dir, "data/s19_forsterite_vc.txt")
    elif mat_id == 401:
        dataloc = os.path.join(this_dir, "data/s20_iron_vc.txt")
    elif mat_id == 402:
        dataloc = os.path.join(this_dir, "data/s20_alloy_vc.txt")
    else:
        raise ValueError(
            "Currently only have SESAME iron (401), SESAME alloy (402) and SESAME forsterite (400) vapour curve data"
        )
    data_PVsl = np.loadtxt(dataloc)

    return data_PVsl


# Fors_PVsl  = load_PSvc_data(mat_id=400)
# Iron_PVsl  = load_PSvc_data(mat_id=401)
# Alloy_PVsl = load_PSvc_data(mat_id=402)
class VapourFrc:

    iron_trip = 0.658993  # iron triple point
    forsterite_trip = 0.159304  # forsterite triple point

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

        if mat_id == 400:
            self.entropy = entropy[pressure < VapourFrc.forsterite_trip]
            self.pressure = pressure[pressure < VapourFrc.forsterite_trip]
        elif mat_id == 401:
            self.entropy = entropy[pressure < VapourFrc.iron_trip]
            self.pressure = pressure[pressure < VapourFrc.iron_trip]
        else:
            raise ValueError(
                "Currently only have iron and forsterite vapour curve loaded"
            )
        self.PVsl = load_PSvc_data(mat_id)

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
            self.pressure, np.flip(self.PVsl[0]), np.flip(self.PVsl[2])
        )
        vapour_side_entropy = np.interp(
            self.pressure, np.flip(self.PVsl[1]), np.flip(self.PVsl[3])
        )

        vapour_frac = self.lever(
            liquid_side_entropy * 1e3, vapour_side_entropy * 1e3, self.entropy
        )
        self.vapour_frac = vapour_frac
        return vapour_frac


def main():
    pass


if __name__ == "__main__":
    main()
