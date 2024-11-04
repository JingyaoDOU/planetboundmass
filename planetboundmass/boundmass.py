import numpy as np
import matplotlib.pyplot as plt
import swiftsimio as sw
import h5py
import unyt
import random
from copy import deepcopy, copy
import planetboundmass.sw_planet_tools as swtools
from planetboundmass.sw_planet_tools import VapourFrc
import pandas as pd
import woma
import seaborn as sns


class Bound:
    G = 6.67408e-11  # m^3 kg^-1 s^-2
    M_earth = 5.97240e24  # kg
    R_earth = 6.371e6  # m
    id_body = 200000000

    def __init__(
        self,
        filename=None,
        minibound=2000,
        num_rem=1,
        total_mass=1e9,
        npt=1e9,
        m_tar=None,
        b=None,
        v=None,
        tolerance=1e-10,
        maxit=1000,
        max_bad_seeds=1000,
        verbose=1,
        hit_direction="pX",
    ):
        if not filename.endswith("hdf5"):
            raise TypeError("Wrong filename, please check the file extension")
        self.filename = filename

        self.filename_check()
        if not self.standard_filename:
            if verbose:
                print("not standard filename")
            # if eosid_list is None:
            #     raise ValueError(
            #         "Non-standard filename, please provide an eos id dictionary"
            #     )
            self.total_mass = total_mass
            self.npt = npt
            # self.ironid = eosid_list["ironid"]
            # self.siid = eosid_list["siid"]
            # self.waterid = eosid_list["waterid"]
            # self.atmosid = eosid_list["atmosid"]
            self.m_tar = m_tar
            self.b = b
            self.v = v
            self.attri = None

        self.num_rem = num_rem
        self.minibound = minibound
        self.max_bad_seeds = max_bad_seeds
        self.tolerance = tolerance
        self.verbose = verbose
        self.snap = Snap(filename, npt=self.npt)
        self.hit_direction = hit_direction
        self.load_data()

    def filename_check(self):
        prefix_name = self.filename.split("/")[-1].split("_")[0]
        if (prefix_name == "snapOUT") and (
            len(self.filename.split("/")[-1].split("_")) > 9
        ):
            variable_list = self.filename.split("/")[-1].split("_")[1:]
            self.m_tar = float(variable_list[2].replace("d", "."))
            self.npt = int(variable_list[3][3:])
            self.total_mass = float(variable_list[4].replace("d", "."))
            self.v = float(variable_list[5][1:].split("k")[0].replace("d", "."))
            self.b = float(variable_list[6].replace("d", ".")[1:])
            self.hit_direction = variable_list[7]
            self.attri = ("_".join(variable_list[8:])).split(".")[0]

            self.standard_filename = True
        else:
            self.standard_filename = False
        #

    def load_data(self):
        # load all the data from self.snap
        self.boxsize = self.snap.boxsize
        self.pos = self.snap.pos
        self.vel = self.snap.vel
        self.rho_mks = self.snap.rho_mks
        self.m = self.snap.m
        self.pid = self.snap.pid
        self.matid = self.snap.matid
        self.pot = self.snap.pot
        self.unique_matid = self.snap.unique_matid
        self.u = self.snap.u
        self.p_mks = self.snap.p_mks
        self.matid_tar_imp = self.snap.matid_tar_imp

        self.Di_id_colour = self.snap.Di_id_colour
        self.Di_id_size = self.snap.Di_id_size
        self.Di_id_mat = self.snap.Di_id_mat
        self.iron_key_list = self.snap.iron_key_list
        self.si_key_list = self.snap.si_key_list

    def find_bound(self, reU=False):
        """Find the bound particles and calculate the bound remnants.

        Args:
            reU (bool, optional): Whether to recalculate the potential energy of particles about a remnant was found. Defaults to False.
        """
        # find bound particles
        bad_seeds = 0
        remnant_id = 1  # intialization

        bound = np.zeros(len(self.pid))
        bound_id = np.zeros(self.num_rem)
        m_rem = np.zeros(self.num_rem)
        num_par_rem = np.zeros(
            self.num_rem, dtype=int
        )  # Numbe of particles for each remnant
        mass_ratio = np.zeros(self.num_rem)  # M_rem / M_total

        element_ratio_array = {}
        element_mass_array = {}

        for mat_id in self.unique_matid:
            array_name = self.Di_id_mat[mat_id] + "_ratio"
            element_ratio_array[array_name] = np.zeros(self.num_rem)
            array_name = self.Di_id_mat[mat_id] + "_mass"
            element_mass_array[array_name] = np.zeros(self.num_rem)

        while True:
            if np.sum(bound == 0) < self.minibound:
                if self.verbose:
                    print("----------break------------")
                    print("No enough particles left to be count as a remnant")
                break

            if bad_seeds > self.max_bad_seeds:
                if self.verbose:
                    print("----------break------------")
                    print("Bad seeds larger than the maximum allowed bad seeds")
                break

            unbound_pid = self.pid[bound == 0]
            unbound_pot = self.pot[bound == 0]
            unbound_pos = self.pos[bound == 0]

            if reU:
                unbound_m = self.m[bound == 0]
                bound_pot = self.pot[bound > 0]
                bound_pos = self.pos[bound > 0]
                bound_m = self.m[bound > 0]
                if remnant_id > 1:  # already found one remnant
                    # re-calculating the potential for the unbound particles
                    for i in range(len(unbound_pid)):
                        unbound_pot[i] -= (unyt.J / unyt.kg) * np.sum(
                            -Bound.G
                            * bound_m
                            / np.hypot(
                                bound_pos[:, 2] - unbound_pos[i, 2],
                                np.hypot(
                                    bound_pos[:, 0] - unbound_pos[i, 0],
                                    bound_pos[:, 1] - unbound_pos[i, 1],
                                ),
                            )
                        )

            arg_init_min_potseed = np.argmin(unbound_pot)
            init_min_pot_pid = unbound_pid[arg_init_min_potseed]

            arg_init_min_potseed = np.where(np.in1d(self.pid, init_min_pot_pid))[0]

            bound[arg_init_min_potseed] = remnant_id

            bnd_m = np.squeeze(self.m[arg_init_min_potseed])

            bnd_pos = np.squeeze(self.pos[arg_init_min_potseed])

            bnd_vel = np.squeeze(self.vel[arg_init_min_potseed])

            oldm = bnd_m / 10.0
            count = 0
            goback = False
            a = 1
            b = 1
            maxit = 1000
            # raise TypeError("check")
            while (count <= maxit) & (np.abs(oldm - bnd_m) / oldm > self.tolerance):
                oldm = bnd_m
                sel = np.where(bound == 0)[0]
                pid_tmp = self.pid[sel]
                # compute kinetic velocities
                ke = 0.5 * self.m[sel] * np.sum((self.vel[sel] - bnd_vel) ** 2, axis=1)
                pe = (
                    -Bound.G
                    * bnd_m
                    * self.m[sel]
                    / np.hypot(
                        self.pos[sel, 2] - bnd_pos[2],
                        np.hypot(
                            self.pos[sel, 0] - bnd_pos[0], self.pos[sel, 1] - bnd_pos[1]
                        ),
                    )
                )

                sel_bound = ke + pe < 0.0

                if (count == 0) and (np.sum(sel_bound) == 0):
                    bound[bound == remnant_id] = (
                        -1
                    )  # through the bad seeds with remnant id -1
                    bad_seeds += 1
                    goback = True
                    # print("Bad starting")
                    break

                if np.sum(sel_bound) > 0:
                    pid_bnd = pid_tmp[sel_bound]
                    arg_bound_in = np.where(np.in1d(self.pid, pid_bnd))[0]
                    bound[arg_bound_in] = remnant_id
                    bnd_m = np.sum(self.m[arg_bound_in])
                    bnd_pos = (
                        np.sum(
                            self.pos[arg_bound_in] * self.m[arg_bound_in, np.newaxis],
                            axis=0,
                        )
                        / bnd_m
                    )
                    bnd_vel = (
                        np.sum(
                            self.vel[arg_bound_in] * self.m[arg_bound_in, np.newaxis],
                            axis=0,
                        )
                        / bnd_m
                    )

                count += 1

            if goback:
                continue

            numbound = np.sum(bound == remnant_id)

            if numbound < self.minibound:
                bound[bound == remnant_id] = -1
                # print("Not enough particles in the bound group")
                continue

            arg_bound_out = bound == remnant_id
            m_bound = self.m[arg_bound_out]
            rem_mass = np.sum(m_bound)  # mass of the remnant in this turn
            matid_bound = self.matid[arg_bound_out]

            for mat_id in self.unique_matid:
                element_mass = np.sum(m_bound[matid_bound == mat_id])

                array_name = self.Di_id_mat[mat_id] + "_mass"
                element_mass_array[array_name][remnant_id - 1] = (
                    element_mass / Bound.M_earth
                )
                array_name = self.Di_id_mat[mat_id] + "_ratio"
                element_ratio_array[array_name][remnant_id - 1] = (
                    element_mass / rem_mass
                )

            bound_id[remnant_id - 1] = remnant_id
            m_rem[remnant_id - 1] = rem_mass / Bound.M_earth
            mass_ratio[remnant_id - 1] = rem_mass / Bound.M_earth / self.total_mass
            num_par_rem[remnant_id - 1] = np.sum(bound == remnant_id)

            remnant_id += 1
            bad_seeds = 0

            if remnant_id > self.num_rem:
                # print('Reach maximum number of remnants')
                # print('')
                break

        bound[bound == -1] = 0

        # reorder the bound mass to print out the largest one first
        arg_sel_desc = np.argsort(m_rem)[::-1]
        bound_id = bound_id[arg_sel_desc]
        mass_ratio = mass_ratio[arg_sel_desc]
        num_par_rem = num_par_rem[arg_sel_desc]

        for mat_id in self.unique_matid:
            array_name = self.Di_id_mat[mat_id] + "_mass"
            element_mass_array[array_name] = element_mass_array[array_name][
                arg_sel_desc
            ]
            array_name = self.Di_id_mat[mat_id] + "_ratio"
            element_ratio_array[array_name] = element_ratio_array[array_name][
                arg_sel_desc
            ]
        m_rem = m_rem[arg_sel_desc]

        cp_bid = deepcopy(bound)  # make a copy
        for i in range(len(bound_id)):
            if bound_id[i] != 0:
                cp_bid[bound == bound_id[i]] = i + 1
                bound_id[i] = i + 1

        bound = cp_bid
        self.bound_id = bound_id
        self.bound = bound
        self.mass_ratio = mass_ratio
        self.m_rem = m_rem
        self.num_par_rem = num_par_rem
        self.element_mass_array = element_mass_array
        self.element_ratio_array = element_ratio_array
        if self.verbose != 0:
            self.print_info()

    def calculate_mass(self):
        """recalculate the mass of the bound particles after re-disribution"""
        m_rem = np.zeros(self.num_rem)
        num_par_rem = np.zeros(self.num_rem, dtype=int)
        mass_ratio = np.zeros(self.num_rem)  # M_rem / M_total

        element_ratio_array = {}
        element_mass_array = {}

        for mat_id in self.unique_matid:
            array_name = self.Di_id_mat[mat_id] + "_ratio"
            element_ratio_array[array_name] = np.zeros(self.num_rem)
            array_name = self.Di_id_mat[mat_id] + "_mass"
            element_mass_array[array_name] = np.zeros(self.num_rem)

        for remnant_id in self.bound_id[:-1]:
            remnant_id = int(remnant_id)
            arg_bound_out = self.bound == remnant_id
            m_bound = self.m[arg_bound_out]
            rem_mass = np.sum(m_bound)  # mass of the remnant in this turn
            matid_bound = self.matid[arg_bound_out]

            for mat_id in self.unique_matid:
                element_mass = np.sum(m_bound[matid_bound == mat_id])

                array_name = self.Di_id_mat[mat_id] + "_mass"
                element_mass_array[array_name][remnant_id - 1] = (
                    element_mass / Bound.M_earth
                )
                array_name = self.Di_id_mat[mat_id] + "_ratio"
                element_ratio_array[array_name][remnant_id - 1] = (
                    element_mass / rem_mass
                )

            m_rem[remnant_id - 1] = rem_mass / Bound.M_earth
            mass_ratio[remnant_id - 1] = rem_mass / Bound.M_earth / self.total_mass
            num_par_rem[remnant_id - 1] = np.sum(self.bound == remnant_id)

        arg_sel_desc = np.argsort(m_rem)[::-1]
        bound_id = self.bound_id[arg_sel_desc]
        mass_ratio = mass_ratio[arg_sel_desc]
        num_par_rem = num_par_rem[arg_sel_desc]

        for mat_id in self.unique_matid:
            array_name = self.Di_id_mat[mat_id] + "_mass"
            element_mass_array[array_name] = element_mass_array[array_name][
                arg_sel_desc
            ]
            array_name = self.Di_id_mat[mat_id] + "_ratio"
            element_ratio_array[array_name] = element_ratio_array[array_name][
                arg_sel_desc
            ]
        m_rem = m_rem[arg_sel_desc]

        cp_bid = deepcopy(self.bound)  # make a copy
        for i in range(len(bound_id)):
            if bound_id[i] != 0:
                cp_bid[self.bound == bound_id[i]] = i + 1
                bound_id[i] = i + 1

        bound = cp_bid
        self.bound_id = bound_id
        self.bound = bound
        self.mass_ratio = mass_ratio
        self.m_rem = m_rem
        self.num_par_rem = num_par_rem
        self.element_mass_array = element_mass_array
        self.element_ratio_array = element_ratio_array
        if self.verbose != 0:
            self.print_info()

    def re_distribute(self, verbose=1):
        """After the bound remnants are found, recalcuate is any particles in a remnant should be redistributed to other remnants.

        Args:
            verbose (int, optional): Whether to print out information verbosely Defaults to 1.
        """
        # bound_cp = copy(self.bound)
        for bid in self.bound_id[::-1]:
            bound_cp = copy(self.bound)
            for rem_bid in self.bound_id[self.bound_id != bid]:
                rem_m = np.sum(self.m[self.bound == rem_bid])

                rem_com = (
                    np.sum(
                        self.pos[self.bound == rem_bid]
                        * self.m[self.bound == rem_bid, np.newaxis],
                        axis=0,
                    )
                    / rem_m
                )
                rem_vel = (
                    np.sum(
                        self.vel[self.bound == rem_bid]
                        * self.m[self.bound == rem_bid, np.newaxis],
                        axis=0,
                    )
                    / rem_m
                )

                c_rem_m = np.sum(self.m[self.bound == bid])

                c_rem_com = (
                    np.sum(
                        self.pos[self.bound == bid]
                        * self.m[self.bound == bid, np.newaxis],
                        axis=0,
                    )
                    / c_rem_m
                )
                c_rem_vel = (
                    np.sum(
                        self.vel[self.bound == bid]
                        * self.m[self.bound == bid, np.newaxis],
                        axis=0,
                    )
                    / c_rem_m
                )

                ke = (
                    0.5
                    * self.m[self.bound == bid]
                    * np.sum((self.vel[self.bound == bid] - rem_vel) ** 2, axis=1)
                )
                pe = (-Bound.G * rem_m * self.m[self.bound == bid]) / np.hypot(
                    self.pos[self.bound == bid, 2] - rem_com[2],
                    np.hypot(
                        self.pos[self.bound == bid, 0] - rem_com[0],
                        self.pos[self.bound == bid, 1] - rem_com[1],
                    ),
                )

                c_ke = (
                    0.5
                    * self.m[self.bound == bid]
                    * np.sum((self.vel[self.bound == bid] - c_rem_vel) ** 2, axis=1)
                )
                c_pe = (-Bound.G * c_rem_m * self.m[self.bound == bid]) / np.hypot(
                    self.pos[self.bound == bid, 2] - c_rem_com[2],
                    np.hypot(
                        self.pos[self.bound == bid, 0] - c_rem_com[0],
                        self.pos[self.bound == bid, 1] - c_rem_com[1],
                    ),
                )

                sel_redis_bound = (ke + pe < 0.0) & (c_pe > pe)

                bid_mask = self.bound == bid
                update_mask = np.zeros_like(bound_cp, dtype=bool)
                update_mask[bid_mask] = sel_redis_bound

                bound_cp[update_mask] = rem_bid

                if verbose:
                    print(
                        "Remnant %d: %d particles are redistributed to remnant %d"
                        % (bid, np.sum(sel_redis_bound), rem_bid)
                    )

            self.bound = bound_cp

        self.calculate_mass()

    def print_info(self):
        i = 0
        print("++++++++++++++++++++++++++++++++++++++++++++++")
        print(
            r"+ ALL BOUND MASS = %.2f %% total initial mass+"
            % (100 * np.sum(self.mass_ratio))
        )
        print("++++++++++++++++++++++++++++++++++++++++++++++")
        print("\n")
        print("------------------------------------------")
        while self.m_rem[i] > 0:
            print("Remnant %d mass = %.5f \u004D\u2295" % ((i + 1), self.m_rem[i]))
            print(r"Number of particles = %d" % self.num_par_rem[i])
            print(r"Mass ratio = %.2f %%" % (100 * self.mass_ratio[i]))
            for mat_id in self.unique_matid:
                array_name = self.Di_id_mat[mat_id] + "_ratio"
                print(
                    "%s ratio = %.2f %%"
                    % (
                        self.Di_id_mat[mat_id],
                        100 * self.element_ratio_array[array_name][i],
                    )
                )
            print("")
            print("------------------------------------------")
            print("")
            i += 1
            if i > self.num_rem - 1:
                break

    def source_track(self, verbose=1):
        # self.matid_tar_imp = self.matid
        # self.matid_tar_imp[self.npt <= self.pid] += Bound.id_body

        element_target_ratio_array = {}
        for mat_id in self.unique_matid:
            array_name = self.Di_id_mat[mat_id] + "_ratio_from_target"
            element_target_ratio_array[array_name] = np.zeros(
                np.count_nonzero(self.bound_id)
            )
        for rem_id in self.bound_id[self.bound_id != 0]:
            for mat_id in self.unique_matid:
                array_name = self.Di_id_mat[mat_id] + "_ratio_from_target"
                mass_array_name = self.Di_id_mat[mat_id] + "_mass"
                element_mass = self.element_mass_array[mass_array_name][int(rem_id) - 1]
                target_ratio = (
                    np.sum(
                        self.m[
                            np.logical_and(
                                self.bound == rem_id, self.matid_tar_imp == mat_id
                            )
                        ]
                    )
                    / Bound.M_earth
                    / element_mass
                )
                element_target_ratio_array[array_name][int(rem_id) - 1] = target_ratio

                if verbose:
                    print(
                        "-----------------------------------------------------------------"
                    )
                    print(
                        "In remnant %d, : ratio = %.2f %% %s is from target"
                        % (
                            int(rem_id),
                            100 * target_ratio,
                            self.Di_id_mat[mat_id],
                        )
                    )
                    print(
                        "-----------------------------------------------------------------"
                    )
        self.element_target_ratio_array = element_target_ratio_array

    def write_bound_id(self):
        f = h5py.File(self.filename, "r+")
        if "GasParticles/boundIDs" in f:
            del f["GasParticles/boundIDs"]
        if "GasParticles/npt" in f:
            del f["GasParticles/npt"]
        f["GasParticles"].create_dataset("boundIDs", data=self.bound, dtype="d")
        f["GasParticles"].create_dataset("npt", data=np.array([self.npt]), dtype="d")

    def all_in_one(
        self,
        loc_tar=None,
        loc_imp=None,
    ):
        assert loc_tar is not None
        assert loc_imp is not None

        Q_R, Q_R_norm, Q_RD_star_prime, M_tot, M_tar, M_imp = swtools.edacm(
            b=self.b, v=self.v * 1e5, loc_tar=loc_tar, loc_imp=loc_imp
        )
        self.m_tar = M_tar
        self.m_imp = M_imp
        self.total_mass = M_tot

        self.find_bound()

        accretion_rate = (self.m_rem[0] - self.m_tar) / (self.total_mass - self.m_tar)

        iron_key_list = np.array([100, 300, 401, 402])
        core_id = np.intersect1d(iron_key_list, self.unique_matid)
        core_array_name = self.Di_id_mat[int(core_id)] + "_ratio"
        self.gamma = (self.total_mass - self.m_tar) / self.m_tar
        d = {
            "M_tar": self.m_tar,
            "M_total": self.total_mass,
            "gamma": self.gamma,
            "npt": self.npt,
            "b": self.b,
            "v": self.v,
            "hit_dir": self.hit_direction,
            "attributes": self.attri,
            "m_lr": self.m_rem[0],
            "ratio_lr": self.mass_ratio[0],
            "targ_ratio_lr": (1 + self.gamma) * self.mass_ratio[0],
            "Z_Fe_lr": self.element_ratio_array[core_array_name][0],
            "np_lr": self.num_par_rem[0],
            "m_slr": self.m_rem[1:2],
            "ratio_slr": self.mass_ratio[1:2],
            "targ_ratio_slr": (1 + self.gamma) * self.mass_ratio[1:2],
            "Z_Fe_slr": self.element_ratio_array[core_array_name][1:2],
            "np_slr": self.num_par_rem[1:2],
            "accretion_rate": accretion_rate,
            "Q_R_norm": Q_R_norm,
            "Q_RD_star_prime": Q_RD_star_prime,
            "Q_R": Q_R,
        }

        return pd.Series(data=d)

    def calculate_entropy(self):
        if not hasattr(self, "entropy"):
            woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "ANEOS_Fe85Si15"])
            self.entropy = woma.eos.eos.A1_s_u_rho(self.u, self.rho_mks, self.matid)

    def calculate_T(self):
        if not hasattr(self, "T"):
            woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "ANEOS_Fe85Si15"])
            self.T = woma.eos.eos.A1_T_u_rho(self.u, self.rho_mks, self.matid)

    def total_vap_fraction(self, calculate_targ=False, verbose=1):
        """
        Calculates the total vapour fraction of core and mantle materials.
        """
        if not hasattr(self, "entropy"):
            woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "ANEOS_Fe85Si15"])
            self.entropy = woma.eos.eos.A1_s_u_rho(self.u, self.rho_mks, self.matid)

        self.particle_vapour_fraction = np.zeros(len(self.pid))

        core_id = np.intersect1d(self.iron_key_list, self.unique_matid)
        core_arg = self.matid == core_id
        mantle_id = np.intersect1d(self.si_key_list, self.unique_matid)
        mantle_arg = self.matid == mantle_id

        core_vf = VapourFrc(core_id, self.entropy[core_arg], self.p_mks[core_arg])
        core_vapour_fraction = core_vf.vapour_fraction()

        mantle_vf = VapourFrc(
            mantle_id, self.entropy[mantle_arg], self.p_mks[mantle_arg]
        )
        mantle_vapour_fraction = mantle_vf.vapour_fraction()

        self.core_vapour_fraction = np.sum(
            self.m[core_arg] * core_vapour_fraction
        ) / np.sum(self.m[core_arg])

        self.mantle_vapour_fraction = np.sum(
            self.m[mantle_arg] * mantle_vapour_fraction
        ) / np.sum(self.m[mantle_arg])

        self.particle_vapour_fraction[core_arg] = core_vapour_fraction
        self.particle_vapour_fraction[mantle_arg] = mantle_vapour_fraction

        sel_super_mantle = mantle_vf.super_critical()
        sel_super_core = core_vf.super_critical()

        self.mantle_super_critical = np.sum(
            self.m[mantle_arg][sel_super_mantle]
        ) / np.sum(self.m[mantle_arg])

        self.core_super_critical = np.sum(self.m[core_arg][sel_super_core]) / np.sum(
            self.m[core_arg]
        )

        self.particle_super_critical_fraction = np.zeros(len(self.pid))

        self.particle_super_critical_fraction[core_arg] = sel_super_core
        self.particle_super_critical_fraction[mantle_arg] = sel_super_mantle

        if calculate_targ:
            self.targ_core_vapour_fraction = np.sum(
                self.particle_vapour_fraction[self.matid_tar_imp == core_id]
                * self.m[self.matid_tar_imp == core_id]
            ) / np.sum(self.m[self.matid_tar_imp == core_id])
            self.targ_mantle_vapour_fraction = np.sum(
                self.particle_vapour_fraction[self.matid_tar_imp == mantle_id]
                * self.m[self.matid_tar_imp == mantle_id]
            ) / np.sum(self.m[self.matid_tar_imp == mantle_id])

            self.targ_core_sup_fraction = np.sum(
                self.particle_super_critical_fraction[self.matid_tar_imp == core_id]
                * self.m[self.matid_tar_imp == core_id]
            ) / np.sum(self.m[self.matid_tar_imp == core_id])
            self.targ_mantle_sup_fraction = np.sum(
                self.particle_super_critical_fraction[self.matid_tar_imp == mantle_id]
                * self.m[self.matid_tar_imp == mantle_id]
            ) / np.sum(self.m[self.matid_tar_imp == mantle_id])

        if verbose:
            print("{:.2f} % of core vapourized".format(100 * self.core_vapour_fraction))
            print(
                "{:.2f} % of mantle vapourized".format(
                    100 * self.mantle_vapour_fraction
                )
            )
            print(
                "-----------------------------------------------------------------------"
            )
            print(
                "{:.2f} % of core in super critical state".format(
                    100 * self.core_super_critical
                )
            )
            print(
                "{:.2f} % of mantle in super critical state".format(
                    100 * self.mantle_super_critical
                )
            )
            print(
                "-----------------------------------------------------------------------"
            )
            if calculate_targ:
                print(
                    "{:.2f} % of core in target is vapourized".format(
                        100 * self.targ_core_vapour_fraction
                    )
                )
                print(
                    "{:.2f} % of mantle in target is vapourized".format(
                        100 * self.targ_mantle_vapour_fraction
                    )
                )
                print(
                    "-----------------------------------------------------------------------"
                )
                print(
                    "{:.2f} % of core in target is super critical".format(
                        100 * self.targ_core_sup_fraction
                    )
                )
                print(
                    "{:.2f} % of mantle in target is super critical".format(
                        100 * self.targ_mantle_sup_fraction
                    )
                )

    def rem_vap_fraction(self, rem_id=1, verbose=1):
        """Given remnant id, calculated the vapour fraction in the remnant

        Args:
            rem_id (int, optional): Id of the remnant. Defaults to 1, which is the largest remnant.
            verbose (int, optional): Print out the results or not. Defaults to 1.
        """
        if not hasattr(self, "entropy"):
            woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "ANEOS_Fe85Si15"])
            self.entropy = woma.eos.eos.A1_s_u_rho(self.u, self.rho_mks, self.matid)

        if not hasattr(self, "rem_core_vapour_fraction"):
            self.rem_core_vapour_fraction = np.zeros(self.num_rem)
        if not hasattr(self, "rem_mantle_vapour_fraction"):
            self.rem_mantle_vapour_fraction = np.zeros(self.num_rem)

        core_id = np.intersect1d(self.iron_key_list, self.unique_matid)
        core_arg = np.logical_and(self.matid == core_id, self.bound == rem_id)

        mantle_id = np.intersect1d(self.si_key_list, self.unique_matid)
        mantle_arg = np.logical_and(self.matid == mantle_id, self.bound == rem_id)

        core_vf = VapourFrc(core_id, self.entropy[core_arg], self.p_mks[core_arg])
        core_vapour_fraction = core_vf.vapour_fraction()

        mantle_vf = VapourFrc(
            mantle_id, self.entropy[mantle_arg], self.p_mks[mantle_arg]
        )
        mantle_vapour_fraction = mantle_vf.vapour_fraction()

        self.rem_core_vapour_fraction[rem_id - 1] = np.sum(
            self.m[core_arg] * core_vapour_fraction
        ) / np.sum(self.m[core_arg])
        self.rem_mantle_vapour_fraction[rem_id - 1] = np.sum(
            self.m[mantle_arg] * mantle_vapour_fraction
        ) / np.sum(self.m[mantle_arg])

        if verbose:
            print(
                "In the remnant {:d}: {:.2f} % of iron vapourized".format(
                    rem_id, 100 * self.rem_core_vapour_fraction[rem_id - 1]
                )
            )
            print(
                "In the remnant {:d}: {:.2f} % of si vapourized".format(
                    rem_id, 100 * self.rem_mantle_vapour_fraction[rem_id - 1]
                )
            )

    def disk_planet_mass(self, Navg=200, verbose=1):
        """seperate the largest bound remnant to disk and planet particles

        Args:
            Navg (int): number of particles to average over
            verbose (int, optional): _description_. Defaults to 1.
        """
        # Need first run find_bound() to get the bound particles
        if not hasattr(self, "bound"):
            raise AttributeError("Need to run find_bound() first")
        # recenter the position and velocity to the cm,cvel of the largest remnant
        pos_bnd = self.pos[self.bound == 1]
        vel_bnd = self.vel[self.bound == 1]
        m_bnd = self.m[self.bound == 1]
        pid_bnd = self.pid[self.bound == 1]
        rho_bnd = self.rho_mks[self.bound == 1]

        pos_bnd_centerM = np.sum(pos_bnd * m_bnd[:, np.newaxis], axis=0) / np.sum(m_bnd)
        pos_bnd -= pos_bnd_centerM
        vel_bnd_centerM = np.sum(vel_bnd * m_bnd[:, np.newaxis], axis=0) / np.sum(m_bnd)
        vel_bnd -= vel_bnd_centerM
        # calcualte the radius of the particles
        r_bnd = np.sqrt(np.sum(pos_bnd**2, axis=1))
        # calculate the kenetic energy of the particles
        ke_bnd = 0.5 * m_bnd * np.sum(vel_bnd**2, axis=1)
        # Using a moving average calculation to find the planet-disk boundary (max particle KE)
        r_bnd_argsort = np.argsort(r_bnd)
        Navg = Navg
        KEmax = 0
        rKEmax = 0
        for i in range(len(r_bnd) - Navg):
            # calculate the average kinetic energy of the particles
            KEavg = np.sum(ke_bnd[r_bnd_argsort[i : i + Navg]]) / Navg
            if KEavg > KEmax:
                KEmax = KEavg
                rKEmax = r_bnd[r_bnd_argsort[i]]
        # particles with r <= rKEmax are planet particles
        self.planet_sel = r_bnd <= rKEmax
        # particles with r > rKEmax are disk particles
        self.disk_sel = r_bnd > rKEmax

        self.planet_pid = pid_bnd[self.planet_sel]
        self.disk_pid = pid_bnd[self.disk_sel]

        self.planet_m = np.sum(m_bnd[self.planet_sel]) / Bound.M_earth  # in earth mass
        self.disk_m = np.sum(m_bnd[self.disk_sel]) / Bound.M_earth  # in earth mass

        self.planet_vel = vel_bnd[self.planet_sel]
        self.disk_vel = vel_bnd[self.disk_sel]

        # calculate the angular momentum of the planet particles along the z axis
        self.planet_Lz = np.sum(
            np.cross(pos_bnd[self.planet_sel], vel_bnd[self.planet_sel])[:, 2]
        )
        # calculate the angular momentum of the disk particles along the z axis
        self.disk_Lz = np.sum(
            np.cross(pos_bnd[self.disk_sel], vel_bnd[self.disk_sel])[:, 2]
        )
        # calculate the total angular momentum of the all the bound particles along the z axis
        self.total_bnd_Lz = np.sum(np.cross(pos_bnd, vel_bnd)[:, 2])
        # calculate the total angular momentum of the all the particles along the z axis
        self.total_Lz = np.sum(np.cross(self.pos, self.vel)[:, 2])

        self.disk_rho_mks = rho_bnd[self.disk_sel]
        self.planet_rho_mks = rho_bnd[self.planet_sel]

    def calculate_disk_vapor_fraction(self, verbose=1):
        """Calculate the vapour fraction of the disk particles

        Args:
            verbose (int, optional): _description_. Defaults to 1.
        """
        if not hasattr(self, "disk_pid"):
            raise AttributeError("Need to run disk_planet_mass() first")

        if not hasattr(self, "entropy"):
            woma.load_eos_tables(["ANEOS_iron", "ANEOS_forsterite", "ANEOS_Fe85Si15"])
            self.entropy = woma.eos.eos.A1_s_u_rho(self.u, self.rho_mks, self.matid)

        if not hasattr(self, "disk_core_vapour_fraction"):
            self.disk_core_vapour_fraction = np.zeros(self.num_rem)
        if not hasattr(self, "disk_mantle_vapour_fraction"):
            self.disk_mantle_vapour_fraction = np.zeros(self.num_rem)

        core_id = np.intersect1d(self.iron_key_list, self.unique_matid)
        core_arg = np.logical_and(
            self.matid == core_id, np.isin(self.pid, self.disk_pid)
        )

        mantle_id = np.intersect1d(self.si_key_list, self.unique_matid)
        mantle_arg = np.logical_and(
            self.matid == mantle_id, np.isin(self.pid, self.disk_pid)
        )

        core_vf = VapourFrc(core_id, self.entropy[core_arg], self.p_mks[core_arg])
        core_vapour_fraction = core_vf.vapour_fraction()

        mantle_vf = VapourFrc(
            mantle_id, self.entropy[mantle_arg], self.p_mks[mantle_arg]
        )
        mantle_vapour_fraction = mantle_vf.vapour_fraction()

        self.disk_core_vapour_fraction[0] = np.sum(
            self.m[core_arg] * core_vapour_fraction
        ) / np.sum(self.m[core_arg])
        self.disk_mantle_vapour_fraction[0] = np.sum(
            self.m[mantle_arg] * mantle_vapour_fraction
        ) / np.sum(self.m[mantle_arg])

        self.disk_vapor_fraction = (
            np.sum(self.m[core_arg] * core_vapour_fraction)
            + np.sum(self.m[mantle_arg] * mantle_vapour_fraction)
        ) / (np.sum(self.m[core_arg]) + np.sum(self.m[mantle_arg]))

    def rotate_along_z(self, rot_bound=True, percentile=5, verbose=1):
        """Rotate coordinate of particles along the z axis.
        Roration angle is calcualte base on percentile% and percentile% of x and y coordinate
        rot_bound: If false, rotate based on the particles coordinate of around all particles
                   if 1, then base on the largest remnant particles
                   if 2, then base on all the bound remnant particles
        percentile:how many particles to select to calculate the rotation angle
        """

        if rot_bound == 1:
            rem_pos = self.pos[self.bound == 1]
            x_p_low = np.percentile(np.sort(rem_pos[:, 0]), percentile)
            x_p_high = np.percentile(np.sort(rem_pos[:, 0]), 100 - percentile)

            y_p_low = np.percentile(np.sort(rem_pos[:, 1]), percentile)
            y_p_high = np.percentile(np.sort(rem_pos[:, 1]), 100 - percentile)

        elif rot_bound == 2:
            rem_pos = self.pos[self.bound > 0]
            x_p_low = np.percentile(np.sort(rem_pos[:, 0]), percentile)
            x_p_high = np.percentile(np.sort(rem_pos[:, 0]), 100 - percentile)

            y_p_low = np.percentile(np.sort(rem_pos[:, 1]), percentile)
            y_p_high = np.percentile(np.sort(rem_pos[:, 1]), 100 - percentile)

        elif rot_bound == 0:
            x_p_low = np.percentile(np.sort(self.pos[:, 0]), percentile)
            x_p_high = np.percentile(np.sort(self.pos[:, 0]), 100 - percentile)

            y_p_low = np.percentile(np.sort(self.pos[:, 1]), percentile)
            y_p_high = np.percentile(np.sort(self.pos[:, 1]), 100 - percentile)
        else:
            raise ValueError("rot_bound should be 0,1,2")

        y_dis = abs(y_p_low - y_p_high)
        x_dis = abs(x_p_low - x_p_high)

        radians = np.arctan(x_dis / y_dis)
        if verbose:
            print("radians : %.4f" % radians)
            print("Degree : %.4f" % np.rad2deg(radians))

        rotation_matrix = np.array(
            [
                [np.cos(radians), -np.sin(radians), 0],
                [np.sin(radians), np.cos(radians), 0],
                [0, 0, 1],
            ]
        )

        self.pos = np.dot(self.pos, rotation_matrix)

    def basic_plot(
        self,
        mode=0,
        matid_plot=-1,
        extent=None,
        equal_axis=False,
        color_mode0=False,
        sel_pid=None,
        selp_size=3,
        selp_color="cyan",
    ):
        """
        Plot the bound particles.
        mode = 0, plot all remnants with different colors
        mode = -2, only plot the second largest remnant
        mode = -1, only plot the largest remnant
        mode = 1 show largest remnant only
        mode = 2 show second largest remnant only
        mode = 3 show third largest remnant only
        ...
        mode =10 show 10th largest remnant only
        matid_plot: if you want to plot only one material, set this to the material id. -1 means all materials.
        extent: set the extent of the plot, [xmin, xmax, ymin, ymax,zmin,zmax], unit in Rearth radius.
        equal_axis: set the axis to be equal or not.
        color_mode0: if you want the particles to be colored by their material id, set this to True.
        sel_pid: if you want to annotate additionally a subset of particles, set this to the particle id array.
        selp_size: size of the annotated particles.
        selp_color: color of the annotated particles.
        """
        colours = np.empty(len(self.pid), dtype=object)
        sizes = np.zeros(len(self.pid))

        for matid in np.unique(self.matid_tar_imp):
            colours[self.matid_tar_imp == matid] = self.Di_id_colour[matid]
            sizes[self.matid_tar_imp == matid] = self.Di_id_size[matid]

        if sel_pid is not None:
            colours[np.in1d(self.pid, sel_pid)] = selp_color
            sizes[np.in1d(self.pid, sel_pid)] = (
                selp_size * sizes[np.in1d(self.pid, sel_pid)]
            )

        fig = plt.figure(figsize=(12, 6))
        ax1, ax2 = fig.subplots(1, 2)

        # if ax1 is None:
        #     output_fig = False
        # else:
        #     ax1 = ax1_in
        #     ax2 = ax2_in
        #     output_fig = True

        if mode < 0:
            if matid_plot == -1:
                sel_lr_rem_arg = self.bound == -mode
            else:
                sel_lr_rem_arg = np.logical_and(
                    self.bound == -mode, self.matid == matid_plot
                )
            # recnter the position to the cm of the largest remnant
            lr_center = np.sum(
                self.pos[sel_lr_rem_arg] * self.m[sel_lr_rem_arg, np.newaxis], axis=0
            ) / np.sum(self.m[sel_lr_rem_arg])
            pos_lr = self.pos[sel_lr_rem_arg] - lr_center

            colours_lr = colours[sel_lr_rem_arg]
            sizes_lr = sizes[sel_lr_rem_arg]

            # plot particles with z or x < 0.1 R_earth and sort them by z or x
            arg_sort_pos_z = np.argsort(pos_lr[:, 2])
            arg_sort_pos_x = np.argsort(pos_lr[:, 0])

            search_sort_z = np.searchsorted(
                pos_lr[arg_sort_pos_z, 2],
                pos_lr[:, 2][pos_lr[:, 2] <= 0.1 * Bound.R_earth],
            )
            search_sort_x = np.searchsorted(
                pos_lr[arg_sort_pos_x, 0],
                pos_lr[:, 0][pos_lr[:, 0] <= 0.1 * Bound.R_earth],
            )

            arg_z = arg_sort_pos_z[search_sort_z]
            arg_x = arg_sort_pos_x[search_sort_x]

            ax1.scatter(
                pos_lr[arg_z, 0] / Bound.R_earth,
                pos_lr[arg_z, 1] / Bound.R_earth,
                c=colours_lr[arg_z],
                s=sizes_lr[arg_z],
            )

            ax2.scatter(
                pos_lr[arg_x, 1] / Bound.R_earth,
                pos_lr[arg_x, 2] / Bound.R_earth,
                c=colours_lr[arg_x],
                s=sizes_lr[arg_x],
            )

        elif mode == 0:
            # if number of remnants is less than 9, then each element's color will be picked here
            # recenterization
            self.generate_rem_colour()

            pos_center = np.sum(self.pos * self.m[:, np.newaxis], axis=0) / np.sum(
                self.m
            )
            self.pos -= pos_center
            if not color_mode0:
                for bnd_id in self.bound_id:
                    if bnd_id != 0:
                        colours[self.bound == bnd_id] = self.rem_colours[
                            int(bnd_id) - 1
                        ]

            rem_labels = ["remnant {:d}".format(int(i)) for i in self.bound_id]
            if matid_plot == -1:
                arg_bound = self.bound != 0
            else:
                arg_bound = np.logical_and(self.bound != 0, self.matid == matid_plot)

            ax1.scatter(
                self.pos[arg_bound, 0] / Bound.R_earth,
                self.pos[arg_bound, 1] / Bound.R_earth,
                s=sizes[arg_bound],
                c=colours[arg_bound],
                # label=rem_labels,
            )

            ax2.scatter(
                self.pos[arg_bound, 1] / Bound.R_earth,
                self.pos[arg_bound, 2] / Bound.R_earth,
                s=sizes[arg_bound],
                c=colours[arg_bound],
                # label=rem_labels,
            )
            # ax1.legend()
            # ax2.legend()

        elif mode > 0:
            self.generate_rem_colour()

            pos_center = np.sum(self.pos * self.m[:, np.newaxis], axis=0) / np.sum(
                self.m
            )
            self.pos -= pos_center

            colours[self.bound == mode] = self.rem_colours[mode - 1]
            sizes[self.bound == mode] = (
                3 * sizes[self.bound == mode]
            )  # make the remnant larger and clearer to see

            arg_sort_pos_z = np.argsort(self.pos[:, 2])
            arg_sort_pos_x = np.argsort(self.pos[:, 0])

            search_sort_z = np.searchsorted(
                self.pos[arg_sort_pos_z, 2],
                self.pos[:, 2][self.pos[:, 2] <= 0.1 * Bound.R_earth],
            )
            search_sort_x = np.searchsorted(
                self.pos[arg_sort_pos_x, 0],
                self.pos[:, 0][self.pos[:, 0] <= 0.1 * Bound.R_earth],
            )

            arg_z = arg_sort_pos_z[search_sort_z]
            arg_x = arg_sort_pos_x[search_sort_x]

            ax1.scatter(
                self.pos[arg_z, 0] / Bound.R_earth,
                self.pos[arg_z, 1] / Bound.R_earth,
                s=sizes[arg_z],
                c=colours[arg_z],
            )

            ax2.scatter(
                self.pos[arg_x, 1] / Bound.R_earth,
                self.pos[arg_x, 2] / Bound.R_earth,
                s=sizes[arg_x],
                c=colours[arg_x],
            )
        else:
            raise ValueError("mode must be an integer.")

        if equal_axis:
            ax1.set_aspect("equal", anchor="C")
            ax2.set_aspect("equal", anchor="C")

        if extent is not None:
            ax1.set_xlim(extent[0], extent[1])
            ax1.set_ylim(extent[2], extent[3])
            ax2.set_xlim(extent[2], extent[3])
            ax2.set_ylim(extent[4], extent[5])

        ax1.set_xlabel(r"x Position ($R_\oplus$)", fontsize=16)
        ax1.set_ylabel(r"y Position ($R_\oplus$)", fontsize=16)
        ax1.set_facecolor("#111111")
        ax2.set_xlabel(r"y Position ($R_\oplus$)", fontsize=16)
        ax2.set_ylabel(r"z Position ($R_\oplus$)", fontsize=16)
        ax2.set_facecolor("#111111")
        fig.tight_layout()

        # if output_fig:
        #     return fig, ax1, ax2
        # else:
        plt.show()
        plt.cla()
        plt.clf()
        plt.close()

    def generate_rem_colour(self):
        assert (
            self.bound_id is not None
        ), "bound_id is not generated yet. Please run find_bound() first."

        default_colours_rem_array = [
            "lime",
            "cyan",
            "yellow",
            "blue",
            "red",
            "green",
            "deeppink",
            "olive",
        ]
        # if number of remnants exceed 8, then random generate some colours
        if np.count_nonzero(self.bound_id) > len(default_colours_rem_array):
            rem_colours = []
            for i in range(np.count_nonzero(self.bound_id)):
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)
                hex_code = "#{:02x}{:02x}{:02x}".format(r, g, b)
                rem_colours.append(hex_code)
        else:
            rem_colours = default_colours_rem_array[: np.count_nonzero(self.bound_id)]
        self.rem_colours = rem_colours


class Snap:
    G = 6.67408e-11  # m^3 kg^-1 s^-2
    M_earth = 5.97240e24  # kg
    R_earth = 6.371e6  # m

    def __init__(self, filename, npt=1e9):
        self.filename = filename
        self.npt = npt
        self.load_data()
        self.init_scatter_colour()
        self.init_scatter_size()
        self.update_material_dictionary()

    def load_data(self):
        # load all the data from snapshot
        data = sw.load(self.filename)

        self.snaptime = data.metadata.time.value
        self.boxsize = data.metadata.boxsize[0]
        box_mid = 0.5 * self.boxsize.to(unyt.m)
        data.gas.coordinates.convert_to_mks()
        pos = data.gas.coordinates - box_mid
        self.pos = np.array(pos)
        data.gas.densities.convert_to_cgs()
        self.rho_cgs = np.array(data.gas.densities)
        data.gas.densities.convert_to_mks()
        self.rho_mks = np.array(data.gas.densities)
        data.gas.internal_energies.convert_to_mks()
        self.u = np.array(data.gas.internal_energies)
        data.gas.pressures.convert_to_mks()
        self.p_mks = np.array(data.gas.pressures)
        data.gas.potentials.convert_to_mks()
        self.pot = data.gas.potentials
        self.matid = np.array(data.gas.material_ids)
        self.pid = np.array(data.gas.particle_ids)
        data.gas.masses.convert_to_mks()
        self.m = np.array(data.gas.masses)
        data.gas.velocities.convert_to_mks()
        self.vel = np.array(data.gas.velocities)
        data.gas.smoothing_lengths.convert_to_mks()
        self.h = np.array(data.gas.smoothing_lengths)
        # set different id for target and impactor materials
        self.unique_matid = np.unique(self.matid)

        self.matid_tar_imp = deepcopy(self.matid)
        self.matid_tar_imp[self.npt <= self.pid] += Bound.id_body

    # plot the density and pressure distribution
    def dis_rho_p(self, sel_matid=None, output_fig=False):
        if sel_matid is None:
            sel = np.ones(len(self.pid), dtype=bool)
        else:
            sel = self.matid == sel_matid

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        ax1 = sns.kdeplot(self.rho_cgs[sel], ax=ax1)
        ax2 = sns.kdeplot(np.log10(self.p_mks[sel]), ax=ax2)
        ax1.set(xlabel="Density (g/cm^3)", ylabel="distribution")
        ax2.set(xlabel="Pressure (log10(Gpa))", ylabel="ditribution")

    def recenter(self, new_center):
        """recenter the position to new_center."""
        self.pos -= new_center

    def pos_com_center(self):
        """center the coordinates of particles by the center of mass."""
        pos_centerM = np.sum(self.pos * self.m[:, np.newaxis], axis=0) / np.sum(self.m)
        self.pos -= pos_centerM

        return 0

    def vel_com_center(self):
        """center the velocities of particles by the center of mass."""
        vel_centerM = np.sum(self.vel * self.m[:, np.newaxis], axis=0) / np.sum(self.m)
        self.vel -= vel_centerM

        return 0

    def dis_to_com(self, if_R_atmos=False):
        """Distance to the center of mass, normall used to calculate the radius
        Set "if_R_atmos" to True if calculate the radius including atmospher layer.
        """
        self.pos_com_center()

        atmos_key_list = np.array([0, 1, 2, 200, 305, 306, 307])
        uniq_mat = np.unique(self.matid)
        atmos_id = np.intersect1d(atmos_key_list, uniq_mat)

        if not if_R_atmos:
            pos_to_use = self.pos[self.matid != atmos_id]
        else:
            pos_to_use = self.pos
        pos_to_use = pos_to_use.squeeze()
        xy = np.hypot(pos_to_use[:, 0], pos_to_use[:, 1])
        self.dis_com_noorder = np.hypot(xy, pos_to_use[:, 2])
        self.dis_com_order = np.sort(self.dis_com_noorder)
        self.dis_com_outer = np.mean(self.dis_com_order[-200:])

        return 0

    def v_rms(self, if_R_atmos=False, verbose=1):
        """Calculate the root mean square velocity in m/s (mks unit). Normally useful when checking cooling snapshot."""
        self.dis_to_com(if_R_atmos=if_R_atmos)
        v_rms = np.sqrt(np.sum(self.vel**2) / len(self.m))
        v_esc = np.sqrt(2 * Snap.G * np.sum(self.m) / self.dis_com_outer)
        if verbose:
            print("Root mean square velocity is : % .4f km/s" % (v_rms / 1e3))
            print("Escape velocity is : % .4f km/s" % (v_esc / 1e3))
            print("V_rms/V_esc = %.4f %%" % (100 * (v_rms / v_esc)))
        return 0

    def splot(
        self,
        aspect="xy",
        extent=None,
        ax=None,
        equal_axis=True,
        sel_pid=None,
        sel_matid=-1,
        selp_size=3,
        selp_color="cyan",
        figsize=(16, 9),
    ):
        """_summary_

        Args:
            aspect (bool, optional): _description_. Defaults to True.
            extent (list, optional): [xmin,xmax,ymin,ymax,zmin,zmax] Defaults to None.
            output_fig (bool, optional): _description_. Defaults to False.
            equal_axis (bool, optional): _description_. Defaults to True.
            sel_pid (_type_, optional): _description_. Defaults to None.
            sel_matid (int, optional): _description_. Defaults to -1.
            selp_size (int, optional): _description_. Defaults to 3.
            selp_color (str, optional): _description_. Defaults to "cyan".

        Returns:
           plt.axes:return particles plot
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            created_ax = True
        else:
            created_ax = False
            fig = ax.figure

        colours = np.empty(len(self.pid), dtype=object)
        sizes = np.zeros(len(self.pid))

        for matid in np.unique(self.matid_tar_imp):
            colours[self.matid_tar_imp == matid] = self.Di_id_colour[matid]
            sizes[self.matid_tar_imp == matid] = self.Di_id_size[matid]

        sel_pos = np.ones_like(self.pid, dtype=bool)

        if sel_pid is not None:
            colours[np.isin(self.pid, sel_pid)] = selp_color
            sizes[np.isin(self.pid, sel_pid)] = (
                selp_size * sizes[np.in1d(self.pid, sel_pid)]
            )
        if sel_matid >= 0:
            sel_pos[(self.matid != sel_matid)] = 0
            # colours = colours[sel_pos]
            # sizes = sizes[sel_pos]

        if extent is not None:
            extent = np.array(extent, dtype=float)
            extent *= Bound.R_earth

            sel_pos[
                (self.pos[:, 0] < extent[0])
                | (self.pos[:, 0] > extent[1])
                | (self.pos[:, 1] < extent[2])
                | (self.pos[:, 1] > extent[3])
                | (self.pos[:, 2] < extent[4])
                | (self.pos[:, 2] > extent[5])
            ] = 0

        plot_pos = self.pos[sel_pos]
        plot_colours = colours[sel_pos]
        plot_sizes = sizes[sel_pos]
        self.sel_pos = sel_pos
        # fig = plt.figure(figsize=(8, 8))
        # ax = fig.add_subplot(111)
        if aspect == "xy":
            arg_sort_pos_z = np.argsort(plot_pos[:, 2])
            search_sort_z = np.searchsorted(
                plot_pos[arg_sort_pos_z, 2],
                plot_pos[:, 2][plot_pos[:, 2] <= 0.1 * Bound.R_earth],
            )
            arg_z = arg_sort_pos_z[search_sort_z]
            ax.scatter(
                plot_pos[arg_z, 0] / Bound.R_earth,
                plot_pos[arg_z, 1] / Bound.R_earth,
                s=plot_sizes[arg_z],
                c=plot_colours[arg_z],
            )

            ax.set_xlabel(r"x Position ($R_\oplus$)", fontsize=16)
            ax.set_ylabel(r"y Position ($R_\oplus$)", fontsize=16)

        elif aspect == "yz":
            arg_sort_pos_x = np.argsort(plot_pos[:, 0])
            search_sort_x = np.searchsorted(
                plot_pos[arg_sort_pos_x, 0],
                plot_pos[:, 0][plot_pos[:, 0] <= 0.1 * Bound.R_earth],
            )
            arg_x = arg_sort_pos_x[search_sort_x]
            ax.scatter(
                plot_pos[arg_x, 1] / Bound.R_earth,
                plot_pos[arg_x, 2] / Bound.R_earth,
                s=plot_sizes[arg_x],
                c=plot_colours[arg_x],
            )
            ax.set_xlabel(r"y Position ($R_\oplus$)", fontsize=16)
            ax.set_ylabel(r"z Position ($R_\oplus$)", fontsize=16)
        else:
            arg_sort_pos_x = np.argsort(plot_pos[:, 1])
            search_sort_x = np.searchsorted(
                plot_pos[arg_sort_pos_x, 1],
                plot_pos[:, 1][plot_pos[:, 1] <= 0.1 * Bound.R_earth],
            )
            arg_x = arg_sort_pos_x[search_sort_x]
            ax.scatter(
                plot_pos[arg_x, 0] / Bound.R_earth,
                plot_pos[arg_x, 2] / Bound.R_earth,
                s=plot_sizes[arg_x],
                c=plot_colours[arg_x],
            )
            ax.set_xlabel(r"x Position ($R_\oplus$)", fontsize=16)
            ax.set_ylabel(r"z Position ($R_\oplus$)", fontsize=16)

        if equal_axis:
            ax.set_aspect("equal", anchor="C")

        ax.set_facecolor("#111111")
        ax.tick_params(axis="both", which="major", labelsize=14)

        return ax

    def update_material_dictionary(self):
        type_factor = 100
        Di_mat_type = {
            "idg": 0,
            "Til": 1,
            "HM80": 2,
            "SESAME": 3,
            "ANEOS": 4,
        }
        Di_mat_id = {
            # Ideal Gas
            "idg_HHe": Di_mat_type["idg"] * type_factor,
            "idg_N2": Di_mat_type["idg"] * type_factor + 1,
            "idg_CO2": Di_mat_type["idg"] * type_factor + 2,
            # Tillotson
            "Til_iron": Di_mat_type["Til"] * type_factor,
            "Til_granite": Di_mat_type["Til"] * type_factor + 1,
            "Til_water": Di_mat_type["Til"] * type_factor + 2,
            "Til_basalt": Di_mat_type["Til"] * type_factor + 3,
            # Hubbard & MacFarlane (1980) Uranus/Neptune
            "HM80_HHe": Di_mat_type["HM80"] * type_factor,  # Hydrogen-helium atmosphere
            "HM80_ice": Di_mat_type["HM80"] * type_factor + 1,  # H20-CH4-NH3 ice mix
            "HM80_rock": Di_mat_type["HM80"] * type_factor
            + 2,  # SiO2-MgO-FeS-FeO rock mix
            # SESAME etc
            "SESAME_iron": Di_mat_type["SESAME"] * type_factor,  # 2140
            "SESAME_basalt": Di_mat_type["SESAME"] * type_factor + 1,  # 7530
            "SESAME_water": Di_mat_type["SESAME"] * type_factor + 2,  # 7154
            "SS08_water": Di_mat_type["SESAME"] * type_factor
            + 3,  # Senft & Stewart (2008)
            "AQUA": Di_mat_type["SESAME"] * type_factor + 4,  # Haldemann+2020
            "CMS19_H": Di_mat_type["SESAME"] * type_factor
            + 5,  # Chabrier+2019 Hydrogen
            "CMS19_He": Di_mat_type["SESAME"] * type_factor + 6,  # Helium
            "CD21_HHe": Di_mat_type["SESAME"] * type_factor + 7,  # H/He mixture Y=0.275
            # ANEOS
            "ANEOS_forsterite": Di_mat_type["ANEOS"]
            * type_factor,  # Stewart et al. (2019)
            "ANEOS_iron": Di_mat_type["ANEOS"] * type_factor + 1,  # Stewart (2020)
            "ANEOS_Fe85Si15": Di_mat_type["ANEOS"] * type_factor + 2,  # Stewart (2020)
        }
        Di_mat_id.update(
            {matname + "_2": mid + Bound.id_body for matname, mid in Di_mat_id.items()}
        )

        # Invert so the ID are the keys
        self.Di_id_mat = {mat_id: mat for mat, mat_id in Di_mat_id.items()}

        atmos_key_list = np.array([0, 1, 2, 200, 305, 306, 307])
        self.atmos_key_list = np.concatenate(
            (atmos_key_list, atmos_key_list + Bound.id_body)
        )
        water_key_list = np.array([102, 201, 302, 303, 304])
        self.water_key_list = np.concatenate(
            (water_key_list, water_key_list + Bound.id_body)
        )
        iron_key_list = np.array([100, 300, 401, 402])
        self.iron_key_list = np.concatenate(
            (iron_key_list, iron_key_list + Bound.id_body)
        )
        si_key_list = np.array([101, 103, 202, 301, 400])
        self.si_key_list = np.concatenate((si_key_list, si_key_list + Bound.id_body))

        Di_id_colour = {}
        Di_id_size = {}

        for key in self.Di_id_mat.keys():
            if key in self.atmos_key_list:
                if key < Bound.id_body:
                    Di_id_colour[key] = self.colour_atmos_tar
                    Di_id_size[key] = self.size_atmos_tar
                else:
                    Di_id_colour[key] = self.colour_atmos_imp
                    Di_id_size[key] = self.size_atmos_imp

            elif key in self.water_key_list:
                if key < Bound.id_body:
                    Di_id_colour[key] = self.colour_water_tar
                    Di_id_size[key] = self.size_water_tar
                else:
                    Di_id_colour[key] = self.colour_water_imp
                    Di_id_size[key] = self.size_water_imp

            elif key in self.iron_key_list:
                if key < Bound.id_body:
                    Di_id_colour[key] = self.colour_iron_tar
                    Di_id_size[key] = self.size_iron_tar
                else:
                    Di_id_colour[key] = self.colour_iron_imp
                    Di_id_size[key] = self.size_iron_imp

            elif key in self.si_key_list:
                if key < Bound.id_body:
                    Di_id_colour[key] = self.colour_si_tar
                    Di_id_size[key] = self.size_si_tar
                else:
                    Di_id_colour[key] = self.colour_si_imp
                    Di_id_size[key] = self.size_si_imp

        self.Di_id_colour = Di_id_colour
        self.Di_id_size = Di_id_size

        return

    def set_scatter_colour(
        self,
        colour_iron_tar="tomato",
        colour_si_tar="mediumseagreen",
        colour_water_tar="skyblue",
        colour_atmos_tar="aliceblue",
        colour_iron_imp="sandybrown",
        colour_si_imp="pink",
        colour_water_imp="skyblue",
        colour_atmos_imp="aliceblue",
    ):
        self.colour_iron_tar = colour_iron_tar
        self.colour_si_tar = colour_si_tar
        self.colour_water_tar = colour_water_tar
        self.colour_atmos_tar = colour_atmos_tar

        self.colour_iron_imp = colour_iron_imp
        self.colour_si_imp = colour_si_imp
        self.colour_water_imp = colour_water_imp
        self.colour_atmos_imp = colour_atmos_imp
        self.update_material_dictionary()

    def set_scatter_size(
        self, size_iron=0.1, size_si=0.1, size_water=0.1, size_atmos=0.1
    ):
        self.size_iron_tar = size_iron
        self.size_si_tar = size_si
        self.size_water_tar = size_water
        self.size_atmos_tar = size_atmos

        self.size_iron_imp = size_iron
        self.size_si_imp = size_si
        self.size_water_imp = size_water
        self.size_atmos_imp = size_atmos

        self.update_material_dictionary()

    def init_scatter_colour(
        self,
        colour_iron_tar="tomato",
        colour_si_tar="mediumseagreen",
        colour_water_tar="skyblue",
        colour_atmos_tar="aliceblue",
        colour_iron_imp="sandybrown",
        colour_si_imp="pink",
        colour_water_imp="skyblue",
        colour_atmos_imp="aliceblue",
    ):
        """
        Initialise the scatter plot colour for different materials.
        """
        self.colour_iron_tar = colour_iron_tar
        self.colour_si_tar = colour_si_tar
        self.colour_water_tar = colour_water_tar
        self.colour_atmos_tar = colour_atmos_tar

        self.colour_iron_imp = colour_iron_imp
        self.colour_si_imp = colour_si_imp
        self.colour_water_imp = colour_water_imp
        self.colour_atmos_imp = colour_atmos_imp

    def init_scatter_size(
        self, size_iron=0.1, size_si=0.1, size_water=0.1, size_atmos=0.1
    ):
        """
        Initialise the scatter plot size for different materials.

        """
        self.size_iron_tar = size_iron
        self.size_si_tar = size_si
        self.size_water_tar = size_water
        self.size_atmos_tar = size_atmos

        self.size_iron_imp = size_iron
        self.size_si_imp = size_si
        self.size_water_imp = size_water
        self.size_atmos_imp = size_atmos


def main():
    import cProfile
    import pstats

    loc = "/Users/qb20321/Desktop/SWIFTother/test_snap/snapOUT_PLANETimpact_0d0h_1d58085_npt173750_3d16170_v36d8387kms_b0d000_pX_EiEf.hdf5"
    x = Bound(filename=loc, verbose=1, num_rem=6)

    with cProfile.Profile() as pr:
        x.find_bound()

    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.dump_stats(
        filename="/Users/qb20321/Desktop/SWIFTother/test_snap/boundmass_profiling.prog"
    )


if __name__ == "__main__":
    main()
