"""Simple teaching library for 3D Lennard-Jones molecular dynamics and analysis."""

from __future__ import annotations
import numpy as np

def make_fcc_positions(n_cells, box_length):
    """
    Build an FCC lattice with 4 atoms per unit cell.

    Parameters
    ----------
    n_cells : int
        Number of FCC unit cells in each dimension.
    box_length : float
        Simulation box length.

    Returns
    -------
    positions : (N, 3) ndarray
        Particle positions.
    """
    basis = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ])
    a = box_length / n_cells
    positions = []
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                cell_origin = np.array([i, j, k], dtype=float) * a
                for b in basis:
                    positions.append(cell_origin + a * b)
    return np.array(positions, dtype=float)

def lj_forces(positions, box_length, epsilon=1.0, sigma=1.0, r_cut=2.5, shifted=True):
    """
    Compute Lennard-Jones forces and potential energy using an O(N^2) algorithm.
    """
    n_atoms = positions.shape[0]
    forces = np.zeros_like(positions)
    potential = 0.0

    r_cut_sigma = r_cut * sigma
    r_cut2 = r_cut_sigma ** 2

    if shifted:
        inv_rc2 = (sigma ** 2) / r_cut2
        inv_rc6 = inv_rc2 ** 3
        inv_rc12 = inv_rc6 ** 2
        u_shift = 4.0 * epsilon * (inv_rc12 - inv_rc6)
    else:
        u_shift = 0.0

    for i in range(n_atoms - 1):
        for j in range(i + 1, n_atoms):
            rij = positions[j] - positions[i]
            rij -= box_length * np.round(rij / box_length)
            r2 = np.dot(rij, rij)

            if 1e-12 < r2 < r_cut2:
                inv_r2 = (sigma ** 2) / r2
                inv_r6 = inv_r2 ** 3
                inv_r12 = inv_r6 ** 2

                pair_potential = 4.0 * epsilon * (inv_r12 - inv_r6) - u_shift
                potential += pair_potential

                force_scalar = 24.0 * epsilon * (2.0 * inv_r12 - inv_r6) / r2
                fij = force_scalar * rij

                forces[i] -= fij
                forces[j] += fij

    return forces, potential

def run_lj_md(
    n_cells=3,
    rho=0.80,
    temperature=1.0,
    n_steps=5000,
    dt=0.005,
    epsilon=1.0,
    sigma=1.0,
    mass=1.0,
    r_cut=2.5,
    sample_every=10,
    thermostat_interval=100,
    random_seed=123,
    shifted=True,
):
    """
    Run a simple 3D LJ molecular dynamics simulation in reduced units.

    Notes
    -----
    This uses velocity-Verlet integration and optional periodic velocity rescaling.
    Setting thermostat_interval=None gives a simple NVE-style production run.
    """
    rng = np.random.default_rng(random_seed)

    n_atoms = 4 * n_cells**3
    number_density = rho / (sigma**3)
    box_length = (n_atoms / number_density) ** (1.0 / 3.0)

    positions = make_fcc_positions(n_cells, box_length)
    velocities = rng.normal(0.0, np.sqrt(temperature / mass), size=(n_atoms, 3))
    velocities -= velocities.mean(axis=0, keepdims=True)

    kinetic = 0.5 * mass * np.sum(velocities**2)
    current_temp = 2.0 * kinetic / (3.0 * n_atoms - 3.0)
    velocities *= np.sqrt(temperature / current_temp)

    forces, potential = lj_forces(
        positions, box_length, epsilon=epsilon, sigma=sigma, r_cut=r_cut, shifted=shifted
    )

    wrapped_traj = []
    unwrapped_traj = []
    times = []
    temperatures = []
    kinetic_energies = []
    potential_energies = []
    total_energies = []

    image = np.zeros_like(positions, dtype=int)

    for step in range(n_steps):
        velocities += 0.5 * dt * forces / mass
        positions += dt * velocities

        crossed = np.floor_divide(positions, box_length).astype(int)
        image += crossed
        positions -= crossed * box_length

        new_forces, potential = lj_forces(
            positions, box_length, epsilon=epsilon, sigma=sigma, r_cut=r_cut, shifted=shifted
        )
        velocities += 0.5 * dt * new_forces / mass
        forces = new_forces

        if thermostat_interval is not None and thermostat_interval > 0 and (step + 1) % thermostat_interval == 0:
            kinetic = 0.5 * mass * np.sum(velocities**2)
            current_temp = 2.0 * kinetic / (3.0 * n_atoms - 3.0)
            velocities *= np.sqrt(temperature / current_temp)

        kinetic = 0.5 * mass * np.sum(velocities**2)
        current_temp = 2.0 * kinetic / (3.0 * n_atoms - 3.0)
        total_energy = kinetic + potential

        if step % sample_every == 0:
            wrapped_traj.append(positions.copy())
            unwrapped_traj.append(positions + image * box_length)
            times.append(step * dt)
            temperatures.append(current_temp)
            kinetic_energies.append(kinetic)
            potential_energies.append(potential)
            total_energies.append(total_energy)

    return {
        "positions": positions,
        "velocities": velocities,
        "wrapped_traj": np.array(wrapped_traj),
        "unwrapped_traj": np.array(unwrapped_traj),
        "times": np.array(times),
        "temperatures": np.array(temperatures),
        "kinetic_energies": np.array(kinetic_energies),
        "potential_energies": np.array(potential_energies),
        "total_energies": np.array(total_energies),
        "box_length": box_length,
        "n_atoms": n_atoms,
        "rho": rho,
        "temperature_target": temperature,
        "dt": dt,
        "sample_every": sample_every,
        "epsilon": epsilon,
        "sigma": sigma,
        "mass": mass,
    }

def write_pdb_topology(filename, positions, box_length, atom_name="Ar", resname="LJ"):
    """
    Write a minimal PDB file for the LJ particles.
    """
    with open(filename, "w") as f:
        f.write("TITLE     LJ FLUID\n")
        f.write(
            f"CRYST1{box_length:9.3f}{box_length:9.3f}{box_length:9.3f}"
            f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1\n"
        )
        for i, xyz in enumerate(positions, start=1):
            x, y, z = xyz
            f.write(
                f"ATOM  {i:5d} {atom_name:<4s} {resname:>3s} A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom_name[:2]:>2s}\n"
            )
        f.write("END\n")

def radial_distribution_function(traj, box_length, r_max=None, n_bins=100):
    """
    Compute g(r) from a wrapped trajectory of a one-component fluid.

    Parameters
    ----------
    traj : ndarray, shape (n_frames, N, 3)
        Wrapped coordinates in a cubic periodic box.
    box_length : float
        Box length.
    r_max : float or None
        Maximum radius for g(r). Defaults to L/2.
    n_bins : int
        Number of histogram bins.
    """
    n_frames, n_atoms, _ = traj.shape
    volume = box_length**3
    number_density = n_atoms / volume

    if r_max is None:
        r_max = box_length / 2.0

    edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.zeros(n_bins, dtype=float)

    for frame in traj:
        for i in range(n_atoms - 1):
            rij = frame[i + 1 :] - frame[i]
            rij -= box_length * np.round(rij / box_length)
            distances = np.linalg.norm(rij, axis=1)
            hist, _ = np.histogram(distances, bins=edges)
            counts += hist

    r_centers = 0.5 * (edges[:-1] + edges[1:])
    shell_volumes = (4.0 / 3.0) * np.pi * (edges[1:]**3 - edges[:-1]**3)

    normalization = n_frames * n_atoms * 0.5 * number_density * shell_volumes
    g_r = counts / normalization
    return r_centers, g_r, counts

def estimate_diffusion_constant(times, msd, fit_start_fraction=0.5, fit_end_fraction=1.0):
    """
    Fit the linear regime of the MSD and return the slope and diffusion constant.
    """
    n = len(times)
    i0 = int(fit_start_fraction * n)
    i1 = max(i0 + 2, int(fit_end_fraction * n))

    x_fit = times[i0:i1]
    y_fit = msd[i0:i1]

    slope, intercept = np.polyfit(x_fit, y_fit, 1)
    diffusion_constant = slope / 6.0
    return {
        "slope": slope,
        "intercept": intercept,
        "D": diffusion_constant,
        "fit_slice": slice(i0, i1),
        "x_fit": x_fit,
        "y_fit": y_fit,
    }

def analyze_lj_with_mdanalysis(
    topology_file="lj_topology.pdb",
    wrapped_traj_file="lj_wrapped.dcd",
    unwrapped_traj_file="lj_unwrapped.dcd",
    rdf_range=None,
    rdf_bins=100,
    msd_fft=True,
):
    """
    Analyze an LJ simulation with MDAnalysis.

    Requires MDAnalysis to be installed in the notebook environment.
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis.rdf import InterRDF
    from MDAnalysis.analysis.msd import EinsteinMSD

    uw = mda.Universe(topology_file, wrapped_traj_file)
    uu = mda.Universe(topology_file, unwrapped_traj_file)

    group_w = uw.atoms
    group_u = uu.atoms

    if rdf_range is None:
        box_length = uw.dimensions[0]
        rdf_range = (0.0, box_length / 2.0)

    rdf = InterRDF(group_w, group_w, nbins=rdf_bins, range=rdf_range, exclusion_block=(1, 1))
    rdf.run()

    msd = EinsteinMSD(group_u, select="all", msd_type="xyz", fft=msd_fft)
    msd.run()

    return {
        "rdf_bins": rdf.results.bins,
        "rdf": rdf.results.rdf,
        "msd_timeseries": msd.results.timeseries,
        "n_frames": len(uw.trajectory),
        "box_length": uw.dimensions[0],
    }
