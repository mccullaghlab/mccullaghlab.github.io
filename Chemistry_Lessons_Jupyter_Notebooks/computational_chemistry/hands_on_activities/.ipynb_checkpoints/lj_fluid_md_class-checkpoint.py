"""Simple class-based teaching library for 3D Lennard-Jones molecular dynamics."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np


@dataclass
class LennardJonesMD:
    """
    Minimal 3D Lennard-Jones molecular dynamics engine in reduced units.

    Parameters
    ----------
    n_cells : int
        Number of FCC unit cells in each dimension. Total atoms = 4 * n_cells**3.
    rho : float
        Reduced number density rho* = N sigma^3 / V.
    temperature : float
        Target reduced temperature.
    n_steps : int
        Number of MD steps.
    dt : float
        Time step in reduced units.
    epsilon, sigma, mass : float
        Lennard-Jones parameters and particle mass.
    r_cut : float
        Cutoff in units of sigma.
    sample_every : int
        Save trajectory/thermo information every this many steps.
    thermostat_interval : int | None
        If not None, apply simple velocity rescaling every thermostat_interval steps.
        Set to None for an NVE-style run after initialization.
    random_seed : int | None
        Seed for random velocity initialization.
    shifted : bool
        If True, use a shifted potential so U(r_cut)=0.
    """

    n_cells: int = 3
    rho: float = 0.80
    temperature: float = 1.0
    n_steps: int = 5000
    dt: float = 0.005
    epsilon: float = 1.0
    sigma: float = 1.0
    mass: float = 1.0
    r_cut: float = 2.5
    sample_every: int = 10
    thermostat_interval: int | None = 100
    random_seed: int | None = 123
    shifted: bool = True

    positions: np.ndarray | None = field(default=None, init=False, repr=False)
    velocities: np.ndarray | None = field(default=None, init=False, repr=False)
    forces: np.ndarray | None = field(default=None, init=False, repr=False)
    box_length: float | None = field(default=None, init=False)
    n_atoms: int | None = field(default=None, init=False)
    image: np.ndarray | None = field(default=None, init=False, repr=False)

    wrapped_traj: list[np.ndarray] = field(default_factory=list, init=False, repr=False)
    unwrapped_traj: list[np.ndarray] = field(default_factory=list, init=False, repr=False)
    times: list[float] = field(default_factory=list, init=False)
    temperatures: list[float] = field(default_factory=list, init=False, repr=False)
    kinetic_energies: list[float] = field(default_factory=list, init=False, repr=False)
    potential_energies: list[float] = field(default_factory=list, init=False, repr=False)
    total_energies: list[float] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.n_cells < 1:
            raise ValueError("n_cells must be >= 1")
        if self.rho <= 0 or self.temperature <= 0 or self.dt <= 0:
            raise ValueError("rho, temperature, and dt must be positive")
        if self.sample_every < 1:
            raise ValueError("sample_every must be >= 1")
        if self.r_cut <= 0:
            raise ValueError("r_cut must be positive")

        self.n_atoms = 4 * self.n_cells**3
        number_density = self.rho / (self.sigma**3)
        self.box_length = (self.n_atoms / number_density) ** (1.0 / 3.0)
        self._initialize_state()

    def _make_fcc_positions(self) -> np.ndarray:
        basis = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ],
            dtype=float,
        )
        a = self.box_length / self.n_cells
        positions = []
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                for k in range(self.n_cells):
                    origin = np.array([i, j, k], dtype=float) * a
                    for b in basis:
                        positions.append(origin + a * b)
        return np.array(positions, dtype=float)

    def _initialize_state(self) -> None:
        rng = np.random.default_rng(self.random_seed)
        self.positions = self._make_fcc_positions()
        self.velocities = rng.normal(
            0.0, np.sqrt(self.temperature / self.mass), size=(self.n_atoms, 3)
        )
        self.velocities -= self.velocities.mean(axis=0, keepdims=True)

        kinetic = 0.5 * self.mass * np.sum(self.velocities**2)
        current_temp = 2.0 * kinetic / (3.0 * self.n_atoms - 3.0)
        self.velocities *= np.sqrt(self.temperature / current_temp)

        self.image = np.zeros_like(self.positions, dtype=int)
        self.forces, _ = self._lj_forces(self.positions)
        self._clear_samples()

    def _clear_samples(self) -> None:
        self.wrapped_traj = []
        self.unwrapped_traj = []
        self.times = []
        self.temperatures = []
        self.kinetic_energies = []
        self.potential_energies = []
        self.total_energies = []

    def _lj_forces(self, positions: np.ndarray) -> tuple[np.ndarray, float]:
        forces = np.zeros_like(positions)
        potential = 0.0

        r_cut_sigma = self.r_cut * self.sigma
        r_cut2 = r_cut_sigma**2

        if self.shifted:
            inv_rc2 = (self.sigma**2) / r_cut2
            inv_rc6 = inv_rc2**3
            inv_rc12 = inv_rc6**2
            u_shift = 4.0 * self.epsilon * (inv_rc12 - inv_rc6)
        else:
            u_shift = 0.0

        for i in range(self.n_atoms - 1):
            for j in range(i + 1, self.n_atoms):
                rij = positions[j] - positions[i]
                rij -= self.box_length * np.round(rij / self.box_length)
                r2 = np.dot(rij, rij)

                if 1e-12 < r2 < r_cut2:
                    inv_r2 = (self.sigma**2) / r2
                    inv_r6 = inv_r2**3
                    inv_r12 = inv_r6**2

                    potential += 4.0 * self.epsilon * (inv_r12 - inv_r6) - u_shift

                    force_scalar = 24.0 * self.epsilon * (2.0 * inv_r12 - inv_r6) / r2
                    fij = force_scalar * rij
                    forces[i] -= fij
                    forces[j] += fij

        return forces, potential

    def _instantaneous_kinetic_temperature(self) -> tuple[float, float]:
        kinetic = 0.5 * self.mass * np.sum(self.velocities**2)
        current_temp = 2.0 * kinetic / (3.0 * self.n_atoms - 3.0)
        return kinetic, current_temp

    def _store_sample(self, step: int, potential: float) -> None:
        kinetic, current_temp = self._instantaneous_kinetic_temperature()
        total_energy = kinetic + potential

        self.wrapped_traj.append(self.positions.copy())
        self.unwrapped_traj.append(self.positions + self.image * self.box_length)
        self.times.append(step * self.dt)
        self.temperatures.append(current_temp)
        self.kinetic_energies.append(kinetic)
        self.potential_energies.append(potential)
        self.total_energies.append(total_energy)

    def run(self, reset: bool = True) -> dict[str, np.ndarray | float | int]:
        """
        Run the MD simulation.

        Parameters
        ----------
        reset : bool
            If True, reinitialize positions and velocities before running.

        Returns
        -------
        dict
            Dictionary with trajectory, energies, and simulation metadata.
        """
        if reset:
            self._initialize_state()

        potential = self._lj_forces(self.positions)[1]
        self._store_sample(step=0, potential=potential)

        for step in range(1, self.n_steps + 1):
            self.velocities += 0.5 * self.dt * self.forces / self.mass
            self.positions += self.dt * self.velocities

            crossed = np.floor_divide(self.positions, self.box_length).astype(int)
            self.image += crossed
            self.positions -= crossed * self.box_length

            new_forces, potential = self._lj_forces(self.positions)
            self.velocities += 0.5 * self.dt * new_forces / self.mass
            self.forces = new_forces

            if (
                self.thermostat_interval is not None
                and self.thermostat_interval > 0
                and step % self.thermostat_interval == 0
            ):
                kinetic, current_temp = self._instantaneous_kinetic_temperature()
                self.velocities *= np.sqrt(self.temperature / current_temp)

            if step % self.sample_every == 0:
                self._store_sample(step=step, potential=potential)

        return self.results()

    def results(self) -> dict[str, np.ndarray | float | int]:
        """Return sampled data as numpy arrays and metadata."""
        return {
            "positions": self.positions.copy(),
            "velocities": self.velocities.copy(),
            "wrapped_traj": np.array(self.wrapped_traj),
            "unwrapped_traj": np.array(self.unwrapped_traj),
            "times": np.array(self.times),
            "temperatures": np.array(self.temperatures),
            "kinetic_energies": np.array(self.kinetic_energies),
            "potential_energies": np.array(self.potential_energies),
            "total_energies": np.array(self.total_energies),
            "box_length": float(self.box_length),
            "n_atoms": int(self.n_atoms),
            "rho": float(self.rho),
            "temperature_target": float(self.temperature),
            "dt": float(self.dt),
            "sample_every": int(self.sample_every),
            "epsilon": float(self.epsilon),
            "sigma": float(self.sigma),
            "mass": float(self.mass),
        }

    def write_pdb(
        self,
        filename: str | Path,
        positions: np.ndarray | None = None,
        atom_name: str = "Ar",
        resname: str = "LJ",
    ) -> None:
        """Write a minimal PDB topology file for the current or supplied positions."""
        xyz = self.positions if positions is None else positions
        filename = Path(filename)
        with filename.open("w") as f:
            f.write("TITLE     LJ FLUID\n")
            f.write(
                f"CRYST1{self.box_length:9.3f}{self.box_length:9.3f}{self.box_length:9.3f}"
                f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1\n"
            )
            for i, (x, y, z) in enumerate(xyz, start=1):
                f.write(
                    f"ATOM  {i:5d} {atom_name:<4s} {resname:>3s} A{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom_name[:2]:>2s}\n"
                )
            f.write("END\n")

    def write_trajectory(
        self,
        filename: str | Path,
        wrapped: bool = True,
        atom_name: str = "Ar",
        fmt: str | None = None,
    ) -> None:
        """
        Write trajectory data.

        Supported formats
        -----------------
        xyz : multi-frame XYZ text trajectory
        npz : compressed NumPy archive containing trajectory, times, and box length

        The format is inferred from the filename extension unless `fmt` is supplied.
        """
        filename = Path(filename)
        traj = np.array(self.wrapped_traj if wrapped else self.unwrapped_traj)
        times = np.array(self.times)

        out_format = (fmt or filename.suffix.lstrip(".")).lower()
        if out_format == "xyz":
            with filename.open("w") as f:
                for frame, time in zip(traj, times):
                    f.write(f"{self.n_atoms}\n")
                    f.write(f"LJ fluid frame; time={time:.6f}; box={self.box_length:.6f}\n")
                    for x, y, z in frame:
                        f.write(f"{atom_name} {x:.8f} {y:.8f} {z:.8f}\n")
        elif out_format == "npz":
            np.savez_compressed(
                filename,
                trajectory=traj,
                times=times,
                box_length=self.box_length,
                wrapped=wrapped,
            )
        else:
            raise ValueError("Unsupported trajectory format. Use '.xyz' or '.npz'.")
