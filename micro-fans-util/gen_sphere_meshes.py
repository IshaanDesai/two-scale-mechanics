import numpy as np
import h5py



def write_h5(grid, filename):
    with h5py.File(filename, "w") as f:
        f.create_dataset("sphere/32x32x32/ms", data=grid)
    return grid

def mark_sphere_cells(N: int, K: int, r: float) -> np.ndarray:
    """
    Discretize the unit cube into N^3 cells and mark cells as 1
    if at least K of their 8 corners lie inside a sphere centered
    at (0.5, 0.5, 0.5) with radius r.

    Parameters
    ----------
    N : int
        Number of cells per dimension.
    K : int
        Minimum number of corners inside the sphere (0 <= K <= 8).
    r : float
        Sphere radius.
    filename : str
        Output HDF5 file path.

    Returns
    -------
    grid : np.ndarray
        Binary array of shape (N, N, N).
    """
    if not (0 <= K <= 8):
        raise ValueError("K must satisfy 0 <= K <= 8")

    # Cell corner coordinates along one axis
    coords = np.linspace(0.0, 1.0, N + 1, dtype=np.float64)

    # Generate 3D grid of corner coordinates
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")

    # Compute squared distance to sphere center
    dx = X - 0.5
    dy = Y - 0.5
    dz = Z - 0.5
    dist2 = dx * dx + dy * dy + dz * dz

    inside = dist2 <= r * r  # Boolean array shape (N+1, N+1, N+1)

    # Count how many of the 8 corners of each cell are inside
    corner_count = (
        inside[:-1, :-1, :-1].astype(np.uint8) +
        inside[1:, :-1, :-1].astype(np.uint8) +
        inside[:-1, 1:, :-1].astype(np.uint8) +
        inside[:-1, :-1, 1:].astype(np.uint8) +
        inside[1:, 1:, :-1].astype(np.uint8) +
        inside[1:, :-1, 1:].astype(np.uint8) +
        inside[:-1, 1:, 1:].astype(np.uint8) +
        inside[1:, 1:, 1:].astype(np.uint8)
    )

    grid = (corner_count >= K).astype(np.uint8)
    return grid

def mark_sphere_cells_blocked(
    N: int,
    K: int,
    r: float,
    block_size: int = 64
) -> np.ndarray:
    """
    Memory-efficient blocked version of mark_sphere_cells.

    Processes the cube in z-blocks to avoid allocating full (N+1)^3 arrays.
    Suitable for N >= 1024.

    Parameters
    ----------
    N : int
        Number of cells per dimension.
    K : int
        Minimum number of corners inside sphere (0 <= K <= 8).
    r : float
        Sphere radius.
    filename : str
        Output HDF5 file.
    block_size : int
        Number of cell layers processed per block (e.g. 32 or 64).

    Returns
    -------
    grid : np.ndarray
        Binary array of shape (N, N, N), dtype uint8.
    """
    if not (0 <= K <= 8):
        raise ValueError("K must satisfy 0 <= K <= 8")

    if N % 1 != 0:
        raise ValueError("N must be integer")

    # Final output buffer (allowed by requirement)
    grid = np.empty((N, N, N), dtype=np.uint8)

    coords = np.linspace(0.0, 1.0, N + 1, dtype=np.float64)
    dx = coords - 0.5
    dx2 = dx * dx
    r2 = r * r

    # Precompute x/y 2D part once (N+1 x N+1)
    x2 = dx2[:, None]
    y2 = dx2[None, :]
    xy2 = x2 + y2  # shape (N+1, N+1)

    for z_start in range(0, N, block_size):
        z_end = min(z_start + block_size, N)

        # We need +1 layer for corners
        z_slice = slice(z_start, z_end + 1)

        dz2 = dx2[z_slice]
        z2 = dz2[None, None, :]  # broadcast

        # inside shape: (N+1, N+1, block+1)
        inside = (xy2[:, :, None] + z2) <= r2

        # Corner aggregation
        corner_count = (
            inside[:-1, :-1, :-1].astype(np.uint8) +
            inside[1:, :-1, :-1].astype(np.uint8) +
            inside[:-1, 1:, :-1].astype(np.uint8) +
            inside[:-1, :-1, 1:].astype(np.uint8) +
            inside[1:, 1:, :-1].astype(np.uint8) +
            inside[1:, :-1, 1:].astype(np.uint8) +
            inside[:-1, 1:, 1:].astype(np.uint8) +
            inside[1:, 1:, 1:].astype(np.uint8)
        )

        grid[:, :, z_start:z_end] = (corner_count >= K).astype(np.uint8)

    return grid

def downsample_majority(grid: np.ndarray) -> np.ndarray:
    """
    Downsample a binary 3D grid by a factor of 2 in each dimension.
    A coarse cell becomes 1 if >= 4 of the corresponding 8 fine cells are 1.

    Parameters
    ----------
    grid : np.ndarray
        Binary array of shape (N, N, N) with N even.

    Returns
    -------
    coarse : np.ndarray
        Downsampled binary array of shape (N//2, N//2, N//2).
    """
    if grid.ndim != 3:
        raise ValueError("Input grid must be 3-dimensional")

    N = grid.shape[0]
    if grid.shape[1] != N or grid.shape[2] != N:
        raise ValueError("Input grid must be cubic (N, N, N)")

    if N % 2 != 0:
        raise ValueError("Grid size N must be even")

    # Ensure integer type for summation
    fine = grid.astype(np.uint8, copy=False)

    # Reshape into (N/2, 2, N/2, 2, N/2, 2)
    coarse_counts = fine.reshape(
        N // 2, 2,
        N // 2, 2,
        N // 2, 2
    ).sum(axis=(1, 3, 5))

    # Majority rule: >= 4 out of 8
    coarse = (coarse_counts >= 4).astype(np.uint8)

    return coarse

def downsample_to_resolution(grid: np.ndarray, N_target: int) -> np.ndarray:
    """
    Downsample a cubic binary grid to a target resolution using majority rule.

    A coarse cell becomes 1 if >= 50% of the corresponding fine cells are 1.

    Parameters
    ----------
    grid : np.ndarray
        Binary array of shape (N, N, N)
    N_target : int
        Desired output resolution (must divide N)

    Returns
    -------
    coarse : np.ndarray
        Binary array of shape (N_target, N_target, N_target)
    """
    if grid.ndim != 3:
        raise ValueError("Input must be 3D")

    N = grid.shape[0]
    if grid.shape != (N, N, N):
        raise ValueError("Grid must be cubic")

    if N % N_target != 0:
        raise ValueError("N_target must divide N")

    factor = N // N_target

    fine = grid.astype(np.uint32, copy=False)

    # reshape into blocks
    reshaped = fine.reshape(
        N_target, factor,
        N_target, factor,
        N_target, factor
    )

    # sum inside each block
    block_sum = reshaped.sum(axis=(1, 3, 5))

    # majority threshold
    threshold = (factor ** 3) / 2

    coarse = (block_sum >= threshold).astype(np.uint8)

    return coarse

def downsample_fractional_volume(grid: np.ndarray, N_target: int) -> np.ndarray:
    """
    Downsample a cubic binary grid to target resolution while exactly preserving total volume.
    """
    if grid.ndim != 3:
        raise ValueError("Input grid must be 3D cubic array")

    N = grid.shape[0]
    if N % N_target != 0:
        raise ValueError("Target resolution must divide grid size")

    factor = N // N_target
    fine = grid.astype(np.uint32, copy=False)

    # Compute sum of ones in each coarse block
    reshaped = fine.reshape(
        N_target, factor,
        N_target, factor,
        N_target, factor
    )
    block_sum = reshaped.sum(axis=(1, 3, 5))

    flat_sum = block_sum.flatten()
    total_ones = int(fine.sum())  # <-- convert to Python int

    coarse_flat = np.zeros_like(flat_sum, dtype=np.uint8)

    if total_ones > 0:
        # ensure total_ones does not exceed number of coarse blocks
        n_blocks = flat_sum.size
        n_select = min(total_ones, n_blocks)

        # select top n_select blocks
        top_indices = np.argpartition(-flat_sum, n_select - 1)[:n_select]
        coarse_flat[top_indices] = 1

    coarse = coarse_flat.reshape((N_target, N_target, N_target))
    return coarse

def select_optimal_blocks(grid: np.ndarray, N_target: int) -> np.ndarray:
    """
    Select coarse blocks to preserve volume using a greedy fraction-error approach.
    """
    N = grid.shape[0]
    factor = N // N_target
    fine = grid.astype(np.uint32, copy=False)

    # Sum of ones per coarse block
    reshaped = fine.reshape(
        N_target, factor,
        N_target, factor,
        N_target, factor
    )
    block_sum = reshaped.sum(axis=(1, 3, 5))

    total_ones = fine.sum()

    # Fraction of ones per block
    fractions = block_sum / total_ones

    # Fraction of zeros (error)
    block_size = factor ** 3
    errors = block_size - block_sum
    error_fractions = errors / total_ones

    # Score = maximize fraction, minimize error fraction
    score = fractions - error_fractions  # = (2*block_sum - block_size)/total_ones

    flat_score = score.flatten()

    # Number of blocks we can set to 1 to preserve total volume
    # total ones in fine grid / block_size, rounded
    n_set = int(round(total_ones / block_size))+1

    # Indices of top scoring blocks
    top_indices = np.argpartition(-flat_score, n_set - 1)[:n_set]

    # Build coarse grid
    coarse_flat = np.zeros_like(flat_score, dtype=np.uint8)
    coarse_flat[top_indices] = 1
    coarse = coarse_flat.reshape((N_target, N_target, N_target))

    return coarse

def select_optimal_blocks_symmetric(grid: np.ndarray, N_target: int) -> np.ndarray:
    """
    Downsample a cubic binary grid to target resolution using
    volume-preserving fraction/error method, but only on the lower quadrant.
    Results are mirrored to enforce symmetry.

    Parameters
    ----------
    grid : np.ndarray
        Binary 3D array of shape (N, N, N)
    N_target : int
        Target coarse resolution

    Returns
    -------
    coarse : np.ndarray
        Coarse binary array of shape (N_target, N_target, N_target)
        Symmetric across all axes.
    """
    N = grid.shape[0]
    factor = N // N_target
    fine = grid.astype(np.uint32, copy=False)

    # Compute coarse block sums for the lower octant
    half = N_target // 2 + N_target % 2  # include center if odd

    reshaped = fine.reshape(
        N_target, factor,
        N_target, factor,
        N_target, factor
    )
    block_sum = reshaped.sum(axis=(1, 3, 5))

    # Take only lower quadrant
    lower_octant = block_sum[:half, :half, :half]

    total_ones = fine.sum()
    block_size = factor ** 3

    # Fraction and error fraction
    fractions = lower_octant / total_ones
    errors = block_size - lower_octant
    error_fractions = errors / total_ones

    score = fractions - error_fractions
    flat_score = score.flatten()

    # Number of coarse blocks to set in octant
    n_set = int(round(total_ones / 8 / block_size)) + 1

    top_indices = np.argpartition(-flat_score, n_set - 1)[:n_set]

    # Build octant coarse grid
    octant_flat = np.zeros_like(flat_score, dtype=np.uint8)
    octant_flat[top_indices] = 1
    octant = octant_flat.reshape(lower_octant.shape)

    # Now mirror along axes to full cube
    coarse = np.zeros((N_target, N_target, N_target), dtype=np.uint8)

    # Fill the lower octant
    coarse[:half, :half, :half] = octant

    # Mirror across x-axis
    coarse[-half:, :half, :half] = octant[::-1, :, :]
    # Mirror across y-axis
    coarse[:, -half:, :half] = coarse[:, :half, :half][:, ::-1, :]
    # Mirror across z-axis
    coarse[:, :, -half:] = coarse[:, :, :half][:, :, ::-1]

    return coarse

def compute_vol_frac(grid):
    N = grid.shape[0]
    cube_vol = np.power(1.0 / N, 3)
    return np.sum(grid) * cube_vol

def compute_true_frac(rad):
    return 4.0 / 3.0 * np.pi * np.power(rad, 3)

# ==================================================================
# Call this to generate the sphere meshes
# can switch select_optimal_blocks_symmetric to other implementation
# ==================================================================
def generate_spheres():
    with h5py.File("./ref_sphere32.h5", "r") as f:
        ref_data = f["sphere"]["32x32x32"]["ms"][...]

    rad = 0.4
    grid1024 = mark_sphere_cells_blocked(1024, 4, rad)
    grid512 = select_optimal_blocks_symmetric(grid1024, 512)
    grid256 = select_optimal_blocks_symmetric(grid512, 256)
    grid128 = select_optimal_blocks_symmetric(grid256, 128)
    grid64 = select_optimal_blocks_symmetric(grid128, 64)
    grid32 = select_optimal_blocks_symmetric(grid64, 32)

    print(f"True Volume: {compute_true_frac(rad)}")
    print(f"1024 Volume: {compute_vol_frac(grid1024)}")
    print(f"512 Volume: {compute_vol_frac(grid512)}")
    print(f"256 Volume: {compute_vol_frac(grid256)}")
    print(f"128 Volume: {compute_vol_frac(grid128)}")
    print(f"64 Volume: {compute_vol_frac(grid64)}")
    print(f"32 Volume: {compute_vol_frac(grid32)}")
    print(f"ref Volume: {compute_vol_frac(ref_data)}")

    write_h5(grid512, "./sphere512.h5")
    write_h5(grid256, "./sphere256.h5")
    write_h5(grid128, "./sphere128.h5")
    write_h5(grid64, "./sphere64.h5")
    write_h5(grid32, "./sphere32.h5")
