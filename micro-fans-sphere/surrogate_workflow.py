import os
import subprocess
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np
import matplotlib.pyplot as plt


def create_snapshots() -> None:
    """
    Create snapshots of the FANS model in micro_fans/ using the Micro Manager.
    The snapshots are saved in the specified directory.
    """
    rng = np.random.default_rng(42)
    std = np.array([
        0.05, 0.02, 0.0015, 0.0020, 0.015, 0.02
    ])
    mean = np.zeros(6)

    # Draw from independent Gaussians per dimension
    samples = rng.normal(loc=mean, scale=std, size=(100, 6))
    e123 = samples[:, 0:3]
    e456 = samples[:, 3:6]

    with h5py.File("input_samples.hdf5", "w") as f:
        f.create_dataset("strains1to3", data=e123)
        f.create_dataset("strains4to6", data=e456)

    # Run the Micro Manager to create snapshot database
    subprocess.call('micro-manager-precice --snapshot micro-manager-snapshot-config.json', shell=True)


def read_snapshots(snapshots_dir: str = None) -> tuple:
    """
    Read the snapshots created by the Micro Manager and return the inputs and outputs.
    The inputs are the concentration samples, and the outputs are the porosity and conductivity matrix values.

    Parameters
    ----------
    snapshots_dir : str
        The directory where the snapshots are stored.

    Returns
    -------
    tuple
        A tuple containing the inputs (strain samples) and outputs (porosity and conductivity values).
    """
    if snapshots_dir is None: snapshots_dir = '.'
    with h5py.File(os.path.join(snapshots_dir, "snapshot_data.hdf5"), "r") as f:
        strains1to3 = f['strains1to3'][:]
        strains4to6 = f['strains4to6'][:]
        stresses1to3 = f['stresses1to3'][:]
        stresses4to6 = f['stresses4to6'][:]
        cmat1 = f['cmat1'][:]
        cmat2 = f['cmat2'][:]
        cmat3 = f['cmat3'][:]
        cmat4 = f['cmat4'][:]
        cmat5 = f['cmat5'][:]
        cmat6 = f['cmat6'][:]
        cmat7 = f['cmat7'][:]

    outputs = {
        'stresses1to3': stresses1to3,
        'stresses4to6': stresses4to6,
        'cmat1': cmat1,
        'cmat2': cmat2,
        'cmat3': cmat3,
        'cmat4': cmat4,
        'cmat5': cmat5,
        'cmat6': cmat6,
        'cmat7': cmat7,
        "x_values": np.array([0])
    }
    merge_buffer = np.zeros((strains1to3.shape[0], strains1to3.shape[1] + strains4to6.shape[1]))
    merge_buffer[:, 0:3] = strains1to3
    merge_buffer[:, 3:6] = strains4to6
    return merge_buffer, outputs


def split_samples(X, y, n_valid):
    """
    Split the samples and evaluations into training and validation/test data.
    The split is performed randomly.

    Parameters
    ----------
    X : np.ndarray
        Samples, shape (#samples, #parameters)
    y : dict
        Corresponding model evaluations. Expected to match BVR output format.
    n_valid : int
        Number of samples to keep for validation.

    Returns
    -------
    X_train, y_train,
    X_valid, y_valid
    """
    n_samples = X.shape[0]
    if n_valid >= n_samples:
        raise AttributeError('The set number of validation points is invalid.')

    # Random split
    n_train = n_samples - n_valid
    choice = np.random.choice(
        range(n_samples), size=(n_train,), replace=False
    )
    ind = np.zeros(n_samples, dtype=bool)
    ind[choice] = True

    # Split samples
    X_train = X[ind]
    X_valid = X[~ind]

    # Split outputs
    y_train = {}
    y_valid = {}
    for key in y:
        if key != "x_values":
            y_train[key] = y[key][ind]
            y_valid[key] = y[key][~ind]

    return X_train, y_train, X_valid, y_valid


def create_surrogate(snapshots_dir: str = None) -> tuple:
    """
    Create a surrogate model from the Micro Manager snapshots.

    """
    # We create a fake model from model.py because we directly provide the
    # input and outputs from the Micro Manager snapshots.
    model = PyLinkForwardModel()
    model.py_file = "model"
    model.name = "micro-fans-surrogate"
    model.link_type = "function"
    model.output_names = ['stresses1to3', 'stresses4to6', 'cmat1', 'cmat2', 'cmat3', 'cmat4', 'cmat5', 'cmat6', 'cmat7']

    x, y = read_snapshots(snapshots_dir)

    # Split the samples into training and validation sets
    n_valid = int(0.4 * x.shape[0])
    x_train, y_train, x_valid, y_valid = split_samples(x, y, n_valid)

    #0.05, 0.02, 0.0015, 0.0020, 0.015, 0.02
    inputs = Input()
    inputs.add_marginals()
    inputs.marginals[0].name = '$X_1$'
    inputs.marginals[0].dist_type = 'norm'
    inputs.marginals[0].parameters = [0, 0.05]
    inputs.add_marginals()
    inputs.marginals[1].name = '$X_2$'
    inputs.marginals[1].dist_type = 'norm'
    inputs.marginals[1].parameters = [0, 0.02]
    inputs.add_marginals()
    inputs.marginals[2].name = '$X_3$'
    inputs.marginals[2].dist_type = 'norm'
    inputs.marginals[2].parameters = [0, 0.0015]
    inputs.add_marginals()
    inputs.marginals[3].name = '$X_4$'
    inputs.marginals[3].dist_type = 'norm'
    inputs.marginals[3].parameters = [0, 0.002]
    inputs.add_marginals()
    inputs.marginals[4].name = '$X_5$'
    inputs.marginals[4].dist_type = 'norm'
    inputs.marginals[4].parameters = [0, 0.015]
    inputs.add_marginals()
    inputs.marginals[5].name = '$X_6$'
    inputs.marginals[5].dist_type = 'norm'
    inputs.marginals[5].parameters = [0, 0.02]

    exp_design = ExpDesigns(inputs)
    exp_design.x = x_train
    exp_design.y = y_train

    # Create the surrogate model
    meta_model = PCE(inputs)
    meta_model.meta_model_type = "aPCE"
    meta_model.pce_reg_method = "FastARD"
    meta_model.pce_deg = 5

    # Train the surrogate model
    engine = Engine(meta_model, model, exp_design)
    engine.train_normal()

    with open(f'{model.name}.pkl', 'wb') as output:
        joblib.dump(engine, output, 2)

    return x_valid, y_valid


def validate_surrogate(x_valid, y_valid, model_name="micro-fans-surrogate.pkl"):
    """
    Validate the surrogate model using the validation samples and outputs.

    Parameters
    ----------
    x_valid : np.ndarray
        Validation samples.
    y_valid : dict
        Corresponding model evaluations for validation samples.
    model_name : str
        Name of the surrogate model file.
    """
    with open(model_name, 'rb') as input:
        engine = joblib.load(input)

    y_metamod, _ = engine.eval_metamodel(x_valid)

    # engine.plot_adapt(y_valid, y_metamod, y_metamod_std, x_valid)

    # Compare predictions with true values
    plt.figure()
    plt.scatter(y_valid["stresses1to3"][:, 0], y_metamod["stresses1to3"][:, 0], c='#FF0000', label='s1')
    plt.scatter(y_valid["stresses1to3"][:, 1], y_metamod["stresses1to3"][:, 1], c='#FF3A00', label='s2')
    plt.scatter(y_valid["stresses1to3"][:, 2], y_metamod["stresses1to3"][:, 2], c='#FF5500', label='s3')
    plt.scatter(y_valid["stresses4to6"][:, 0], y_metamod["stresses4to6"][:, 0], c='#FF7700', label='s4')
    plt.scatter(y_valid["stresses4to6"][:, 1], y_metamod["stresses4to6"][:, 1], c='#FF9900', label='s5')
    plt.scatter(y_valid["stresses4to6"][:, 2], y_metamod["stresses4to6"][:, 2], c='#FFB300', label='s6')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: stresses")
    plt.show()

    plt.figure()
    plt.scatter(y_valid["cmat1"][:, 0], y_metamod["cmat1"][:, 0], c='#0000FF', label='c1')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat1"][:, 1], y_metamod["cmat1"][:, 1], c='#0A00F5', label='c2')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat1"][:, 2], y_metamod["cmat1"][:, 2], c='#1400EB', label='c3')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat2"][:, 0], y_metamod["cmat2"][:, 0], c='#1E00E1', label='c4')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat2"][:, 1], y_metamod["cmat2"][:, 1], c='#2800D7', label='c5')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat2"][:, 2], y_metamod["cmat2"][:, 2], c='#3200CD', label='c6')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat3"][:, 0], y_metamod["cmat3"][:, 0], c='#3C00C3', label='c7')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat3"][:, 1], y_metamod["cmat3"][:, 1], c='#4600B9', label='c8')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat3"][:, 2], y_metamod["cmat3"][:, 2], c='#5000AF', label='c9')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat4"][:, 0], y_metamod["cmat4"][:, 0], c='#5A00A5', label='c10')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat4"][:, 1], y_metamod["cmat4"][:, 1], c='#64009B', label='c11')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat4"][:, 2], y_metamod["cmat4"][:, 2], c='#6E0091', label='c12')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat5"][:, 0], y_metamod["cmat5"][:, 0], c='#780087', label='c13')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat5"][:, 1], y_metamod["cmat5"][:, 1], c='#82007D', label='c14')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat5"][:, 2], y_metamod["cmat5"][:, 2], c='#8C0073', label='c15')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat6"][:, 0], y_metamod["cmat6"][:, 0], c='#960069', label='c16')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat6"][:, 1], y_metamod["cmat6"][:, 1], c='#A0005F', label='c17')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat6"][:, 2], y_metamod["cmat6"][:, 2], c='#AA0055', label='c18')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat7"][:, 0], y_metamod["cmat7"][:, 0], c='#B4004B', label='c19')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat7"][:, 1], y_metamod["cmat7"][:, 1], c='#BE0041', label='c20')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()
    plt.figure()
    plt.scatter(y_valid["cmat7"][:, 2], y_metamod["cmat7"][:, 2], c='#C80037', label='c21')
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: c-mat")
    plt.show()

def main():
    create_snapshots()
    x_valid, y_valid = create_surrogate()
    validate_surrogate(x_valid, y_valid)


if __name__ == "__main__":
    main()
