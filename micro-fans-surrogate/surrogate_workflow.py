import os
import subprocess
import json
import h5py
import joblib
import numpy as np
import shutil

from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import matplotlib.pyplot as plt


def run_meso():
    meso_config = {
        "output_path": "output/surr_wf_notch",
        "mesh": {
            "path": "examples/coupled-notch/notch.msh",
            "bc_dc_use_tag": False,
            "bc_dc_locator": "plane_xz_low",
            "bc_nm_use_tag": False,
            "bc_nm_locator": "plane_xz_high",
            "bc_nm_dim": 1,
            "bc_nm_value": 0.01,
        },
        "problem": {
            "lambda": 10.0,
            "mu": 5.0,
            "alpha": 100.0,
            "strain_type": "small_strain",
            "elem_degree": 1,
        },
        "simulation": {
            "type": "CoupledSim",
            "micro_type": "PYFANS",
            "write_state": "output/surr_wf_notch",
            "write_state_type": ["E"],
            "init_with_micro": True,
            "precice_xml_path": "../../precice-config-fans-small-strain.xml",
            "slurm_id": "default",
        },
    }

    meso_path = "../meso-fenics"
    config_name = "config-surrogate-meso.json"
    meso_config_path = f"{meso_path}/examples/{config_name}"
    if os.path.exists(meso_config_path):
        os.remove(meso_config_path)
    with open(meso_config_path, "w") as f:
        json.dump(meso_config, f)

    subprocess.Popen(
        "(cd ../micro-fans-notch-sphere/ && mpiexec -n 8 micro-manager-precice micro-manager-pyfans-config-stateless.json)",
        shell=True,
    )
    cmd = f"(cd {meso_path} && python macro.py examples/{config_name})"
    proc = subprocess.run(cmd, shell=True)
    proc.check_returncode()

    if os.path.exists(meso_config_path):
        os.remove(meso_config_path)

    with h5py.File(f"{meso_path}/output/surr_wf_notch_0_1.h5", "r") as f:
        strains = f["strain_data"][...]

    def safe_remove(path):
        if os.path.exists(path):
            if os.path.isdir(path):
                subprocess.run(f"rm -rf {path}", shell=True)
            else:
                os.remove(path)

    safe_remove(f"{meso_path}/output/surr_wf_notch_0_0.h5")
    safe_remove(f"{meso_path}/output/surr_wf_notch_0_0.bp")
    safe_remove(f"{meso_path}/output/surr_wf_notch_0_1.h5")
    safe_remove(f"{meso_path}/output/surr_wf_notch_0_1.bp")
    safe_remove(f"{meso_path}/output/surr_wf_notch_0_init-micro-sol.h5")
    safe_remove(f"{meso_path}/output/surr_wf_notch_0_init-micro-sol.bp")

    return strains


def generate_snapshot_data(meso_strain, n_samples):
    mean = np.mean(meso_strain, axis=0)
    std = np.std(meso_strain, axis=0)
    cov = np.cov(meso_strain.T)

    rand = np.random.default_rng(42)
    samples = rand.multivariate_normal(mean, cov, n_samples)

    with h5py.File(f"snapshot_input.h5", "w") as f:
        f.create_dataset("strains1to3", data=samples[:, 0:3])
        f.create_dataset("strains4to6", data=samples[:, 3:6])

    with open("distributions.json", "w") as f:
        json.dump({"mean": mean.tolist(), "std": std.tolist()}, f)


def run_snapshot():
    snapshot_config = {
        "micro_file_name": "PyFANS",
        "coupling_params": {
            "parameter_file_name": "snapshot_input.h5",
            "read_data_names": ["strains1to3", "strains4to6"],
            "write_data_names": [
                "stresses1to3",
                "stresses4to6",
                "cmat1",
                "cmat2",
                "cmat3",
                "cmat4",
                "cmat5",
                "cmat6",
                "cmat7",
            ],
        },
        "simulation_params": {"micro_dt": 0.01},
        "snapshot_params": {
            "initialize_once": True,
            "output_file_name": "snapshot_output.h5",
        },
        "output_directory": ".",
    }

    with open("snapshot_config.json", "w") as f:
        json.dump(snapshot_config, f)

    if not os.path.exists("./PyFANS.so"):
        shutil.copy("../micro-fans-notch-sphere/PyFANS.so", "./PyFANS.so")

    pyfans_config = {
        "microstructure": {
            "filepath": "sphere32.h5",
            "datasetname": "/sphere/32x32x32/ms",
            "L": [1.0, 1.0, 1.0],
        },
        "problem_type": "mechanical",
        "matmodel": "LinearElasticIsotropic",
        "material_properties": {
            "bulk_modulus": [62.5000, 222.222],
            "shear_modulus": [28.8462, 166.6667],
        },
        "no_mpi": True,
        "method": "cg",
        "error_parameters": {
            "measure": "Linfinity",
            "type": "absolute",
            "tolerance": 1e-10,
        },
        "n_it": 100,
        "macroscale_loading": [[[0, 0, 0, 0, 0, 0]]],
        "results": [],
    }

    with open("input.json", "w") as f:
        json.dump(pyfans_config, f)

    if not os.path.exists("./sphere32.h5"):
        shutil.copy("../micro-fans-notch-sphere/sphere32.h5", "./sphere32.h5")

    proc = subprocess.run(
        "mpiexec -n 8 micro-manager-precice --snapshot snapshot_config.json", shell=True
    )
    proc.check_returncode()

    os.remove("./input.json")
    os.remove("./sphere32.h5")
    os.remove("./PyFANS.so")
    os.remove("./snapshot_config.json")

    shutil.move("./snapshot_data.hdf5", "./snapshot_output.h5")


def generate_model():
    with h5py.File("./snapshot_output.h5", "r") as f:
        strains1to3 = f["strains1to3"][:]
        strains4to6 = f["strains4to6"][:]
        stresses1to3 = f["stresses1to3"][:]
        stresses4to6 = f["stresses4to6"][:]
        cmats = {k: f[k][:] for k in [f"cmat{i}" for i in range(1, 8)]}
    strains = np.zeros((len(strains1to3), 6))
    stresses = np.zeros((len(stresses1to3), 6))
    strains[:, 0:3] = strains1to3
    strains[:, 3:6] = strains4to6
    stresses[:, 0:3] = stresses1to3
    stresses[:, 3:6] = stresses4to6
    data_y = {"stresses": stresses}
    data_y.update(cmats)
    data_x = strains

    def split_samples(X, y, n_valid):
        n_samples = X.shape[0]
        if n_valid >= n_samples:
            raise AttributeError("The set number of validation points is invalid.")

        # Random split
        n_train = n_samples - n_valid
        choice = np.random.choice(range(n_samples), size=(n_train,), replace=False)
        ind = np.zeros(n_samples, dtype=bool)
        ind[choice] = True

        # Split samples
        X_train = X[ind]
        X_valid = X[~ind]

        # Split outputs
        y_train = {}
        y_valid = {}
        for key in y:
            if key == "x_values":
                continue
            y_train[key] = y[key][ind]
            y_valid[key] = y[key][~ind]

        return X_train, y_train, X_valid, y_valid

    data_t_x, data_t_y, data_v_x, data_v_y = split_samples(
        data_x, data_y, int(len(data_x) * 0.2)
    )

    if not os.path.exists("./model.py"):
        with open("./model.py", "w") as f:
            f.writelines(["def model(samples):\n", "\treturn None"])

    with open("./distributions.json", "r") as f:
        dist = json.load(f)
        std = dist["std"]
        mean = dist["mean"]

    model = PyLinkForwardModel()
    model.py_file = "model"
    model.name = "pyfans-surrogate"
    model.link_type = "function"
    model.output_names = ["stresses"]
    model.output_names.extend(list(cmats.keys()))

    inputs = Input()
    for i in range(6):
        inputs.add_marginals()
        inputs.marginals[i].name = f"$X_{i+1}"
        inputs.marginals[i].dist_type = "norm"
        inputs.marginals[i].parameters = [mean[i], std[i]]
    exp_design = ExpDesigns(inputs)
    exp_design.x = data_t_x
    exp_design.y = data_t_y

    meta_model = PCE(inputs)
    meta_model.meta_model_type = "aPCE"
    meta_model.pce_reg_method = "FastARD"
    meta_model.pce_deg = 5

    engine = Engine(meta_model, model, exp_design)
    engine.train_normal()

    with open(f"{model.name}.pkl", "wb") as f:
        joblib.dump(engine, f, 2)

    return data_v_x, data_v_y


def validate_surrogate(x_valid, y_valid, model_name="pyfans-surrogate.pkl"):
    with open(model_name, "rb") as input:
        engine = joblib.load(input)

    y_metamod, _ = engine.eval_metamodel(x_valid)
    # Compare predictions with true values
    def plot_line(y_true, y_pred):
        min = np.minimum(np.min(y_true), np.min(y_pred))
        max = np.maximum(np.max(y_true), np.max(y_pred))
        if min != 0:
            order = np.log10(min)
            diff_order = np.log10(np.abs(max - min))
            if diff_order < order:
                min -= np.power(10, order)
                max += np.power(10, order)

        plt.plot([min, max], [min, max], c="#000000", label="ideal")

    plt.figure()
    plt.scatter(
        y_valid["stresses"][:, 0], y_metamod["stresses"][:, 0], c="#FF0000", label="s1"
    )
    plt.scatter(
        y_valid["stresses"][:, 1], y_metamod["stresses"][:, 1], c="#FF3A00", label="s2"
    )
    plt.scatter(
        y_valid["stresses"][:, 2], y_metamod["stresses"][:, 2], c="#FF5500", label="s3"
    )
    plt.scatter(
        y_valid["stresses"][:, 3], y_metamod["stresses"][:, 3], c="#FF7700", label="s4"
    )
    plt.scatter(
        y_valid["stresses"][:, 4], y_metamod["stresses"][:, 4], c="#FF9900", label="s5"
    )
    plt.scatter(
        y_valid["stresses"][:, 5], y_metamod["stresses"][:, 5], c="#FFB300", label="s6"
    )
    plot_line(y_valid["stresses"], y_metamod["stresses"])
    plt.legend()
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: stresses")
    plt.show()

    colors = {
        "cmat1": ["#0000FF", "#0A00F5", "#1400EB"],
        "cmat2": ["#1E00E1", "#2800D7", "#3200CD"],
        "cmat3": ["#3C00C3", "#4600B9", "#5000AF"],
        "cmat4": ["#5A00A5", "#64009B", "#6E0091"],
        "cmat5": ["#780087", "#82007D", "#8C0073"],
        "cmat6": ["#960069", "#A0005F", "#AA0055"],
        "cmat7": ["#B4004B", "#BE0041", "#C80037"],
    }

    for key in colors.keys():
        for i in range(3):
            plt.figure()
            plt.scatter(
                y_valid[key][:, i],
                y_metamod[key][:, i],
                c=colors[key][i],
                label=f"{key}-{i}",
            )
            plot_line(y_valid[key][:, i], y_metamod[key][:, i])
            plt.legend()
            plt.xlabel("True Values")
            plt.ylabel("Predictions")
            plt.title(f"Validation: {key}-{i}")
            plt.show()


if not os.path.exists("./snapshot_input.h5"):
    meso_strain = run_meso()
    generate_snapshot_data(meso_strain, 100)

if not os.path.exists("./snapshot_output.h5"):
    run_snapshot()

if os.path.exists("./snapshot_output.h5"):
    val_x, val_y = generate_model()
    if os.path.exists("./pyfans-surrogate.pkl"):
        validate_surrogate(val_x, val_y)
