import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from probeinterface import ProbeGroup
from probeinterface.generator import generate_multi_columns_probe
from probeinterface.io import write_probeinterface, read_probeinterface
from probeinterface.plotting import plot_probegroup

def make_probe(
    probes_path,
    probe_name,
    num_columns=1,
    num_channels=24,
    ypitch=150,
    yshift=[300],
    contact_radius=7.5,
    show_fig=False
):
    """
    Generate a probe layout, save it as JSON, and export a figure preview.

    Parameters
    ----------
    probes_path : Path or str
        Directory where the probe files will be saved.
    probe_name : str
        Name of the probe (used for filenames).
    num_columns : int, optional
        Number of columns in the probe (default = 1).
    num_channels : int, optional
        Total number of recording sites (default = 24).
    ypitch : float, optional
        Vertical distance between adjacent contacts (µm) (default = 150).
    yshift : list of float, optional
        Vertical offset per column (µm) (default = [300]).
    contact_radius : float, optional
        Radius of each contact (µm) (default = 7.5).
    """

    probes_path = Path(probes_path)
    probes_path.mkdir(parents=True, exist_ok=True)

    # Generate the probe
    probe = generate_multi_columns_probe(
        num_columns=num_columns,
        num_contact_per_column=int(num_channels / num_columns),
        ypitch=ypitch,
        y_shift_per_column=yshift,
        contact_shapes='circle',
        contact_shape_params={'radius': contact_radius}
    )

    # Add probe to group and write JSON
    probegroup = ProbeGroup()
    probegroup.add_probe(probe)
    json_path = probes_path / f"{probe_name}.json"
    write_probeinterface(json_path, probegroup)

    # Plot and save figure
    fig, ax = plt.subplots(figsize=(25, 10))
    plt.tight_layout()
    plot_probegroup(probegroup, with_contact_id=True, same_axes=True, ax=ax)
    ax.set_title(probe_name)
    fig.savefig(probes_path / f"{probe_name}.png", dpi=300, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)

    print(f"Saved probe JSON to {json_path}")

    return probegroup

def combine_probes(session_path,
                  probes_path,
                  show_fig=False):
    """
    Generate a combined probe layout for a given session, save it as JSON, and export a figure preview.

    Parameters
    ----------
    session_path : Path or str
        Directory where the metadata is contained and where probegroup will be saved.
    probes_path : Path or str
        Directory where individual probe files are saved.
    """
    
    session_path = Path(session_path)
    with open(session_path / "metadata.json", "r") as f:
        metadata = json.load(f)
    
    # Extract probes collected via Ripple system
    sess_name = metadata["sess_name"]
    probe_configs = metadata["probe_config"]
    hardware_configs = metadata["hardware_config"]
    
    probe_names = [
        p for p, h in zip(probe_configs, hardware_configs)
        if "elec" in h]
    
    # Read individual probe groups
    probegroups = []
    for name in probe_names:
        probe_path = Path(probes_path) / f"{name}.json"
        if not probe_path.exists():
            print(f"Probe map '{name}.json' does not exist in {probes_path}. It needs to be generated first.")
            return
        probegroups.append(read_probeinterface(probe_path))

    # If none match, exit the function
    if len(probe_names) == 0:
        return
    
    # Combine into a single ProbeGroup
    combined_pg = ProbeGroup()
    contact_num = 0
    for i, pg in enumerate(probegroups):
        for probe in pg.probes:
            this_probe = probe.copy()
            this_probe.move([i*500, 0])
            num_contacts = this_probe.get_contact_count()
            combined_pg.add_probe(this_probe)
            this_probe.set_device_channel_indices(
                np.arange(num_contacts) + contact_num)
            contact_num = contact_num + num_contacts
    
    print(combined_pg.get_global_device_channel_indices())
    
    json_path = session_path / f"{sess_name}_prbMap.json"
    write_probeinterface(json_path, combined_pg)
    
    # Plot and save figure
    fig, ax = plt.subplots(figsize=(25, 10))
    plt.tight_layout()
    plot_probegroup(combined_pg, with_contact_id=True, with_device_index=True, same_axes=True, ax=ax)
    ax.set_title(metadata["sess_name"])
    fig.savefig(session_path / f"{sess_name}_prbMap.png", dpi=300, bbox_inches="tight")
    if show_fig:
        plt.show()
    else:
        plt.close(fig)
    
    print(f"Saved probe JSON to {json_path}")

