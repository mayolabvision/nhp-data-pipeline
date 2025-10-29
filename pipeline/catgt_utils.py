import subprocess
from pathlib import Path
from config import CATGT_PATH

def run_catgt(session: str, raw_data_root: Path):
    """
    Run CatGT if it hasn't been run yet for the given session.
    """
    run_name = session.removesuffix("_g0")
    runit_path = Path(CATGT_PATH) / "runit.sh"
    catgt_dir = raw_data_root / f"{session}" / f"catgt_{session}"

    if catgt_dir.is_dir():
        return

    dest_folder = raw_data_root / session
    cmd = [
        str(runit_path),
        f"-dir={raw_data_root}",
        f"-run={run_name}",
        "-g=0",
        "-t=0,0",
        "-t_miss_ok",
        "-ni",
        "-prb=0",
        "-bf=0,0,-1,0,9,1",
        f"-dest={dest_folder}"
    ]

    print("Running CatGT command:")
    print(" ".join(map(str, cmd)))

    # capture output so we can show helpful diagnostics on failure
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        print("CatGT stdout:")
        print(result.stdout)
        print("CatGT stderr:")
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)

    # optional: verify that catgt_dir now exists after successful run
    if not catgt_dir.is_dir():
        raise RuntimeError(f"CatGT finished without creating expected output folder: {catgt_dir}")

