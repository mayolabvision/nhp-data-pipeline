import subprocess
from pathlib import Path

def run_catgt(session: str, raw_data_root: Path):
    """
    Run CatGT if it hasn't been run yet for the given session.
    """
    run_name = session.removesuffix("_g0")
    runit_path = Path("/ix1/pmayo/packages/CatGT-linux/runit.sh")
    catgt_dir = raw_data_root / f"{session}" / f"catgt_{session}"

    if not catgt_dir.is_dir():
        print(f"CatGT output not found for {session}. Running CatGT...")
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
            f"-dest={raw_data_root / session}"
        ]
        subprocess.run(cmd, check=True)
    else:
        print(f"CatGT output already exists for {session}, skipping.")
