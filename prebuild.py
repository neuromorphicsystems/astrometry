import pathlib
import subprocess

dirname = pathlib.Path(__file__).resolve().parent

git_revision = (
    subprocess.run(
        ["git", "describe"],
        capture_output=True,
        check=True,
        cwd=dirname / "astrometry.net",
    )
    .stdout.decode()
    .strip()
)

git_date = (
    subprocess.run(
        ["git", "log", "-n", "1", "--format=%cd"],
        capture_output=True,
        check=True,
        cwd=dirname / "astrometry.net",
    )
    .stdout.decode()
    .strip()
    .replace(" ", "_")
)

with open(
    dirname / "astrometry.net" / "include" / "astrometry" / "os-features-config.h",
    "w",
) as output:
    output.write("#define HAVE_NETPBM 0\n")
    output.write(f'#define AN_GIT_REVISION "{git_revision}"\n')
    output.write(f'#define AN_GIT_DATE "{git_date}"\n')
    output.write('#define AN_GIT_URL "https://github.com/dstndstn/astrometry.net"\n')
