import sys
import click
import os
import shutil
import subprocess

from packageinfo import BUILD, VERSION, NAME

if "WM_PROJECT" not in os.environ:
    print("To run this command you must source edmenv.sh first")
    sys.exit(1)

# The version of the buildcommon to checkout.
BUILDCOMMONS_VERSION = "v0.2"


def bootstrap_devenv():
    try:
        os.makedirs(".devenv")
    except OSError:
        pass

    if not os.path.exists(".devenv/buildrecipes-common"):
        subprocess.check_call([
            "git", "clone", "-b", BUILDCOMMONS_VERSION,
            "http://github.com/simphony/buildrecipes-common.git",
            ".devenv/buildrecipes-common"
            ])
    sys.path.insert(0, ".devenv/buildrecipes-common")


bootstrap_devenv()
import buildcommons as common  # noqa


workspace = common.workspace()
common.edmenv_setup()


def write_endist_dat():
    with open("endist.dat", "w") as f:
        f.write("""
packages = [
  # A list of strings in the format
  # "dependency_name dependency_version_at_patchlevel_detail"
]
add_files=[
    ("{}", ".*", "EGG-INFO/usr/bin/"),
    ("{}", ".*", "EGG-INFO/usr/lib/"),
    (".", "post_egginst.py", "EGG-INFO/"),
]
""".format(os.environ["FOAM_USER_APPBIN"], os.environ["FOAM_USER_LIBBIN"]))


@click.group()
def cli():
    pass


@cli.command()
def egg():
    write_endist_dat()
    common.edmenv_run("python setup.py bdist_egg")

    try:
        os.makedirs("endist")
    except OSError:
        pass

    common.run(
        "edm repack-egg -b {build} -a rh6-x86_64 "
        "dist/{name}-{version}-py2.7-linux-x86_64.egg".format(
            name=NAME, version=VERSION, build=BUILD))

    edm_egg_filename = "{name}-{version}-{build}.egg".format(
        name=NAME, version=VERSION, build=BUILD)

    shutil.move(
        os.path.join("dist", edm_egg_filename),
        os.path.join("endist", edm_egg_filename)
        )

    return os.path.join("endist", edm_egg_filename)


@cli.command()
def upload_egg():
    egg_path = "endist/{NAME}-{VERSION}-{BUILD}.egg".format(
        NAME=NAME,
        VERSION=VERSION,
        BUILD=BUILD)
    click.echo("Uploading {} to EDM repo".format(egg_path))
    common.upload_egg(egg_path)
    click.echo("Done")


@cli.command()
def clean():
    click.echo("Cleaning")
    common.clean(["endist", ".devenv"])


cli()
