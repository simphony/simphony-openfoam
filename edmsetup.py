import sys
import click
import os
import subprocess

from packageinfo import BUILD, VERSION, NAME

if "WM_PROJECT" not in os.environ:
    print("To run this command you must source edmenv.sh first")
    sys.exit(1)

# The version of the buildcommon to checkout.
BUILDCOMMONS_VERSION="v0.2"


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


@click.group()
def cli():
    pass


@cli.command()
def egg():
    common.local_repo_to_edm_egg(".", name=NAME, version=VERSION, build=BUILD)
    with common.cd("openfoam-interface/internal-interface/wrapper"):
        common.run("python edmsetup.py egg")



@cli.command()
def upload_egg():
    egg_path = "endist/{NAME}-{VERSION}-{BUILD}.egg".format(
        NAME=NAME,
        VERSION=VERSION,
        BUILD=BUILD)
    click.echo("Uploading {} to EDM repo".format(egg_path))
    common.upload_egg(egg_path)

    with common.cd("openfoam-interface/internal-interface/wrapper"):
        try:
            common.run("python edmsetup.py upload_egg")
        except subprocess.CalledProcessError as e:
            print("Error during egg upload of submodule: {}. Continuing.".format(e))

    click.echo("Done")


@cli.command()
def clean():
    click.echo("Cleaning")
    common.clean(["endist", ".devenv"])


cli()
