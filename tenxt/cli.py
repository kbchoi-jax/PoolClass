# -*- coding: utf-8 -*-

"""Console script for tenxt."""
import sys
import click
from . import tenxt
from . import utils
from . import __logo__, __version__
from . import get_data


@click.group()
@click.version_option(version=__version__, message=__logo__)
def main(args=None):
    """Console script for tenxt"""
    return 0


@main.command()
@click.argument('cntfile', metavar='<cntfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('-v', '--verbose', count=True, help='-v Level 1, -vv Level 2 verbosity')
def score_test(cntfile, verbose):
    """
    Tests if excess zeros exist
    """
    utils.configure_logging(verbose)
    tenxt.score_test(cntfile)


if __name__ == "__main__":
    sys.exit(main())
