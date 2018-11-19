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
@click.argument('npzfile', metavar='<npzfile>', type=click.Path(exists=True, resolve_path=True, dir_okay=False))
@click.option('-r', '--common-scale', metavar='<common_scale>', type=float, default=10000,
help='Read counts per cell after scaliing (default: 10000)')
@click.option('-o', '--outfile', metavar='<outfile>', type=click.Path(dir_okay=False), default=None, help='Name of output file')
@click.option('-v', '--verbose', count=True, help='-v Level 1, -vv Level 2 verbosity')
def extra_zero_test(npzfile, common_scale, outfile, verbose):
    """
    Tests if excess zeros exist
    """
    utils.configure_logging(verbose)
    tenxt.extra_zero_test(npzfile, common_scale, outfile)


@main.command()
@click.argument('loomfile', metavar='<loomfile>', type=click.Path(exists=True, dir_okay=False))
@click.option('-c', '--chunk', metavar='<chunk_size>', type=int, default=25, help='Number of genes in each chunk')
@click.option('-o', '--outdir', metavar='<outdir>', type=click.Path(exists=True, resolve_path=True, file_okay=False), default='.',
help='Folder name to store parameter files')
@click.option('--email', metavar='<email>', type=str, default=None, help='Notification E-mail')
@click.option('--queue', metavar='<queue>', type=str, default=None, help='Queue name')
@click.option('--mem', metavar='<mem>', type=int, default=0, help='Memory in GB (default: 16GB)')
@click.option('--walltime', metavar='<walltime>', type=int, default=0, help='Walltime in hours (default: 24h)')
@click.option('--systype', metavar='<systype>', default='pbs', help='Type of HPC cluster system (default: pbs)')
@click.option('--dryrun', is_flag=True, help='Use this when you want to rehearse your submit commands')
@click.option('-v', '--verbose', count=True, help='\'-v\' is Level 1 and \'-vv\' is Level 2')
def submit(loomfile, chunk, outdir, email, queue, mem, walltime, systype, dryrun, verbose):
    """
    Submits scBASE fitting jobs to HPC clusters
    """
    utils.configure_logging(verbose)
    tenxt.submit(loomfile, chunk, outdir, email, queue, mem, walltime, systype, dryrun)


if __name__ == "__main__":
    sys.exit(main())
