import click

from psite.train import train
from psite.predict import predict
from psite.coverage import coverage
from psite.pbam import pbam

@click.group()
@click.help_option("--help", "-h")
def psite():
    """main interface"""
    pass

psite.add_command(train)
psite.add_command(predict)
psite.add_command(coverage)
psite.add_command(pbam)

if __name__ == '__main__':
    psite()
