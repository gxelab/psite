import click

from train import train
from predict import predict
from coverage import coverage
from pbam import pbam

@click.group()
def psite():
    """main interface"""
    pass

psite.add_command(train)
psite.add_command(predict)
psite.add_command(coverage)
psite.add_command(pbam)

if __name__ == '__main__':
    psite()
