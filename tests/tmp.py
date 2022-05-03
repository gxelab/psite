import click

@click.command()
@click.option('-v', '--verbose', is_flag=True)
def main(verbose):
    print(verbose)


if __name__ == '__main__':
    main()