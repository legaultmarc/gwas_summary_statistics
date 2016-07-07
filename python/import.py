#!/usr/bin/env python


from __future__ import print_function

import logging
import argparse
from getpass import getpass

import pyfaidx
import psycopg2

from settings import DATABASE as DBINFO


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("import")


# The complement
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def main():
    """The main."""
    args = parse_args()

    # Reading the reference
    human_ref = pyfaidx.Fasta(args.reference)

    # Creating the connection and getting the metadata
    conn, cur = db_connect()

    # If array
    if args.command == "array":
        import_array_name(conn, cur, args)
        import_array_variants(conn, cur, human_ref, args)


def import_array_name(db_conn, db_cur, args):
    """Adds the array to the DB (if absent)."""
    # Checking if the array name is already in the DB
    db_cur.execute(
        "SELECT COUNT(name) "
        "FROM gwas_results.GenotypingArray "
        "WHERE name = '{}'".format(args.name),
    )
    if db_cur.fetchone()[0] > 0:
        logger.info("Already an array named {} in the "
                    "database".format(args.name))
        return

    # Inserting the array into the DB
    db_cur.execute(
        "INSERT INTO gwas_results.GenotypingArray (name) VALUES (%s)",
        (args.name, ),
    )
    db_conn.commit()


def import_array_variants(db_conn, db_cur, human_ref, args):
    """Parses the array's annotation file and add to the DB."""
    # We read the data from the annotation file
    with open(args.annotation, "r") as f:
        # Reading until the header
        line = next(f)
        while not line.startswith("[Assay]"):
            line = next(f)

        # Reading the header
        header = {
            name: i
            for i, name in enumerate(next(f).rstrip("\r\n").split(args.sep))
        }

        # Reading the rest of the file
        chromosomes = []
        positions = []
        references = []
        alternatives = []
        logger.info("Inserting into database")
        for i, line in enumerate(f):
            if line.startswith("[Controls]"):
                break

            row = line.rstrip("\r\n").split(args.sep)

            # Getting the required information
            chrom = encode_chromosome(row[header[args.chrom]])
            pos = int(row[header[args.pos]])
            genotypes = row[header[args.geno]]

            # I've seen N/A in the past...
            if genotypes == "[N/A]":
                continue

            # Skipping if 0/0
            if chrom == "0" or pos == 0:
                continue

            # Separating the alleles
            a1 = genotypes.split("/")[0][1:]
            a2 = genotypes.split("/")[1][:-1]

            # Getting the reference allele
            human_ref_chrom = "chr" + chrom
            if human_ref_chrom == "chrXY":
                human_ref_chrom = "chrX"
            ref = str(human_ref[human_ref_chrom][pos-1]).upper()

            # Finding the alternative allele
            alt = None
            if a1 == "D" or a2 == "D":
                # This is an indel
                alt = "-"
            elif a1 == ref:
                alt = a2
            elif _COMP[a1] == ref:
                alt = _COMP[a2]
            elif a2 == ref:
                alt = a1
            elif _COMP[a2] == ref:
                alt = _COMP[a1]

            chromosomes.append(chrom)
            positions.append(pos)
            references.append(ref)
            alternatives.append(alt)

            if len(chromosomes) >= 100000:
                # Inserting
                db_cur.callproc(
                    "gwas_results.InsertVariantArray",
                    [chromosomes, positions, references, alternatives,
                     args.name],
                )
                db_conn.commit()

                # Resetting
                chromosomes = []
                positions = []
                references = []
                alternatives = []
                logger.info("Processed {:,d} markers".format(i+1))

    # Checking if we have more to process
    if len(chromosomes) > 0:
        # Inserting
        db_cur.callproc(
            "gwas_results.InsertVariantArray",
            [chromosomes, positions, references, alternatives, args.name],
        )
        db_conn.commit()


def encode_chromosome(chrom):
    """Encodes the chromosomes."""
    if (chrom.upper() == "MT") or (chrom == "26"):
        return "M"
    if chrom == "23":
        return "X"
    if chrom == "24":
        return "Y"
    if (chrom.upper()) == "YX" or (chrom == "25"):
        return "XY"
    return chrom.upper()


def db_connect():
    """Performs database connection using a database settings."""
    # Asking for the password
    pw = getpass("Password for {}: ".format(DBINFO["username"]))

    # The connection
    conn = psycopg2.connect(
        "dbname='{database}' user='{username}' host='{hostname}' "
        "password='{password}'".format(password=pw, **DBINFO)
    )
    cur = conn.cursor()

    return conn, cur


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(description="Imports data to PSQL")

    parent_parser = argparse.ArgumentParser(add_help=False)

    subparser = parser.add_subparsers(
        dest="command",
        help=("What kind of information to import")
    )
    subparser.required = True

    # The chip import
    parse_parser = subparser.add_parser(
        "array",
        help="Import an array information.",
        parents=[parent_parser]
    )

    group = parse_parser.add_argument_group("Input options")
    group.add_argument(
        "--name", type=str, metavar="NAME", required=True,
        help="The name of the array.",
    )
    group.add_argument(
        "--annotation", type=str, metavar="FILE", required=True,
        help="The annotation provided by the chip.",
    )
    group.add_argument(
        "--sep", type=str, metavar="SEP", default=",",
        help="The field separator of the annotation file. [%(default)s]",
    )
    group.add_argument(
        "--chrom", type=str, metavar="FIELD", default="chr",
        help="The field containing the chromosome. [%(default)s]",
    )
    group.add_argument(
        "--pos", type=str, metavar="FIELD", default="pos",
        help="The field containing the position. [%(default)s]",
    )
    group.add_argument(
        "--geno", type=str, metavar="FIELD", default="geno",
        help="The field containing the genotype. [%(default)s]",
    )
    group.add_argument(
        "--reference", type=str, metavar="FILE", required=True,
        help="The human reference.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
