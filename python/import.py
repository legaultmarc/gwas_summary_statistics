#!/usr/bin/env python


from __future__ import print_function

import sys
import logging
import argparse
from getpass import getpass

import pyfaidx
import psycopg2


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

    # Creating the connection and getting the metadata
    conn, cur = db_connect(args.host_name, args.database_name, args.user_name)

    # If array
    if args.command == "array":
        # Reading the reference
        human_ref = pyfaidx.Fasta(args.reference)

        import_array_name(conn, cur, args)
        import_array_variants(conn, cur, human_ref, args)

    # If stydy
    if args.command == "study":
        import_study(conn, cur, args)


def import_study(db_conn, db_cur, args):
    """Adds a study to the DB."""
    try:
        db_cur.execute(
            "INSERT INTO gwas_results.Study (name, parent, genotyping_array) "
            "    VALUES (%s, %s, %s)",
            (args.name, args.parent_study, args.array),
        )

    except psycopg2.IntegrityError as exception:
        if "duplicate key value violates unique constraint" in str(exception):
            if "Key (name)=" in str(exception):
                logger.error("Study '{}' already exists".format(args.name))
                sys.exit(1)

        elif "violates foreign key constraint" in str(exception):
            if "Key (genotyping_array)=" in str(exception):
                logger.error("Genotyping array '{}' doesn't "
                             "exists".format(args.array))
                sys.exit(1)

            elif "Key (parent)=" in str(exception):
                logger.error("Parent study '{}' doesn't "
                             "exists".format(args.parent_study))
                sys.exit(1)

        raise

    db_conn.commit()


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


def db_connect(hostname, database, username):
    """Performs database connection using a database settings.

    Args:
        hostname (str): The host name of the database server.
        database (str): The name of the database.
        username (str): The name of the user.

    """
    # Asking for the password
    pw = getpass("Password for {}: ".format(username))

    # The connection
    conn = psycopg2.connect(
        "dbname='{database}' user='{username}' host='{hostname}' "
        "password='{password}'".format(
            database=database, username=username, hostname=hostname,
            password=pw,
        )
    )
    cur = conn.cursor()

    return conn, cur


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(description="Imports data to PSQL")

    parent_parser = argparse.ArgumentParser(add_help=False)

    # The database parameters
    group = parent_parser.add_argument_group("Database options")
    group.add_argument(
        "-H", "--host-name", type=str, metavar="HOST", required=True,
        help="The host name of the database server.",
    )
    group.add_argument(
        "-U", "--user-name", type=str, metavar="USER", required=True,
        help="The user name for the database.",
    )
    group.add_argument(
        "-D", "--database-name", type=str, metavar="NAME", required=True,
        help="The name of the database.",
    )

    subparser = parser.add_subparsers(
        dest="command",
        help="What kind of information to import",
    )
    subparser.required = True

    # The chip import
    import_array_parser = subparser.add_parser(
        "array",
        help="Import an array information.",
        parents=[parent_parser]
    )

    group = import_array_parser.add_argument_group("Genotyping array options")
    group.add_argument(
        "--name", type=str, metavar="NAME", required=True,
        help="The name of the array.",
    )

    group = import_array_parser.add_argument_group("Annotation file options")
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

    group = import_array_parser.add_argument_group("Reference options")
    group.add_argument(
        "--reference", type=str, metavar="FILE", required=True,
        help="The human reference fasta file.",
    )

    # The chip import
    import_study_parser = subparser.add_parser(
        "study",
        help="Import a study.",
        parents=[parent_parser]
    )

    group = import_study_parser.add_argument_group("Study options")
    group.add_argument(
        "--name", type=str, metavar="STUDY", required=True,
        help="The name of the study.",
    )
    group.add_argument(
        "--array", type=str, metavar="ARRAY", required=True,
        help="The name of the genotyping arry.",
    )
    group.add_argument(
        "--parent-study", type=str, metavar="STUDY",
        help="The name of the parent study (if any).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
