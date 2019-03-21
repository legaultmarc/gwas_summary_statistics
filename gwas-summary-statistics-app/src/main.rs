mod gwasss;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate enum_display_derive;

extern crate clap;
extern crate csv;

extern crate genepa_rs;

use std::error::Error;
use clap::{Arg, ArgGroup, ArgMatches, App, SubCommand, AppSettings};
use crate::gwasss::{Dataset};
use genepa_rs::Variant;


#[derive(Debug)]
enum VariantsFormat {
    VCF,
    BIM,
    STAT
}


#[derive(Debug)]
struct VariantFile {
    variants_filename: String,
    variants_format: VariantsFormat
}


#[derive(Debug)]
struct ComponentsPair {
    x: String,
    y: String
}



fn _parse_variant_from_args(args: &clap::ArgMatches)
    -> Result<Variant, Box<Error>> {

    Ok(Variant::new(
        String::from(""),  // We don't care about the name and it's not optional yet.
        String::from(args.value_of("chrom").unwrap()),
        String::from(args.value_of("pos").unwrap()).parse()?,
        (String::from(args.value_of("ref_allele").unwrap()),
         String::from(args.value_of("coded_allele").unwrap())),
    ))
}


// Write a CSV entry for a statistics result on a given dataset and component.
fn _write_csv_row(writer: &mut csv::Writer<std::fs::File>,
                  dataset: &gwasss::Dataset,
                  component: &gwasss::Component,
                  stat: &gwasss::AssociationStat) {

    writer.write_record(&[
        &stat.variant.name,
        &stat.variant.chrom.name,
        &stat.variant.position.to_string(),
        stat.get_reference_allele(),
        stat.get_coded_allele(),
        &dataset.name,
        &component.trait_name,
        &component.population.to_string(),
        &component.sex.to_string(),
        &component.effect_type.to_string(),
        &stat.effect.to_string(),
        &stat.se.to_string(),
        &stat.p.to_string()
    ]).expect("Could not write variant to output file.");

}


fn cmd_extract_region(datasets: Vec<Dataset>, args: &clap::ArgMatches) {
    let region = args.value_of("region").unwrap().to_string();
    let output = args.value_of("output").unwrap();

    let mut writer = csv::WriterBuilder::new().from_path(output)
        .expect("Could not open file for writing");

    writer.write_record(&[
        "dataset_variant_name", "chrom", "pos",
        "reference_allele", "coded_allele",
        "dataset_name", "component_name", "population", "sex", "effect_type",
        "effect", "se", "p"
    ]).expect("Could not write header");

    for dataset in datasets.iter() {
        for component in dataset.components.iter() {
            for mut result in component.get_stats_for_region(&region) {
                match result {
                    Ok(ref mut stat) => {
                        _write_csv_row(&mut writer, &dataset, &component,
                                       stat);
                    },
                    Err(e) => println!("{} :: {:?}", e, component)
                }
            }
        }
    }

}


fn cmd_extract_variant(datasets: Vec<Dataset>, args: &clap::ArgMatches) {
    // Parse the variant from the parameters.
    let v = _parse_variant_from_args(args)
        .expect("Could not parse variant from command arguments.");

    let output = args.value_of("output").unwrap();

    let mut writer = csv::WriterBuilder::new().from_path(output)
        .expect("Could not open file for writing");

    // Header
    writer.write_record(&[
        "dataset_variant_name", "chrom", "pos",
        "reference_allele", "coded_allele",
        "dataset_name", "component_name", "population", "sex", "effect_type",
        "effect", "se", "p"
    ]).expect("Could not write header");

    // Extract variant if possible for every dataset.
    for dataset in datasets.iter() {
        for component in dataset.components.iter() {
            match component.get_stats_for_variant(
                &v, args.value_of("coded_allele").unwrap()
            ) {
                Ok(ref mut stat) => {
                    // The variant was found.
                    _write_csv_row(&mut writer, &dataset, &component, stat)
                },
                Err(e) => println!("{} :: {:?}", e, component)
            }
        }
    }

    writer.flush().expect("Broken flush");

}


fn main() {
    // Here is the CLI I would like:
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   variant-effect-matrix \
    //   --variants-filename my_file \
    //   --variants-format {bim, vcf, stats}
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   extract-variant \
    //   --chrom 3 \
    //   --pos 123456 \
    //   --reference_allele A \
    //   --coded_allele G
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   extract-region \
    //   --chrom 3 \
    //   --start 123456 \
    //   --end 127456
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   correlation \
    //   --x megastroke_2018:any_stroke \
    //   --y cardiogram:cad \
    //   --all-variants
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   correlation \
    //   --x megastroke_2018:any_stroke \
    //   --y cardiogram:cad \
    //   --chrom 3 \
    //   --start 123456 \
    //   --end 127456
    //
    // summary-stats \
    //   --root /data/project/summary_statistics \
    //   correlation \
    //   --x megastroke_2018:any_stroke \
    //   --y cardiogram:cad \
    //   --variants-filename my_file \
    //   --variants-format {bim, vcf, stats} \

    let matches = App::new("Gwas Summary Statistics CLI")
        .version("0.9")
        .author("Marc-André Legault")
        .about("Interface to query summary statistics from large GWAS \
                consortia.")

        .setting(AppSettings::SubcommandRequired)

        .arg(Arg::with_name("root")
            .long("root")
            .value_name("root_path")
            .help("Root directory from which manifests will be searched.")
            .takes_value(true)
            .required(true))

        .subcommand(SubCommand::with_name("extract-region")
            .about("Extract a genomic region.")
            .arg(Arg::with_name("region")
                .long("region")
                .help("Region of the form chrom:start-end (e.g. \
                      15:1234567-1253192")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("output")
                .long("output")
                .short("o")
                .help("Output filename (csv format)")
                .takes_value(true)
                .default_value("extracted_region.csv"))
        )

        .subcommand(SubCommand::with_name("extract-variant")
            .about("Extract the effect of a single variant on all components.")
            .arg(Arg::with_name("chrom")
                .long("chrom")
                .help("Chromosome of the variant")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("pos")
                .long("pos")
                .help("Position of the variant")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("ref_allele")
                .long("reference_allele")
                .short("ref")
                .help("Non-coded allele of the variant")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("coded_allele")
                .long("coded_allele")
                .short("coded")
                .help("Coded allele of the variant")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("output")
                .long("output")
                .short("o")
                .help("Output filename (csv format)")
                .takes_value(true)
                .default_value("extracted_variant.csv"))
        )

        .get_matches();

    // Find and parse datasets from the root.
    let datasets = gwasss::load_datasets_from_manifests(
        gwasss::find_manifests(matches.value_of("root").unwrap())
    );

    match matches.subcommand_name() {
        Some("extract-variant") => {
            cmd_extract_variant(
                datasets,
                matches.subcommand_matches("extract-variant").unwrap()
            );
        },
        Some("extract-region") => {
            cmd_extract_region(
                datasets,
                matches.subcommand_matches("extract-region").unwrap()
            );
        },
        Some(cmd) => println!("Command '{}' isn't supported yet.", cmd),
        None => panic!("No subcommand provided")
    }

}
