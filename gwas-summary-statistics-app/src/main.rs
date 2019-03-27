mod gwasss;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate enum_display_derive;

extern crate clap;
extern crate csv;

extern crate genepa_rs;

use std::error::Error;
use std::io::{Read};
use std::fs::File;
use std::collections::{HashMap, HashSet};
use clap::{Arg, App, SubCommand, AppSettings};
use crate::gwasss::{Dataset};
use genepa_rs::{Variant, OrderedAllelesVariant};


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

    // Format small p-values.
    let fmt_p = if stat.p < 0.05 {
        format!("{:e}", stat.p)
    } else { stat.p.to_string() };

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
        &fmt_p
    ]).expect("Could not write variant to output file.");

}


fn _init_writer<I, T>(path: &str, header: I) -> csv::Writer<std::fs::File>
    where I: IntoIterator<Item=T>, T: AsRef<[u8]>
{
    let mut writer = csv::WriterBuilder::new().from_path(path)
        .expect("Could not open file for writing");

    writer.write_record(header).expect("Could not write header");

    writer
}

fn _init_writer_default(path: &str) -> csv::Writer<std::fs::File> {
    _init_writer(path, &[
        "dataset_variant_name", "chrom", "pos",
        "reference_allele", "coded_allele",
        "dataset_name", "component_name", "population", "sex",
        "effect_type", "effect", "se", "p"
    ])
}


fn cmd_extract_region(datasets: Vec<Dataset>, args: &clap::ArgMatches) {
    let region = args.value_of("region").unwrap().to_string();
    let mut writer = _init_writer_default(args.value_of("output").unwrap());

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

    let mut writer = _init_writer_default(args.value_of("output").unwrap());

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


fn cmd_extract_variants(datasets: Vec<Dataset>, args: &clap::ArgMatches) {
    // Read the variants using the right reader.
    let format = args.value_of("variants_format").unwrap();
    let filename = args.value_of("variants_filename").unwrap();

    let variants: Vec<OrderedAllelesVariant> = if format == "vcf" {
        panic!("Can't read VCF files yet.");
    }
    else if format == "bim" {
        // Assume coded is A1 (minor allele).
        genepa_rs::plink::BimReader::new(filename).collect()
    }
    else if format == "stat" {
        // If we read the variants from a summary statistics file, we align
        // the allele with respect to the coded allele.
        gwasss::SummaryStatsFile::read_file(filename)
            .map(|stat| {
                let stat = stat.unwrap();

                match stat.coded_allele {
                    gwasss::CodedAllele::A1Coded => {
                        OrderedAllelesVariant {
                            variant: stat.variant, a1_idx: 0
                        }
                    },
                    gwasss::CodedAllele::A2Coded => {
                        OrderedAllelesVariant {
                            variant: stat.variant, a1_idx: 1
                        }
                    }
                }
            }).collect()
    }
    else {
        panic!(format!("Unknown format '{}'", format));
    };

    // Get an output writer.
    let mut writer = _init_writer_default(args.value_of("output").unwrap());

    for dataset in datasets.iter() {
        for component in dataset.components.iter() {
            // Iterate over ordered allele variants.
            // The "coded" allele is always 'A1'.
            for oav in variants.iter() {
                let coded_allele = if oav.a1_idx == 0 {
                    &oav.variant.alleles.0
                }
                else if oav.a1_idx == 1 {
                    &oav.variant.alleles.1
                }
                else {
                    panic!("Bad alleles.");
                };

                match component.get_stats_for_variant(
                    &oav.variant, coded_allele
                ) {
                    Ok(ref mut stat) => {
                        _write_csv_row(&mut writer, &dataset, &component,
                                       stat);
                    },
                    Err(e) => println!("WARN: {} - {:?}", &oav.variant, e)
                }
            }
        }
    }

}


#[derive(Debug, Deserialize)]
struct _ComponentFilter {
    trait_name: String,
    population: Option<gwasss::Population>,
    sex: Option<gwasss::Sex>,
}


impl _ComponentFilter {
    fn matches(&self, component: &gwasss::Component) -> bool {

        if self.trait_name != component.trait_name {
            return false;
        }

        if let Some(pop) = &self.population {
            if pop != &component.population {
                return false;
            }
        }

        if let Some(sex) = &self.sex {
            if sex != &component.sex {
                return false;
            }
        }

        true
    }
}


type _DatasetFilter = HashMap<String, Vec<_ComponentFilter>>;
fn load_filters(component_set: &str)
    -> Result<_DatasetFilter, serde_yaml::Error>
{
    let mut f = File::open(component_set).expect("Unable to open set file.");
    let mut contents = String::new();
    f.read_to_string(&mut contents).expect("Unable to read set file.");

    let dataset_filter: _DatasetFilter = serde_yaml::from_str(&contents)?;

    Ok(dataset_filter)
}


fn filter_datasets(set_file: &str, datasets: Vec<Dataset>)
    -> Vec<Dataset>
{
    // Load filters from yaml.
    let filters = load_filters(set_file).unwrap();

    // Datasets to keep (ones that are in set file)
    let keep_datasets: HashSet<&String> = filters.keys().collect();

    // New vector to hold the filtered datasets.
    let mut clean_datasets: Vec<Dataset> = Vec::new();

    for mut dataset in datasets {
        if keep_datasets.contains(&dataset.name) {
            // The dataset was requested apply filter to components.
            let component_filters = filters.get(&dataset.name).unwrap();

            // If there is filtering to be done on components.
            if component_filters.len() > 0 {
                let mut new_components: Vec<gwasss::Component> = Vec::new();

                for ref component in dataset.components {
                    if component_filters.iter().any(|f| f.matches(component)) {
                        new_components.push(component.to_owned());
                    }
                }

                dataset.components = new_components;
            }

            clean_datasets.push(dataset);
        }
    }

    clean_datasets
}


fn main() {
    // Here is the CLI I would like:
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
        .arg(Arg::with_name("trait_set")
            .long("trait-set")
            .value_name("trait_set")
            .help("YAML file containing sets of traits to extract.")
            .takes_value(true))

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

        .subcommand(SubCommand::with_name("extract-variants")
            .about("Extract a discrete set of variants represented in a file \
                    fith the 'summary statistics' format.")
            .arg(Arg::with_name("variants_filename")
                .long("variants-filename")
                .short("i")
                .help("Filename containing the input variants.")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("variants_format")
                .long("variants-format")
                .short("f")
                .help("File format containing the variants.")
                .takes_value(true)
                .possible_values(&["stat", "vcf", "bim"])
                .required(true))
            .arg(Arg::with_name("output")
                .long("output")
                .short("o")
                .help("Output filename (csv format)")
                .takes_value(true)
                .default_value("extracted_variants.csv"))
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
                .long("reference-allele")
                .short("ref")
                .help("Non-coded allele of the variant")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("coded_allele")
                .long("coded-allele")
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
    let mut datasets = gwasss::load_datasets_from_manifests(
        gwasss::find_manifests(matches.value_of("root").unwrap())
    );

    // Filter if need with respect to the trait set.
    if let Some(set_file) = matches.value_of("trait_set") {
        datasets = filter_datasets(set_file, datasets);
    }

    match matches.subcommand_name() {
        Some("extract-variant") => {
            cmd_extract_variant(
                datasets,
                matches.subcommand_matches("extract-variant").unwrap()
            );
        },
        Some("extract-variants") => {
            cmd_extract_variants(
                datasets,
                matches.subcommand_matches("extract-variants").unwrap()
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
