mod gwasss;

use std::env;

#[macro_use]
extern crate serde_derive;


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


fn parse_variantfile(args: &Vec<String>) -> Result<VariantFile, &'static str> {
    // Out of bound by default so we can know if they were not found.
    let mut variants_filename_idx = args.len();
    let mut variants_format_idx = args.len();

    for (i, arg) in args.iter().enumerate() {
        match arg.as_str() {
            "--variants-filename" => variants_filename_idx = i + 1,
            "--variants-format" => variants_format_idx = i + 1,
            _ => ()
        }
    }

    // Make sure both arguments were provided
    if variants_filename_idx >= args.len() {
        return Err("Provide --variants-filename");
    }

    else if variants_format_idx >= args.len() {
        return Err("Provide --variants-format");
    }

    // Match the format to an enum
    let format = match args[variants_format_idx].as_str() {
        "vcf" => VariantsFormat::VCF,
        "bim" => VariantsFormat::BIM,
        "stat" => VariantsFormat::STAT,
        _ => { return Err("Invalid --variants-format (use vcf, bim, stat)") }
    };

    Ok(VariantFile {
        variants_filename: args[variants_filename_idx].to_string(),
        variants_format: format
    })
}


fn cmd_variant_effect_matrix(datasets: Vec<gwasss::Dataset>, vf: VariantFile) {
    println!("{:#?}", datasets);
    println!("{:#?}", vf);

    // TODO
    // let variants: Iterator<Type=Variant> = vf.load_variants()
    // For all datasets, extract the variants into a Hash of:
    // Variant -> {"dataset:component": stat}
    // This will likely be a reusable stucture


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

    let mut args: Vec<String> = env::args().skip(1).collect();

    // Parse the root argument which is always required.
    let root_fn_idx = match args.iter().position(|x| x.as_str() == "--root") {
        Some(idx) => {
            // Check if we can parse the root.
            if idx + 1 >= args.len() {
                panic!("Need to provide a filename after --root");
            }
            idx + 1
        },
        None => panic!("Need to provide a --root")
    };

    let root_fn = &args[root_fn_idx].to_string();

    // Load the datasets from the root.
    let datasets = gwasss::load_datasets_from_manifests(
        gwasss::find_manifests(root_fn)
    );

    // Log.
    println!("Found and loaded {} datasets.", datasets.len());

    args.drain(..(root_fn_idx + 1));

    // Next up after --root is the command name.
    let command: Vec<String> = args.drain(..1).collect();
    let command = &command[0];

    match command.as_str() {
        "variant-effect-matrix" => {
            // variant-effect-matrix takes VariantFile
            let vf = parse_variantfile(&args).unwrap();
            cmd_variant_effect_matrix(datasets, vf);
        },
        _ => println!("Unknown command '{:?}'", command)
    }

}
