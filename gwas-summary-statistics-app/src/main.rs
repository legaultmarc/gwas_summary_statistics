extern crate walkdir;
extern crate statrs;
#[macro_use]
extern crate serde_derive;
extern crate serde_yaml;


extern crate genepa_rs;

use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Read};
use std::fs::File;
use std::result::Result;
use std::iter::{Iterator, FromIterator};
use std::boxed::Box;
use std::path::Path;

use walkdir::WalkDir;

use genepa_rs::Variant;


#[derive(Debug, Deserialize)]
enum EffectType {
    Beta,
    OR,
    HR,
    Other
}


#[derive(Debug, Deserialize)]
enum Sex {
    Male,
    Female,
    Both
}

impl Default for Sex { fn default() -> Self { Sex::Both } }


#[derive(Debug, Deserialize)]
enum Population {
    EUR,
    AIS,
    AFR,
    TRANS
}

impl Default for Population { fn default() -> Self { Population::EUR } }


#[derive(Debug, Deserialize)]
enum CodedAllele {
    A1Coded,
    A2Coded
}

#[derive(Debug, Deserialize)]
struct Dataset {
    name: String,
    description: String,
    pmid: Option<u32>,
    url: Option<String>,
    components: Vec<Component>
}


#[derive(Debug, Deserialize)]
struct Component {
    trait_name: String,
    raw_url: Option<String>,
    formatted_file: String,

    #[serde(default)]
    population: Population,

    #[serde(default)]
    sex: Sex,

    effect_type: EffectType,
    n_cases: Option<u32>,
    n_controls: Option<u32>,
    n: Option<u32>,
}


#[allow(dead_code)]
impl Component {
    // Get association statistics for a single variant.
    fn get_stats_for_variant(&self, v: &Variant)
        -> Result<AssociationStat, &'static str> {

        let region = format!(
            "{}:{}-{}", v.chrom.name, v.position, v.position
        );

        let tabix = TabixSummaryStats::new(&self.formatted_file, &region);

        // Convert to vector, filter out errors and 
        let mut v: Vec<AssociationStat> = tabix
            .into_iter()
            .filter_map(|result| {
                // Keep only statistics matching the variant.
                match result {
                    Ok(stat) => {
                        if &stat.variant == v { Some(stat) } else { None }
                    },
                    Err(_) => None
                }
            })
            .collect();


        if v.len() == 0 {
            return Err("Could not find variant in statistics file.");
        }

        else if v.len() == 1 {
            return Ok(v.pop().unwrap());
        }

        else {
            return Err("Variant has multiple entries in the summary \
                        statistics file");
        }
    }

    fn get_stats_for_region(&self, region: &str) -> TabixSummaryStats{
        TabixSummaryStats::new(&self.formatted_file, region)
    }
}


struct TabixSummaryStats {
    iter: Box<Iterator<Item=std::io::Result<String>>>
}


impl TabixSummaryStats {
    fn new(filename: &str, region: &str) -> TabixSummaryStats {
        // Spawn the process for the iterator.
        let output = Command::new("tabix")
            .arg(filename)
            .arg(region)
            .stdout(Stdio::piped())
            .spawn()
            .expect("Tabix failed");

        TabixSummaryStats {
            iter: Box::new(BufReader::new(output.stdout.unwrap()).lines())
        }
    }
}

impl Iterator for TabixSummaryStats {
    type Item = Result<AssociationStat, &'static str>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.iter.next() {
            // Split the string.
            let line = s.unwrap();
            let str_vec = Vec::from_iter(line.split('\t'));

            let reference_allele = str_vec[3].to_string();
            let coded_allele = str_vec[4].to_string();

            let v = Variant::new(
                str_vec[0].to_string(),  // name
                str_vec[1].to_string(),  // chrom
                str_vec[2].to_string().parse().unwrap(),  // pos
                (reference_allele.clone(), coded_allele.clone())  // alleles
            );

            // Check if coded allele is first or second.
            let coded_allele_pos = if (reference_allele == v.alleles.0) &&
                                      (coded_allele == v.alleles.1) {
                CodedAllele::A2Coded
            }
            else if (reference_allele == v.alleles.1) &&
                    (coded_allele == v.alleles.0) {
                CodedAllele::A1Coded
            }
            else {
                return Some(Err("Bad allele parsing in summary stats file."));
            };

            let assoc = AssociationStat {
                variant: v,
                coded_allele: coded_allele_pos,
                effect: str_vec[5].parse().unwrap(),
                se: str_vec[6].parse().unwrap(),
            };

            return Some(Ok(assoc));
        }

        None
    }
}


#[derive(Debug)]
struct AssociationStat {
    variant: Variant,
    coded_allele: CodedAllele,
    effect: f32,
    se: f32
}


fn find_manifests(root_dir: &str) -> Vec<String> {
    WalkDir::new(root_dir)
        .into_iter()
        .filter_map(|x| {
            let path = x.as_ref().unwrap();

            if path.file_name() == "GWAS_MANIFEST.yaml" {
                Some(String::from(path.path().to_str().unwrap()))
            }
            else {
                None
            }
        })
        .collect()
}


fn load_datasets_from_manifests(manifests: Vec<String>) -> Vec<Dataset> {
    let mut v = Vec::new();

    for manifest in manifests {
        match load_dataset_from_manifest(&manifest) {
            Ok(dataset) => v.push(dataset),
            Err(e) => {
                println!("Ignoring malformed manifest '{}'", manifest);
                println!("Error: {:?}", e);
            },
        }
    }

    v
}


fn sanitize_dataset(dataset: &mut Dataset, meta_path: &Path) {
    // For now, the only thing we may want to do is to expand the formatted
    // file path for components and make sure they exist.
    //
    // Eventually, this could be used to inherite properties from the
    // dataset (for all components).
    for component in dataset.components.iter_mut() {
        if component.formatted_file.contains("${DATASET_ROOT}") {
            component.formatted_file = component.formatted_file.replace(
                "${DATASET_ROOT}",
                meta_path.parent().unwrap().to_str().unwrap()
            );
        }
    }
}


fn load_dataset_from_manifest(manifest: &str)
    -> Result<Dataset, serde_yaml::Error> {

    let mut f = File::open(manifest).expect("Unable to open the file");
    let mut contents = String::new();
    f.read_to_string(&mut contents).expect("Unable to read the file");

    let mut dataset: Dataset = serde_yaml::from_str(&contents)?;

    sanitize_dataset(&mut dataset, &Path::new(manifest));

    Ok(dataset)
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

    let datasets = load_datasets_from_manifests(
        find_manifests("/data/projects/summary_statistics/")
    );

    for dataset in datasets {
        println!("{:#?}", dataset);
    }
}
