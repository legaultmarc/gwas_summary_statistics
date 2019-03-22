extern crate walkdir;

extern crate serde_derive;
extern crate serde_yaml;
extern crate genepa_rs;

use std::fmt::Display;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Read};
use std::fs::File;
use std::result::Result;
use std::iter::{Iterator, FromIterator};
use std::boxed::Box;
use std::path::Path;

use walkdir::WalkDir;

use genepa_rs::Variant;


#[derive(Debug, Display, Deserialize)]
pub enum EffectType {
    Beta,
    OR,
    HR,
    Other
}


#[derive(Debug, Display, Deserialize)]
pub enum Sex {
    Male,
    Female,
    Both
}

impl Default for Sex { fn default() -> Self { Sex::Both } }


#[derive(Debug, Display, Deserialize)]
pub enum Population {
    EUR,
    AIS,
    AFR,
    TRANS
}

impl Default for Population { fn default() -> Self { Population::EUR } }


#[derive(Debug, Display, Deserialize)]
pub enum CodedAllele {
    A1Coded,
    A2Coded
}

#[derive(Debug, Deserialize)]
pub struct Dataset {
    pub name: String,
    pub description: String,
    pub pmid: Option<u32>,
    pub url: Option<String>,
    pub components: Vec<Component>
}


#[derive(Debug, Deserialize)]
pub struct Component {
    pub trait_name: String,
    pub raw_url: Option<String>,
    pub formatted_file: String,

    #[serde(default)]
    pub population: Population,

    #[serde(default)]
    pub sex: Sex,

    pub effect_type: EffectType,
    pub n_cases: Option<u32>,
    pub n_controls: Option<u32>,
    pub n: Option<u32>,
}


#[allow(dead_code)]
impl Component {
    // Get association statistics for a single variant.
    // the coded_allele argument makes sure that the results are flipped
    // as required.
    pub fn get_stats_for_variant(&self, v: &Variant, coded_allele: &str)
        -> Result<AssociationStat, &'static str> {

        // Check that the required coded allele is an allelic form of the
        // variant.
        if coded_allele != v.alleles.0 && coded_allele != v.alleles.1 {
            return Err("Provided coded allele is not an allele of the \
                        variant.");
        }

        let region = format!(
            "{}:{}-{}", v.chrom.name, v.position, v.position
        );

        let tabix = SummaryStatsFile::tabix(&self.formatted_file, &region);

        // Convert the tabix results to a vector and filter irrelevant entries.
        let mut vec: Vec<AssociationStat> = tabix
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

        if vec.len() == 0 {
            return Err("Could not find variant in statistics file.");
        }

        else if vec.len() == 1 {
            let mut stat = vec.pop().unwrap();

            // Express the stats according to coded_allele.
            let current_coded_allele = match stat.coded_allele {
                CodedAllele::A1Coded => &stat.variant.alleles.0,
                CodedAllele::A2Coded => &stat.variant.alleles.1
            };

            if current_coded_allele != coded_allele {
                stat.flip_coded_allele(&self.effect_type);
            }

            return Ok(stat);
        }

        else {
            return Err("Variant has multiple entries in the summary \
                        statistics file");
        }
    }

    pub fn get_stats_for_region(&self, region: &str) -> SummaryStatsFile {
        SummaryStatsFile::tabix(&self.formatted_file, region)
    }
}


pub struct SummaryStatsFile {
    iter: Box<Iterator<Item=std::io::Result<String>>>
}


impl SummaryStatsFile {
    pub fn tabix(filename: &str, region: &str) -> SummaryStatsFile {
        // Spawn the process for the iterator.
        let output = Command::new("tabix")
            .arg(filename)
            .arg(region)
            .stdout(Stdio::piped())
            .spawn()
            .expect("Tabix failed");

        SummaryStatsFile {
            iter: Box::new(BufReader::new(output.stdout.unwrap()).lines())
        }
    }

    pub fn read_file(filename: &str) -> SummaryStatsFile {
        let f = File::open(filename)
            .expect(&format!("Couldn't open file: {:?}", filename));

        let mut iter = BufReader::new(f).lines();

        // Because there is a header, we skip it and assume the columns
        // are defined as per the spec.
        iter.next();

        SummaryStatsFile { iter: Box::new(iter) }
    }
}


fn _split_line_to_variant(fields: &Vec<&str>) -> (Variant, CodedAllele) {
    let reference_allele = fields[3].to_string().to_uppercase();
    let coded_allele = fields[4].to_string().to_uppercase();

    let v = Variant::new(
        fields[0].to_string(),  // name
        fields[1].to_string(),  // chrom
        fields[2].to_string().parse().unwrap(),  // pos
        (reference_allele, coded_allele.clone())  // alleles
    );

    let code = if coded_allele == v.alleles.0 { CodedAllele::A1Coded }
               else { CodedAllele::A2Coded };

    (v, code)
}


impl Iterator for SummaryStatsFile {
    type Item = Result<AssociationStat, Box<str>>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.iter.next() {
            // Split the string.
            let line = s.unwrap();
            let str_vec = Vec::from_iter(line.split('\t'));

            let (v, code) = _split_line_to_variant(&str_vec);

            let assoc = AssociationStat {
                variant: v,
                coded_allele: code,
                effect: str_vec[5].parse().unwrap(),
                se: str_vec[6].parse().unwrap(),
                p: str_vec[7].parse().unwrap()
            };

            return Some(Ok(assoc));
        }

        None
    }
}


#[derive(Debug)]
pub struct AssociationStat {
    pub variant: Variant,
    pub coded_allele: CodedAllele,
    pub effect: f32,
    pub se: f32,
    pub p: f32
}


impl AssociationStat {
    pub fn get_reference_allele(&self) -> &str {
        match self.coded_allele {
            CodedAllele::A1Coded => &self.variant.alleles.1,
            CodedAllele::A2Coded => &self.variant.alleles.0,
        }
    }

    pub fn get_coded_allele(&self) -> &str {
        match self.coded_allele {
            CodedAllele::A1Coded => &self.variant.alleles.0,
            CodedAllele::A2Coded => &self.variant.alleles.1,
        }
    }

    pub fn flip_coded_allele(&mut self, effect_type: &EffectType) {
        match effect_type {
            EffectType::OR | EffectType::HR => {
                let beta = -self.effect.ln();
                self.effect = beta.exp();
            },
            EffectType::Beta => self.effect = -self.effect,
            EffectType::Other => panic!(
                "Could not flip coded allele for unknown effect type."
            )
        };

        match self.coded_allele {
            CodedAllele::A1Coded => self.coded_allele = CodedAllele::A2Coded,
            CodedAllele::A2Coded => self.coded_allele = CodedAllele::A1Coded,
        };
    }
}


pub fn find_manifests(root_dir: &str) -> Vec<String> {
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


pub fn load_datasets_from_manifests(manifests: Vec<String>) -> Vec<Dataset> {
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


fn _sanitize_dataset(dataset: &mut Dataset, meta_path: &Path) {
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


pub fn load_dataset_from_manifest(manifest: &str)
    -> Result<Dataset, serde_yaml::Error> {

    let mut f = File::open(manifest).expect("Unable to open the file");
    let mut contents = String::new();
    f.read_to_string(&mut contents).expect("Unable to read the file");

    let mut dataset: Dataset = serde_yaml::from_str(&contents)?;

    _sanitize_dataset(&mut dataset, &Path::new(manifest));

    Ok(dataset)
}