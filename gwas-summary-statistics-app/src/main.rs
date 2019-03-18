extern crate statrs;
extern crate genepa_rs;

use std::process::{Command, Stdio};
use std::marker::PhantomData;
use std::io::{BufReader, BufRead, Result};
use std::iter::{Iterator, FromIterator};
use std::boxed::Box;

use genepa_rs::Variant;


#[derive(Debug)]
enum EffectType {
    Beta,
    OR,
    HR,
    Other
}


#[derive(Debug)]
enum Sex {
    Male,
    Female,
    Both
}


#[derive(Debug)]
enum Population {
    EUR,
    AIS,
    AFR,
    TRANS
}


#[derive(Debug)]
enum CodedAllele {
    A1Coded,
    A2Coded
}

#[derive(Debug)]
struct Dataset {
    name: String,
    description: String,
    pmid: Option<u32>,
    url: Option<String>,
    components: Vec<Component>
}


#[derive(Debug)]
struct Component {
    trait_name: String,
    raw_url: Option<String>,
    formatted_file: String,
    population: Population,
    sex: Sex,
    effect_type: EffectType,
    n_cases: Option<u32>,
    n_controls: Option<u32>,
    n: Option<u32>,
}


impl Component {
    // Get association statistics for a single variant.
    fn get_stats_for_variant(&self, v: Variant) -> AssociationStat {
        // TODO
        AssociationStat {
            variant: v,
            coded_allele: CodedAllele::A1Coded,
            effect: -1.2,
            se: 0.13
        }
    }

    // -> StatIterator
    fn get_stats_for_region(&self, region: String) {
        let tabix = TabixSummaryStats::new(&self.formatted_file, region);

        for stat in tabix {
            println!("{:?}", stat);
        }
        println!("DONE");
    }
}


struct TabixSummaryStats {
    pub filename: String,
    pub region: String,
    iter: Box<Iterator<Item=Result<String>>>
}


impl TabixSummaryStats {
    fn new(filename: &String, region: String) -> TabixSummaryStats {
        // Spawn the process for the iterator.
        let output = Command::new("tabix")
            .arg(filename)
            .arg(&region)
            .stdout(Stdio::piped())
            .spawn()
            .expect("Tabix failed");

        TabixSummaryStats {
            filename: filename.clone(), region: region,
            iter: Box::new(BufReader::new(output.stdout.unwrap()).lines())
        }
    }
}

impl Iterator for TabixSummaryStats {
    type Item = AssociationStat;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.iter.next() {
            // Split the string.
            let line = s.unwrap();

            let str_vec = Vec::from_iter(line.split('\t'));

            let v = Variant::new(
                str_vec[0].to_string(),  // name
                str_vec[1].to_string(),  // chrom
                str_vec[2].to_string().parse().unwrap(),  // pos
                (str_vec[3].to_string(), str_vec[4].to_string())  // alleles
            );

            let assoc = AssociationStat {
                variant: v,
                coded_allele: CodedAllele::A1Coded, // FIXME
                effect: str_vec[5].parse().unwrap(),
                se: str_vec[6].parse().unwrap(),
            };

            return Some(assoc);
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

impl AssociationStat {
    // Get an instance.
    fn tst() -> AssociationStat {
        let v = Variant::new(
            String::from("rs12345"),
            String::from("3"),
            12345,
            (String::from("G"), String::from("C"))
        );

        AssociationStat {
            variant: v,
            coded_allele: CodedAllele::A1Coded,
            effect: -0.1,
            se: 0.23,
        }
    }
}


fn find_manifests(root_dir: &str) -> Vec<String> {
    // TODO
    Vec::new()
}


fn load_datasets_from_manifests(manifests: Vec<String>) -> Vec<Dataset> {
    // TODO
    Vec::new()
}


fn main() {
    let v = Variant::new(
        String::from("rs12345"),
        String::from("3"),
        12345,
        (String::from("G"), String::from("C"))
    );

    let mut dataset = Dataset {
        name: String::from("cardiogram"),
        description: String::from("test"),
        pmid: None,
        url: None,
        components: Vec::new()
    };

    let cad = Component {
        trait_name: String::from("CAD"),
        raw_url: None,
        formatted_file: String::from("/Users/legaultmarc/projects/StatGen/HCN4/data/eppinga_study/_temp/meta_analysis_results.formatted.tsv.gz"),
        population: Population::EUR,
        sex: Sex::Both,
        effect_type: EffectType::OR,
        n: None,
        n_cases: None,
        n_controls: None,
    };

    let mi = Component {
        trait_name: String::from("mi"),
        raw_url: None,
        formatted_file: String::from("test_mi.txt"),
        population: Population::EUR,
        sex: Sex::Both,
        effect_type: EffectType::OR,
        n: None,
        n_cases: None,
        n_controls: None,
    };

    dataset.components.push(cad);
    dataset.components.push(mi);

    println!("{}", v);
    println!("{:?}", dataset);

    // Call the get_stats_for_variant on a component.
    dataset.components[0].get_stats_for_region(
        String::from("3:38621237-179172979")
    );

}
