mod segtp;

use std::error::Error;
use clap::{Command, Arg};
use std::path::Path;

fn run() -> Result<(),Box<dyn Error>>{
    let mut s1 = segtp::SEGTP::new(11,10)?;
    let mut s2 = segtp::SEGTP::default();

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>>{

    let matches = build_cli();

    let query_fname = matches.try_get_one::<String>("query")?.expect("Query file is required");
    check_path(query_fname)?;

    let trail: Vec<_> = matches.get_many::<String>("templates").unwrap().collect();
    for tmpl_fname in trail{
        check_path(tmpl_fname)?;
    }

    let output_fname = matches.get_one::<String>("output").unwrap();

    run()?;

    Ok(())
}

fn check_path(path_fname: &str) -> Result<(),Box<dyn Error>>{
    let path_exists = Path::new(path_fname).try_exists().unwrap();
    match path_exists {
        true => Ok(()),
        false => {
            let err_msg = format!("Template path does not exist: {}", path_fname);
            Err(err_msg.into())
        },
    }
}

fn build_cli() -> clap::ArgMatches{
    Command::new("ORFanage")
        .version("1.0")
        .author("Ales Varabyou <ales.varabyou[at]jhu[dot]edu>")
        .about("Annotating Open Reading Frames based on reference, phylogeny and sequence similarity.")
        .arg(
            Arg::new("templates")
                .num_args(1..)
                .trailing_var_arg(true)
                .help("One or more GFF/GTF files with coding exons to be used as templates.")
                .required_unless_present_any(["help", "version"]),
        )
        .arg(
            Arg::new("query")
                .short('q')
                .long("query")
                .help("Path to a GTF query file with transcripts to which CDSs are to be ported")
                .required_unless_present_any(["help", "version"]),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Basename for all output files generated by this software")
                .required_unless_present_any(["help", "version"]),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .help("Filename to write stats file into")
                .required(false),
        )
        .arg(
            Arg::new("reference")
                .short('g')
                .long("reference")
                .help("Path to the reference genome file in FASTA format. This parameter is required \
                       when the following parameters are used: 1. cleanq; 2. cleant; 3. pd.")
                .required(false),
        )
        .arg(
            Arg::new("cleanq")
                .long("cleanq")
                .help("If enabled - will ensure all transcripts in the output file will have a \
                       valid start and end codons. This option requires the use of --reference parameter")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("cleant")
                .long("cleant")
                .help("If enabled - will ensure all ORFs in the reference annotations start with a \
                      valid start codon and end with the first available stop codon. This option \
                      requires the use of --reference parameter")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("rescue")
                .long("rescue")
                .help("If enabled - will attempt rescuing the broken ORFs in the reference annotations. \
                       This option requires the use of --reference parameter")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("lpd")
                .long("lpd")
                .help("Percent difference by length between the original and reference transcripts. \
                      If -1 (default) is set - the check will not be performed.")
                .required(false)
                .default_value("0")
                .value_parser(clap::value_parser!(u16).range(0..100)),
        )
        .arg(
            Arg::new("ilpd")
                .long("ilpd")
                .help("Percent difference by length of bases in frame of the reference transcript. \
                       If -1 (default) is set - the check will not be performed.")
                .required(false)
                .default_value("0")
                .value_parser(clap::value_parser!(u16).range(0..100)),
        )
        .arg(
            Arg::new("mlpd")
                .long("mlpd")
                .help("Percent difference by length of bases that are in both query and reference. \
                       If -1 (default) is set - the check will not be performed.")
                .required(false)
                .default_value("0")
                .value_parser(clap::value_parser!(u16).range(0..100)),
        )
        .arg(
            Arg::new("pi")
                .long("pi")
                .help("Percent identity between the query and template sequences. This
                      option requires --reference parameter to be set. If enabled -
                      will run alignment between passing pairs.")
                .required(false)
                .default_value("1")
                .value_parser(clap::value_parser!(u16).range(0..100)),
        )
        .arg(
            Arg::new("gapo")
                .long("gapo")
                .help("Gap-open penalty")
                .required(false)
                .default_value("4")
                .value_parser(clap::value_parser!(u16).range(0..)),
        )
        .arg(
            Arg::new("gape")
                .long("gape")
                .help("Gap-extension penalty")
                .required(false)
                .default_value("2")
                .value_parser(clap::value_parser!(u16).range(0..)),
        )
        .arg(
            Arg::new("minlen")
                .long("minlen")
                .help("Minimum length of an open reading frame to consider for the analysis")
                .required(false)
                .default_value("0")
                .value_parser(clap::value_parser!(u16).range(0..)),
        )
        .arg(
            Arg::new("mode")
                .short('m')
                .long("mode")
                .help("Which CDS to report")
                .required(false)
                .value_parser(["ALL","LONGEST","LONGEST_MATCH","BEST"])
                .default_value("LONGEST_MATCH"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use in the run")
                .required(false)
                .default_value("1")
                .value_parser(clap::value_parser!(u16).range(1..)),
        )
        .arg(
            Arg::new("useid")
                .long("useid")
                .help("If enabled, only transcripts with the same gene ID from the
                      query file will be used to form a bundle. In this mode the same
                      template transcript may be used in several bundles, if overlaps
                      transcripts with different gene_ids.")
                .required(false),
        )
        .arg(
            Arg::new("nonaug")
                .long("nonaug")
                .help("If enabled, non-AUG start codons in reference transcripts will
                      not be discarded and will be considered in overlapping query
                      transcripts on equal grounds with the AUG start codon.")
                .required(false),
        )
        .arg(
            Arg::new("keepcds")
                .long("keepcds")
                .help("If enabled, any CDS already presernt in the query will be kept
                      unmodified.")
                .required(false),
        ).get_matches()
}