#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by Ales Varabyou
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <gff_utils.h>
#include <unordered_set>

#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #include <omp.h>
    #endif
#endif

#include "arg_parse.hpp"
#include "Transcriptome.h"

struct Parameters{
    int num_threads = 1;
    bool use_id = false; // if enabled will use query gene ids to form bundles. the same reference id may be evaluated in multiple bundles if overlaps queries with different gene ids

    std::string query_fname;
    std::string template_fname;
    std::string reference_fasta_fname;
    std::string output_fname;
    std::ofstream out_fp;
} global_params;

std::string get_score_string(Score& s) {

}

int run(){
    #ifdef DEBUG
        std::cout<<"Running in DEBUG"<<std::endl;
    #endif
    Transcriptome transcriptome;

    if(!global_params.reference_fasta_fname.empty()){
        std::cerr<<"loading reference genome"<<std::endl;
        transcriptome.set_ref(global_params.reference_fasta_fname);
    }

    std::cerr<<"loading reference transcriptomes"<<std::endl;
    transcriptome.add(global_params.template_fname,true,true);
    transcriptome.add(global_params.query_fname,false,true);

    transcriptome.build_cds_chains();

    std::cerr<<"sorting reference transcriptome"<<std::endl;
    transcriptome.set_cds_as_exons();
    transcriptome.remove_non_coding();

    // sort again since now cds is set as exons
    transcriptome.sort();

    std::cerr<<"start removing duplicates"<<std::endl;
    int template_duplicates = transcriptome.deduplicate(global_params.use_id);

    std::cerr<<"bundling transcriptome"<<std::endl;
    transcriptome.bundleup(global_params.use_id);

    std::cerr<<"starting main evaluation"<<std::endl;

#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #pragma omp parallel for schedule(dynamic)
    #endif
#endif
    for(int bi=0;bi<transcriptome.bsize();bi++){ // iterate over bundles
        std::vector<Bundle>::iterator bundle_it = transcriptome.bbegin();
        bundle_it+=bi;
        TX *q,*t;

        std::string cur_seqid;
        transcriptome.seqid2name(bundle_it->get_seqid(),cur_seqid);

        if(!bundle_it->has_template() && !bundle_it->has_query()){
            continue;
        }
        else{
            // iterate over queries
            for(int qi=0;qi<bundle_it->size();qi++) {
                q = bundle_it->operator[](qi);
                if(q->is_template()){continue;}
                
                q->build_cds();

                std::string q_aa = q->get_aa();
                char q_start_codon = '-';
                char q_stop_codon = '-';
                if(!global_params.reference_fasta_fname.empty()){
                    q->load_seq();
                    q_aa = q->get_aa();
                    q_start_codon = q_aa.front();
                    q_stop_codon = q_aa.back();
                }

                if (!bundle_it->has_template()) {
                    global_params.out_fp<<
                                        q->get_tid()<<"\t"
                                        "-\t"
                                        <<q->cds_len()<<"\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        "-\t"
                                        <<q_start_codon<<"\t"
                                        "-\t"
                                        <<q_stop_codon<<"\t"
                                        "-\t"
                                        <<std::endl;
                    continue;
                }
                
                // iterate over templates
                for(int ti=0;ti<bundle_it->size();ti++){
                    t=bundle_it->operator[](ti);
                    if(!t->is_template()){continue;}

                    t->build_cds();

                    std::string t_aa = t->get_aa();
                    char t_start_codon = '-';
                    char t_stop_codon = '-';
                    if(!global_params.reference_fasta_fname.empty()){
                        t->load_seq();
                        t_aa = t->get_aa();
                        t_start_codon = t_aa.front();
                        t_stop_codon = t_aa.back();
                    }

                    Score score = q->score(*t);

                    // we now have two transcripts with CDSs
                    // compare to each other and write stats to the output
                    global_params.out_fp<<q->get_tid()<<"\t"
                                        <<t->get_tid()<<"\t"
                                        <<q->cds_len()<<"\t"
                                        <<t->cds_len()<<"\t"
                                        <<score.num_bp_match<<"\t"
                                        <<score.num_bp_inframe<<"\t"
                                        <<score.num_bp_outframe<<"\t"
                                        <<score.num_bp_extra<<"\t"
                                        <<score.num_bp_missing<<"\t"
                                        <<score.lpi<<"\t"
                                        <<score.mlpi<<"\t"
                                        <<score.ilpi<<"\t"
                                        <<q_start_codon<<"\t"
                                        <<t_start_codon<<"\t"
                                        <<q_stop_codon<<"\t"
                                        <<t_stop_codon<<"\t"
                                        <<std::endl;
                }
            }
        }
    }

    std::cerr<<"Done"<<std::endl;

    return 0;
}

int main(int argc, char** argv) {
    ArgParse args("orfcompare",
                  "run comparison of ORFs between two files.");

    args.add_option("query",ArgParse::Type::STRING, "Path to a GTF query file with transcripts to which CDSs are to be ported", ArgParse::Level::GENERAL,true);
    args.add_option("template",ArgParse::Type::STRING,"Path to the reference genome file in FASTA format. This parameter is required when the following parameters are used: 1. cleanq; 2. cleant; 3. pi.",ArgParse::Level::GENERAL,true);
    args.add_option("output", ArgParse::Type::STRING, "Basename for all output files generated by this software", ArgParse::Level::GENERAL,true);
    args.add_option("threads",ArgParse::Type::INT,"Number of threads to run in parallel",ArgParse::Level::GENERAL,false);
    args.add_option("use_id",ArgParse::Type::FLAG,"If enabled, only transcripts with the same gene ID from the query file will be used to form a bundle. In this mode the same template transcript may be used in several bundles, if overlaps transcripts with different gene_ids.",ArgParse::Level::GENERAL,false);
    args.add_option("reference",ArgParse::Type::STRING,"Path to the reference genome file in FASTA format. This parameter is required when the following parameters are used: 1. cleanq; 2. cleant; 3. pi.",ArgParse::Level::GENERAL,false);
    
    args.add_option("help", ArgParse::Type::FLAG, "Prints this help message.", ArgParse::Level::HELP, false);

    args.parse_args(argc, argv);
    args.check_args();

    // first create the execution string
    std::string cl = "orfcompare ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    // check that all files exist
    global_params.query_fname = args.get_string("query");
    if(!file_exists(global_params.query_fname)){
        std::cerr << "Query file does not exist! "<<args.get_string("query")<<std::endl;
        exit(2);
    }

    global_params.template_fname = args.get_string("template");
    if(!file_exists(global_params.template_fname)){
        std::cerr << "Template file does not exist! "<<args.get_string("template")<<std::endl;
        exit(2);
    }

    global_params.num_threads = args.is_set("threads") ? args.get_int("threads") : 1;
#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        omp_set_num_threads(global_params.num_threads);
    #endif
#endif

    if(args.is_set("use_id")){
        global_params.use_id = true;
    }

    if(args.is_set("reference")){ // much from gpertea bamcons.cpp
        global_params.reference_fasta_fname = args.get_string("reference");
        if(!file_exists(global_params.reference_fasta_fname)){
            std::cerr << "Reference FASTA file does not exist!\n";
            exit(2);
        }

        // get potential fasta index file name
        std::string fa_idx_fname = args.get_string("reference")+".fai";
        GFastaIndex faIdx(args.get_string("reference").c_str(), fa_idx_fname.c_str());
        if (!faIdx.hasIndex()){
            std::cerr<<"No fasta index found for "<<fa_idx_fname<<". Building now"<<std::endl;
            faIdx.buildIndex();
            if (faIdx.getCount() == 0){
                std::cerr<<"Error: no fasta records found!"<<std::endl;
                exit(2);
            }
            FILE* fcreate = fopen(fa_idx_fname.c_str(), "w");
            if (fcreate == nullptr){
                std::cerr<<"Error: cannot create fasta index: "<<fa_idx_fname<<std::endl;
                exit(2);
            }
            if (faIdx.storeIndex(fcreate) < faIdx.getCount()){
                std::cerr<<"Warning: error writing the index file!"<<std::endl;
                exit(2);
            }
            std::cerr<<"FASTA index rebuilt."<<std::endl;
        }
    }

    global_params.output_fname = args.get_string("output");
    global_params.out_fp.open(global_params.output_fname);
    Score s;
    global_params.out_fp<<"query_id\t"
                          "template_id\t"
                          "query_len\t"
                          "template_len\t"
                          "len_match\t"
                          "len_inframe\t"
                          "len_outframe\t"
                          "len_extra\t"
                          "len_missing\t"
                          "lpi\t"
                          "mlpi\t"
                          "ilpi\t"
                          "query_start_codon\t"
                          "template_start_codon\t"
                          "query_stop_codon\t"
                          "template_stop_codon"<<std::endl;

    run();

    global_params.out_fp.close();

    return 0;
}
#pragma clang diagnostic pop
