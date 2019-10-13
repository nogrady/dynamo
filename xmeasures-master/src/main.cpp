//! \brief Extrinsic measures evaluation for overlapping multi-resolution clusterings
//! with possible unequal node base.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-13

#include <cstdio>
#include <sstream>
#include "cmdline.h"  // Arguments parsing
#include "macrodef.h"
#include "interface.hpp"

using std::stringstream;


int main(int argc, char **argv)
{
	gengetopt_args_info  args_info;
	auto  err = cmdline_parser(argc, argv, &args_info);
	if(err)
		return err;

	// Validate required xmeasure
	if(!args_info.omega_flag && !args_info.nmi_flag && !args_info.f1_given && !args_info.label_given) {
		fputs("WARNING, no any measures to evaluate are specified\n", stderr);
		cmdline_parser_print_help();
		return EINVAL;
	}

	if(args_info.membership_arg <= 0) {
		fprintf(stderr, "ERROR, positive membership is expected: %G\n", args_info.membership_arg);
		return EDOM;
	}

	{	// Validate the number of input files
		// Note: sync_arg is specified if sync_given
		const auto  inpfiles = args_info.inputs_num + (args_info.sync_given || args_info.label_given);  // The number of input files
		if(inpfiles < 2 || inpfiles > 2 + args_info.sync_given + args_info.label_given) {
			fputs("ERROR, 2 input clusterings are required with possibly additional"
				" node base and clusters labels, i.e. 2-4 input files in total\n", stderr);
			cmdline_parser_print_help();
			return EINVAL;
		}
	}

	// Verify that labeled clusters correspond to the node base if any of them is specified
	if(args_info.sync_given && args_info.label_given && (strcmp(args_info.sync_arg, args_info.label_arg)
	|| (args_info.inputs_num == 2 && strcmp(args_info.sync_arg, args_info.inputs[0]))))
		throw invalid_argument("ERROR, node base file should correspond to the labeled clusters and"
			" represent the first evaluating collection if both are specified\n");

	// Load node base if required
	NodeBase  ndbase;
	::AggHash  nbhash;
	// Note: if label_given then either inputs_num < 2 or inputs_num[0] = sync_arg = label_arg
	if(args_info.sync_given && args_info.inputs_num == 2 && !args_info.label_given)
		ndbase = NodeBase::load(args_info.sync_arg, args_info.membership_arg
			, &nbhash, 0, 0, args_info.detailed_flag);

	auto process = [&](auto evaluation) -> int {
		using Count = decltype(evaluation);
		using Collection = Collection<Count>;
		// Load collections as relations
		::AggHash  cn1hash, cn2hash;
		// Note: cn1 is nodebase if specified and not in the separated file
		const bool  cn1base = (args_info.sync_given || args_info.label_given) && args_info.inputs_num < 2;
		//const char*  nbfile = args_info.sync_given
		auto cn1 = Collection::load(cn1base ? args_info.sync_given ? args_info.sync_arg
			: args_info.label_arg : args_info.inputs[0]
			, args_info.unique_flag, args_info.membership_arg, &cn1hash
			, ndbase ? &ndbase : nullptr, nullptr, args_info.detailed_flag);
		if(ndbase) {
			if(nbhash != cn1hash) {
				fprintf(stderr, "ERROR, nodebase hash %lu (%lu nodes) != filtered"
					" collection nodes hash %lu (%lu)\n", nbhash.hash(), nbhash.size()
					, cn1hash.hash(), cn1hash.size());
				return EINVAL;
			}
			ndbase.clear();
		}
		RawIds  lostcls;
		auto cn2 = Collection::load(args_info.inputs[!cn1base]
			, args_info.unique_flag, args_info.membership_arg, &cn2hash
			, args_info.sync_given ? &cn1 : nullptr
			, args_info.sync_given && args_info.label_given ? &lostcls : nullptr
			, args_info.detailed_flag);

		if(!cn1.ndsnum() || ! cn2.ndsnum()) {
			fprintf(stderr, "WARNING, at least one of the collections is empty, there is nothing"
				" to evaluate. Collection nodes sizes: %u, %u\n", cn1.ndsnum(), cn2.ndsnum());
			return EINVAL;
		}

		// Check the collections' nodebase
		if(cn1hash != cn2hash) {
			fprintf(stderr, "WARNING, the nodes in the collections differ (the quality will be penalized)"
				": %u nodes with hash %lu, size: %lu, ids: %lu, id2s: %lu) !="
				" %u nodes with hash %lu, size: %lu, ids: %lu, id2s: %lu);  synchronize: %s, label: %s\n"
				, cn1.ndsnum(), cn1hash.hash(), cn1hash.size(), cn1hash.idsum(), cn1hash.id2sum()
				, cn2.ndsnum(), cn2hash.hash(), cn2hash.size(), cn2hash.idsum(), cn2hash.id2sum()
				, daoc::toYesNo(args_info.sync_given), daoc::toYesNo(args_info.label_given));
			//if(args_info.sync_given) {
			//	fputs("ERROR, the nodes base should be synchronized\n", stderr);
			//	return EINVAL;
			//}
		}

		// The number of outputting measures (1 .. 4)
		uint8_t  outsnum = args_info.omega_flag + args_info.nmi_flag
			 + args_info.f1_given + args_info.label_given;
		stringstream  aggouts;  // Aggregated outputs
		// Evaluate and output measures
		// Note: evaluation of overlapping F1 after NMI allows to reuse some
		// calculations, for other cases the order of evaluations does not matter
		puts(string("= ").append(is_floating_point<Count>::value
			? "Overlaps" : "Multi-resolution").append(" Evaluation =").c_str());
		if(args_info.nmi_flag) {
			auto rnmi = Collection::nmi(cn1, cn2, args_info.ln_flag, args_info.detailed_flag);
			// Set NMI to NULL if collections have no any mutual information
			// ATTENTION: for some cases, for example when one of the collections is a single cluster,
			// NMI will always yield 0 for any clusters in the second collection, which is limitation
			// of the original NMI measure. Similar issues possible in more complex configurations.
			if(rnmi.mi <= precision_limit<decltype(rnmi.mi)>()) {  // Note: strict ! is fine here
				throw domain_error("NMI is not applicable to the specified collections: 0, which says nothing about the similarity\n");
				rnmi.h1 = rnmi.h2 = 1;
			}
			const auto  nmix = rnmi.mi / std::max(rnmi.h1, rnmi.h2);
			if(args_info.all_flag) {
				printf("NMI_max: %G, NMI_sqrt: %G, NMI_avg: %G, NMI_min: %G\n"
					, nmix, rnmi.mi / sqrt(rnmi.h1 * rnmi.h2)
					, 2 * rnmi.mi / (rnmi.h1 + rnmi.h2)
					, rnmi.mi / std::min(rnmi.h1, rnmi.h2));
				if(--outsnum || aggouts.tellp())
					aggouts << "NMI_max: " << nmix
						<< ", NMI_sqrt: " << rnmi.mi / sqrt(rnmi.h1 * rnmi.h2)
						<< ", NMI_avg: " << 2 * rnmi.mi / (rnmi.h1 + rnmi.h2)
						<< ", NMI_min: " << rnmi.mi / std::min(rnmi.h1, rnmi.h2);
			} else {
				printf("NMI_max:\n%G\n", nmix);
				if(--outsnum || aggouts.tellp())
					aggouts << "NMI_max: " << nmix;
			}
		}
		if(args_info.f1_given) {
			// Assign required F1 type
			F1  f1kind = F1::NONE;
			// Note: args_info.f1_orig is empty if default value is used
			char  f1suf = '-';  // Suffix char of the selected F1 measure
			switch(args_info.f1_arg) {
			case f1_arg_partprob:
				f1kind = F1::PARTPROB;
				f1suf = 'p';
				break;
			case f1_arg_harmonic:
				f1kind = F1::HARMONIC;
				f1suf = 'h';
				break;
			case f1_arg_average:
				f1kind = F1::AVERAGE;  // Suggested by Leskovec
				f1suf = 'a';
				break;
			default:
				throw invalid_argument("main(), UNKNOWN F1 policy specified\n");
			}
			// Assign matching kind
			Match  mkind = Match::NONE;
			// Note: args_info.kind_orig is empty if default value is used
			char  kindsuf = '-';  // Suffix char of the selected F1 measure
			switch(args_info.kind_arg) {
			case kind_arg_weighted:
				mkind = Match::WEIGHTED;
				kindsuf = 'w';
				break;
			case kind_arg_unweighed:
				mkind = Match::UNWEIGHTED;
				kindsuf = 'u';
				break;
			case kind_arg_combined:
				mkind = Match::COMBINED;
				kindsuf = 'c';
				break;
			default:
				throw invalid_argument("main(), UNKNOWN Matching policy specified\n");
			}

			//if(args_info.nmi_flag)
			//	fputs("; ", stdout);
			Prob  prc, rec;  // Precision and recall of cn2 relative to ground-truth cn1
			const auto  f1val = Collection::f1(cn1, cn2, f1kind, rec, prc, mkind, args_info.detailed_flag);
			printf("MF1%c_%c (%s, %s):\n%G", f1suf, kindsuf, to_string(f1kind).c_str()
				, to_string(mkind).c_str(), f1val);
			if(prc || rec)
				printf(" (Prc: %G, Rec: %G)", prc, rec);
			fputc('\n', stdout);
			if(--outsnum || aggouts.tellp()) {
				if(aggouts.tellp())
					aggouts << "; ";
				aggouts << "MF1" << f1suf << '_' << kindsuf << ": " << f1val;
				// Note: prc and rec are zeroized if the matching strategy does not support them
				if(prc || rec)
					aggouts << " (Prc: " << prc << ", Rec: " << rec << ')';
			}
		}
		// Label clusters with the ground-truth clusters indices and output F1 for the labels if required
		if(args_info.label_given) {
			if(args_info.policy_arg == policy__NULL) {
				fputs("WARNING f1(), labels matching policy is not specified, the evaluation is skipped\n", stderr);
				return 0;
			}
			// Reset cluster counters if they were set (could be set only by F1)
			if(args_info.f1_given) {
				cn1.clearcounts();
				cn2.clearcounts();
			}
			const bool  prob = args_info.policy_arg == policy_arg_partprob;  // Partial Probabilities matching policy
			const bool  weighted = !args_info.unweighted_flag;
			PrcRec pr = Collection::label(cn1, cn2, lostcls, prob, weighted
				, args_info.identifiers_arg, args_info.detailed_flag);
			// Note: each measure name should form a single world to be properly parsed in a uniform way (see Clubmark),
			// that is why doubled underscore is used rather than a single space.
			printf("F1%c_%c__labels: %G (Prc: %G, Rec: %G)\n"
				, prob ? 'p' : 'h', weighted ? 'w' : 'u'
				, hmean(pr.prc, pr.rec), pr.prc, pr.rec);
			if(--outsnum || aggouts.tellp()) {
				if(aggouts.tellp())
					aggouts << "; ";
				aggouts << "F1" << (prob ? 'p' : 'h') << '_' << (weighted ? 'w' : 'u')
					<< "__labels: " << hmean(pr.prc, pr.rec)
					<< " (Prc: " << pr.prc << ", Rec: " << pr.rec << ')';
			}
		}
		if(args_info.omega_flag) {
			// Transform loaded and pre-processed collection to the representation
			// suitable for Omega Index evaluation
			RawClusters  cls1;
			RawClusters  cls2;
			NodeRClusters  ndrcs;

			cn1.template transfer<true>(cls1, ndrcs);
			cn2.template transfer<false>(cls2, ndrcs);
			const auto oi = args_info.extended_flag
				? omega<true>(ndrcs, cls1, cls2)
				: omega<false>(ndrcs, cls1, cls2)
				;
			printf("OI%s:\n%G\n", args_info.extended_flag ? "x" : "", oi);
			if(--outsnum || aggouts.tellp()) {
				if(aggouts.tellp())
					aggouts << "; ";
				aggouts << "OI" << (args_info.extended_flag ? "x" : "") << ": " << oi;
			}
		}
		if(aggouts.tellp())
			puts(aggouts.str().c_str());

		return 0;
	};


    return args_info.ovp_flag ? process(AccProb()) : process(Id());
}
