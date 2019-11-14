#!/usr/bin/python3
############################################################
# For popgen, 10.19
# Takes a transformed FASTA file and computes some stats.
############################################################

import sys, os, argparse, math, re, lib.poplib as poplib
import multiprocessing as mp
from collections import defaultdict

# striatus,amoenus,cratericus,minimus,ruficaudus
# 8678588

# Current definitions:
# Total sites:          # of unambiguous sites in file (Ambiguous sites are skipped completely in all subsequent counts).
# Invariant sites:      Sites where all sequences in all groups share the same base (ie. TTT TTTT TTT TTTT TTT).
# Polymorphic sites:    Sites in which at least one sequence differs from another (ie AAA TTTA TTT TCCG TTT).
# Fixed differences:    Sites where all samples within this group share one allele while all other samples share a different allele (ie AAA TTTT TTT TTTT TTT).
# Polymorphic group:    Sites where at least one sample within a group has a different base than the rest, regardless of the alleles in the other groups (ie AAA TTTA TTT TCCG TTT is polymorphic in both groups 2 and 4).
# Shared fixed sites:   All samples in these two groups share the same allele, while all other samples share another allele (ie AAA AAAA TTT TTTT TTT).

# Currently not considering phylogeny.

############################################################

def dsum(*dicts):
# Given a list of dictionaries, this function merges them into a single dictionary, summing values of common keys.
    ret = defaultdict(int);
    for d in dicts:
        for k, v in d.items():
            ret[k] += v;
    return dict(ret);

################

def chunks(l, n):
# Splits a list l into even chunks of size n.
    n = max(1, n);
    return (l[i:i+n] for i in range(0, len(l), n));

################

def outMat(outdict, glist):
    outmat = [];
    headers = " " * pad;
    done = [];
    for group in glist:
        headers += poplib.spacedOut(group, pad);
        outline = poplib.spacedOut(group, pad);
        for gcomp in glist:
            if gcomp == group:
                outline += poplib.spacedOut("-", pad);
                outlines[group].append("NA");
                continue;
            
            outlines[group].append(str(outdict[group][gcomp]));
            p1 = group + "-" + gcomp;
            p2 = gcomp + "-" + group;
            if p1 in done or p2 in done:
                outline += poplib.spacedOut("", pad);
            else:
                outline += poplib.spacedOut(str(outdict[group][gcomp]), pad);
                
                done.append(p1);
                done.append(p2);
        outmat.append(outline);

    print(headers);
    for o in outmat:
        print(o);

################

def siteParse(chunk_item):
# For a given set of sites, this function parses each one and checks for site categories for each group (see Definitions above).
    chunk_lines, globs = chunk_item;

    gcounts = { g1 : { 'fixed' : 0, 'polymorphic' : 0, "biallelic" : 0, "triallelic" : 0, "quadallelic" : 0, 'h-sum' : 0.0 } for g1 in globs['group-list'] };
    sfixed = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    pfixed = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    ppi = { g1 : { g2 : 0.0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    spoly = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    tcounts = { 'sites' : 0, 'invariant' : 0, 'polymorphic' : 0, 'ambiguous' : 0, 'h-sum' : 0.0 };

    for site in chunk_lines:
        site = site.strip();
        if any(c in site for c in globs['ambiguous']):
            tcounts['ambiguous'] += 1;
            continue;
            #return (tcounts, gcounts, sfixed, spoly);
        #print(site);
        tcounts['sites'] += 1;
        site = [s.replace(s, globs['het-codes'][s]) if s in globs['het-codes'] else s for s in list(site)];
        #print(site);

        if site.count(site[0]) == len(site):
            tcounts['invariant'] += 1;
            continue;
            #return (tcounts, gcounts, sfixed, spoly);

        tcounts['polymorphic'] += 1;
        tcounts['h-sum'] += hetCalc(site, globs['het']);

        group_sites = { g : [] for g in globs['group-list'] };
        for group in globs['group-list']:
            for ind in globs['groups'][group]:
                #group_sites[group] += site[ind];
                group_sites[group].append(site[ind]);

        group_alleles = { g : ''.join(set(group_sites[g])) for g in globs['group-list'] };

        # print(site);
        # print(group_sites);
        # print(group_alleles);
        # sys.exit();

        for group in globs['group-list']:
            if len(group_alleles[group]) == 1:
                fixed = True;
                shared = [];
                for gcomp in globs['group-list']:
                    if gcomp == group:
                        continue;

                    paired_site = group_sites[group] + group_sites[gcomp];

                    if group_alleles[group] == group_alleles[gcomp]:
                        pfixed[group][gcomp] += 1;

                    if group_alleles[group] in group_alleles[gcomp]:
                        fixed = False; 
                        if group_alleles[group] == group_alleles[gcomp]:
                            shared.append(gcomp);

                    ppi[group][gcomp] += hetCalc(paired_site, globs['het']);

                if fixed:
                    gcounts[group]['fixed'] += 1;

                if len(shared) == 1:
                    sfixed[group][shared[0]] +=1;

            elif len(group_alleles[group]) > 1:
                gcounts[group]['polymorphic'] += 1;
                gcounts[group]['h-sum'] += hetCalc(group_sites[group], globs['het']);

    return (tcounts, gcounts, sfixed, pfixed, ppi, spoly);

################

def hetCalc(site, het_flag):
    #print(site);
    num_seqs = len(site);
    if het_flag:
        num_hets = len([ h for h in site if len(h) == 2 ]);
        cur_h = num_hets / num_seqs;
    else:
        #base_props_sq = { b : (site.count(b) / num_seqs)**2 for b in "ATCG" };
        base_sq_sum = sum( [ (site.count(b) / num_seqs)**2 for b in "ATCG" ] )
        cur_h = (num_seqs / (num_seqs -1)) * (1 - base_sq_sum);
    return cur_h;

def piCalc(site):
    bases = "ATCG";


#def popFunc():


############################################################


if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.
    globs = {
        'ambiguous' : "BDHVN.-",
        'het-codes' : { "R" : "AG", "Y" : "CT", "S" : "GC", "W" : "AT", "K" : "GT", "M" : "AC" }
    };
    globs['het-pattern'] = re.compile('|'.join(globs['het-codes'].keys()));
    # Some global variables.

    print("PROGRAM CALL: " + " ".join(sys.argv));

    parser = argparse.ArgumentParser(description="Quick population statistics for a FASTA file.");
    parser.add_argument("-i", dest="input", help="The path to a transformed FASTA alignment file (.tfa) with sequences from different groups.", default=False);
    parser.add_argument("-t", dest="titles", help="The path to the file that contains the FASTA titles/headers corresponding to the sequences in the .tfa file. Default: -i with .title extension.", default=False);
    parser.add_argument("-g", dest="groups", help="A comma separated list of group labels in the original FASTA file.", default=False);
    parser.add_argument("--het", dest="het_flag", help="Use this option if your data contains IUPAC ambiguity codes. Heterozygosity will be observed directly from the data rather than calculated as expected heterozygosity.", action="store_true", default=False);
    #parser.add_argument("-r", dest="outgroup", help="The group from -g that should be considered the outgroup. Defalt: none", default=False);
    # Add outgroup option for polarizing?!
    parser.add_argument("-p", dest="procs", help="The number of processes the script should use. Default: 1.", type=int, default=1);
    parser.add_argument("-o", dest="output", help="The output file name.", default=False);
    args = parser.parse_args();
    # Input options.

    globs['het'] = args.het_flag;

    if not os.path.isfile(args.input):
        sys.exit(poplib.errorOut(1, "Please provide a valid transformed FASTA file as input (-i)."));
    # if not os.path.isfile(args.output):
    #     sys.exit(poplib.errorOut(2, "Please provide a an output file name (-o)."));
    if not args.titles:
        args.titles = os.path.splitext(args.input)[0] + ".title";
    if not os.path.isfile(args.titles):
        sys.exit(poplib.errorOut(3, "Please provide a valid .title (-t) file with FASTA titles/headers corresponding to those in the input file."));
    if not args.groups:
        sys.exit(poplib.errorOut(4, "Please provide a comma separated list of sequence groups (-g)."));
    if args.procs < 1 or type(args.procs) != int:
        sys.exit(poplib.errorOut(5, "Number of processes (-p) must be an integer > 0."));
    # Make sure specified number of procs is a positive integer.

    titles = open(args.titles, "r").read().split("\n");
    titles = list(filter(None, titles));
    # Get the list of titles from the titles file.

    globs['group-list'] = args.groups.split(",");
    # Get the list of groups from input.
    #print(globs);
    #sys.exit();
    # if args.outgroup:
    #     if args.outgroup not in group_list:
    #         sys.exit(poplib.errorOut(6, "Could not find the specified outgroup (-r) in the group list (-g)."));
    #     else:
    #         group_list[group_list.index(args.outgroup)] = args.outgroup + "OUTGROUP32891";
    # Parse the outgroup (not currently supported).

    globs['groups'] = defaultdict(list);
    title_found = [];
    group_counts = { group : 0 for group in globs['group-list'] };
    title_counts = { title : 0 for title in titles };
    for group in globs['group-list']:
        for x in range(len(titles)):
            if group in titles[x]:
                globs['groups'][group].append(x);
                title_found.append(titles[x]);
                title_counts[titles[x]] += 1;
                group_counts[group] += 1;
    #print(group_counts);
    #print(title_counts);
    # Match up the titles in the title file to the groups in the group input.

    if any(g == 0 for g in group_counts.values()):
        print("** WARNING: The following groups were not found in any of the FASTA titles:");
        for group in group_counts:
            if group_counts[group] == 0:
                print("\t" + group);

    if any(t == 0 for t in title_counts.values()):
        print("** WARNING: The following titles contained none of the specified groups:");
        for title in title_counts:
            if title_counts[title] == 0:
                print("\t" + title);

    if any(t > 1 for t in title_counts.values()):
        print("** Error 7: The following titles contained multiple group labels:");
        for title in title_counts:
            if title_counts[title] > 1:
                print("\t" + title);
        sys.exit();
    #print(groups);
    # Checking group and title overlap

    group_counts = { g1 : { 'fixed' : 0, 'polymorphic' : 0, "biallelic" : 0, "triallelic" : 0, "quadallelic" : 0, 'h-sum' : 0.0 } for g1 in globs['group-list'] };
    shared_fixed = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    paired_fixed = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    paired_pi = { g1 : { g2 : 0.0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    shared_polymorphic = { g1 : { g2 : 0 for g2 in globs['group-list'] if g2 != g1 } for g1 in globs['group-list'] };
    totals = { 'sites' : 0, 'invariant' : 0, 'polymorphic' : 0, 'ambiguous' : 0, 'h-sum' : 0.0 };
    # All the counts of different site categories.

    lines_per_proc = 10000;
    pool = mp.Pool(processes=args.procs);
    cur_lines = [];
    i = 0;
    i_start = 1;
    for line in open(args.input, "r"):
        i += 1;
        cur_lines.append(line);
        if len(cur_lines) == args.procs*lines_per_proc:
            print("Processing sites " + str(i_start) + "-" + str(i));
            i_start = i + 1;
            line_chunks = list(chunks(cur_lines, lines_per_proc));
            for result in pool.imap(siteParse, ((line_chunk, globs) for line_chunk in line_chunks)):
                totals = dsum(totals, result[0]);
                for group in globs['group-list']:
                    group_counts[group] = dsum(group_counts[group], result[1][group]);
                    shared_fixed[group] = dsum(shared_fixed[group], result[2][group]);
                    paired_fixed[group] = dsum(paired_fixed[group], result[3][group]);
                    paired_pi[group] = dsum(paired_pi[group], result[4][group]);
                    shared_polymorphic[group] = dsum(shared_polymorphic[group], result[5][group]);
            cur_lines = [];   
    # Read the input file line by line. Once a certain number of lines have been read, pass them to siteParse in parallel.

    if cur_lines != []:
        print("Processing sites " + str(i_start) + "-" + str(i));
        i_start = i + 1;
        line_chunks = list(chunks(cur_lines, lines_per_proc));
        for result in pool.imap(siteParse, ((line_chunk, globs) for line_chunk in line_chunks)):
            totals = dsum(totals, result[0]);
            for group in globs['group-list']:
                group_counts[group] = dsum(group_counts[group], result[1][group]);
                shared_fixed[group] = dsum(shared_fixed[group], result[2][group]);
                paired_fixed[group] = dsum(paired_fixed[group], result[3][group]);
                paired_pi[group] = dsum(paired_pi[group], result[4][group]);
                shared_polymorphic[group] = dsum(shared_polymorphic[group], result[5][group]);
    # Count the last chunk of lines if necessary.

    outlines = { g : [] for g in globs['group-list'] };
    outlines["total"] = [];
    # Initialize output for output file

    total_pad = 20;
    print();
    print(poplib.spacedOut("Total sites:", total_pad) + str(totals['sites']));
    print(poplib.spacedOut("Invariant sites:", total_pad) + str(totals['invariant']));
    print(poplib.spacedOut("Polymorphic sites:", total_pad) + str(totals['polymorphic']));
    print(poplib.spacedOut("Pi:", total_pad) + str(round(totals['h-sum'], 3)));
    print(poplib.spacedOut("Pi per site:", total_pad) + str(round(totals['h-sum'] / (totals['sites'] - totals['ambiguous']), 5)));
    # Print the stats for total sites.

    pad = 12;
    print();
    headers = poplib.spacedOut("Group", pad) + poplib.spacedOut("Fixed", pad) + poplib.spacedOut("Polymorphic", pad) + poplib.spacedOut("Pi", pad) + "Pi per site";
    print(headers);
    for group in globs['group-list']:
        outlines[group].append(group);
        outline = poplib.spacedOut(group, pad) + poplib.spacedOut(str(group_counts[group]['fixed']), pad) + poplib.spacedOut(str(group_counts[group]['polymorphic']), pad);
        outline += poplib.spacedOut(str(round(group_counts[group]['h-sum'], 3)), pad);
        outline += poplib.spacedOut(str(round(group_counts[group]['h-sum'] / (totals['sites'] - totals['ambiguous']), 5)), pad);
        # Gotta count sites per group!!
        outlines[group].append(str(group_counts[group]['fixed']));
        outlines[group].append(str(group_counts[group]['polymorphic']));



        print(outline);
    print("");
    # Print the stats for individual groups.

    print("Shared fixed:");
    outMat(shared_fixed, globs['group-list']);

    print("\nPaired fixed:");
    outMat(paired_fixed, globs['group-list']);

    print("\nPaired pi:");
    paired_pi_site = { group : {} for group in globs['group-list'] };
    for group in globs['group-list']:
        for gcomp in globs['group-list']:
            if group == gcomp:
                continue;
            paired_pi_site[group][gcomp] = round(paired_pi[group][gcomp] / (totals['sites'] - totals['ambiguous']), 5);
            paired_pi[group][gcomp] = round(paired_pi[group][gcomp], 5);
    outMat(paired_pi, globs['group-list']);
    print("\nPaired pi per site:");
    outMat(paired_pi_site, globs['group-list']);
    # Construct and print the pairwise shared matrix string.

    print()


    if args.output:
        outheaders = ",".join(["Group","Fixed","Polymorphic"] + globs['group-list']);
        with open(args.output, "w") as outfile:
            outfile.write(outheaders + "\n");
            for g in globs['group-list']:
                outfile.write(",".join(outlines[g]) + "\n");
    # Output to the output file if specified.





'''
line_bytes = [];
print("Reading file...");
lines = open(args.input, "r").read().split("\n");
lines = list(filter(None, lines));

lines_per_chunk = int(math.ceil(len(lines) / args.procs));
print("# " + poplib.getTime() + " Splitting into " + str(args.procs) + " chunks (" + str(lines_per_chunk) + " reads per chunk)");
line_chunks = list(chunks(lines, lines_per_chunk));
del(lines);

i = 0;
pool = mp.Pool(processes=args.procs);
# Declare the processor pool.
#for result in pool.imap(siteParse, ((line, group_list) for line in open(args.input))):
for result in pool.imap(siteParse, ((line_chunk, group_list) for line_chunk in line_chunks)):
#for line in open(args.input):
    # print(line.strip());
    # print(totals);
    #result = siteParse((line, group_list));
    totals = dsum(totals, result[0]);
    # print(totals);
    for group in group_list:
        # if totals['polymorphic'] != 0:
        #     print(group);
        #     print(group_counts[group]);
        #     print(result[1][group]);
        #     print("------");
        group_counts[group] = dsum(group_counts[group], result[1][group]);
        shared_fixed[group] = dsum(shared_fixed[group], result[2][group]);
        shared_polymorphic[group] = dsum(shared_polymorphic[group], result[3][group]);
    #     if totals['polymorphic'] != 0:    
    #         print(group_counts[group]);
    #         print(result[1][group]);
    # if totals['polymorphic'] != 0:
    #     break;

#     #print(type(group_counts));
    #group_counts = dsum(group_counts, result[1]);

# THIS WORKS BUT READS THE WHOLE FILE INTO MEMORY

'''



'''
if i % 10000 == 0:
        print(i);
    i += 1;
    line = line.strip();
    if any(c in line for c in "RYSWKMBDHVN.-"):
        totals['ambiguous'] += 1;
        continue;

    numsites += 1;
    
    # print(line.count(line[0]));
    # print(len(line));
    if line.count(line[0]) == len(line):
        # print(line);
        # sys.exit();
        totals['invariant'] += 1;
        continue;

    totals['polymorphic'] += 1;
    #totals['pi'] += piCalc(line);

    group_sites = { g : "" for g in group_list };
    for group in group_list:
        for ind in groups[group]:
            group_sites[group] += line[ind];


    group_alleles = { g : ''.join(set(group_sites[g])) for g in group_list };

    # print(group_sites);
    # print(group_alleles);


    for group in group_list:

        if len(group_alleles[group]) == 1:
            fixed = True;
            shared = [];
            for gcomp in group_list:
                if gcomp == group:
                    continue;
                if group_alleles[group] in group_alleles[gcomp]:
                    fixed = False; 
                    if group_alleles[group] == group_alleles[gcomp]:
                        shared.append(gcomp);


            if fixed:
                group_counts[group]['fixed'] += 1;

            if len(shared) == 1:
                group_counts[group]['shared-fixed'][shared[0]] +=1;





        elif len(group_alleles[group]) > 1:
            group_counts[group]['polymorphic'] += 1;
            # print(line);
            # print(group_sites);
            # print(group_alleles);
            # sys.exit();
            # if len(group_alleles[group]) == 2:
            #     group_counts[group]['biallelic'] += 1;
            # elif len(group_alleles[group]) == 3:
            #     group_counts[group]['triallelic'] += 1;
            # elif len(group_alleles[group]) == 4:
            #     group_counts[group]['quadallelic'] += 1;

    #print(group_counts);


# for group in groups:
#     print(group);
#     print(group_counts[group]);

# print(totals);

with open(args.input, "r") as infile:
    first = True;
    line = "sadf;";
    while line:
        if line == "":
            continue;
        line = infile.readline();
        line_start = infile.tell() - len(line);
        line_end = infile.tell() - 1;

        line_bytes.append((line_start, line_end));

    
i = 1;
pool = mp.Pool(processes=args.procs);
for result in pool.imap(getLine, ((line_byte, args.input) for line_byte in line_bytes)):
    #print(getLine(line_byte, infile));
    print(result);
    #i += 1;
    #if i > 4:
sys.exit();

################

def getLine(line_info):
    lb, infile = line_info;
    f = open(infile, "r");
    f.seek(lb[0]);
    line = f.read(lb[1] - lb[0]).strip();
    f.close();
    return line;


'''