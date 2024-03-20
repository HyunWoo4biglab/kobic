"""A script for processing vcf files
    contains get_variants : making a variant dictionary
    get_genotype : making a genotype dictionary
    variants_stats: a dictionary of number of variants by each SV type
"""


# A script for processing VCF file
# 03.14.2023



def get_variants(vcf, caller, trio=False):
    """Return a dictionary containing all variants from a vcf file
    in the form of {variant_id : (chromosome, sv_filter_falg, sv_start, sv_end, sv_type)}
    if trio=True: return a genotype dictionary {variant_id : {sample1 : gt, sample2 : gt}}
    if sv_type == BND: sv_end is set to sv_start + 1"""

    variants = dict(); genotype = dict(); bps= dict()
    with open(vcf, 'r') as infile:
        lines = [variant for variant in infile if not variant.startswith('#')]

        for line in lines:
            col = line.strip().split('\t')
            contig = col[0];
            if not contig.startswith('chr'):
                contig = 'chr' + contig
            #record_match = re.search('(chr[\dXY]+):(\d+)[\[\]]', col[4])
            record_match = re.search('[\[\]]([a-zA-Z]+.+):(\d+)[\[\]]', col[4])
            #record_match = re.search('(chr.+):(\d+)[\[\]]', col[4])
            if record_match:
                mate_contig = record_match.group(1)
            else:
                mate_contig = contig
            if re.match('(chr[1-9][0-9]?|chrX|chrY)\\b', contig) and re.match('(chr[1-9][0-9]?|chrX|chrY)\\b', mate_contig):
                sv_start = int(col[1]); 
                sv_id = col[2]
                sv_filter = col[6]
                info = col[7]
                if caller == 'lumpy':
                    sv_filter = 'PASS'
                    #if 'IMPRECISE' in info:
                    #    sv_filter = 'IMPRECISE'
                    #else:
                    #    sv_filter = 'PASS'

                #genotype[sv_id] = dict()

                #sv_type = info.split("SIMPLE_TYPE=")[-1].split(';')[0]
                if caller == 'gridss':
                    if 'SIMPLE_TYPE' in info:
                        sv_type = info.split("SIMPLE_TYPE=")[-1]
                    else:
                        sv_type = 'NA'
                    if 'MATEID' in info:
                        mate_id = info.split('MATEID=')[-1].split(';')[0]
                    else:
                        mate_id = 'NA'
                    #sv_id += ';' + mate_id
                elif caller == 'svaba' or caller == 'svABA':
                    sv_type = col[-1]
                    if 'MATEID' in info:
                        mate_id = info.split('MATEID=')[-1].split(';')[0]
                    else:
                        mate_id = 'NA'
                else:
                    sv_type = info.split("SVTYPE=")[-1].split(';')[0]
                if sv_type != "BND":
                    if caller == 'gridss' or caller == 'svaba':
                        pass
                    #sv_end = int(info.split("END=")[1].split(';')[0])
                    else:
                        if caller == 'manta':
                            sv_end = int(info.split(';')[0].split("END=")[-1])
                        else:
                            sv_end = int(info.split(";END=")[-1].split(';')[0])
                else:
                    sv_end = sv_start + 1

                if trio:
                    if caller == 'svaba':
                        child_gt = col[-4].split(':')[0]
                        father_gt = col[-3].split(':')[0]
                        mother_gt = col[-2].split(':')[0]
                    else:
                        child_gt = col[-3].split(':')[0]
                        father_gt = col[-2].split(':')[0]
                        mother_gt = col[-1].split(':')[0]
                    variants[sv_id] = (contig, sv_filter, sv_start, sv_end, sv_type, child_gt, father_gt, mother_gt)

                else:
                    if caller == 'gridss' or caller == 'svaba':
                        bps[sv_id] = (contig, sv_filter, sv_start, mate_id, sv_type)
                    elif caller == 'svaba':
                        bps[sv_id] = (contig, sv_filter, sv_start, mate_id, sv_type)
                    else:
                        variants[sv_id] = (contig, sv_filter, sv_start, sv_end, sv_type, 'NA', 'NA', 'NA')

            if caller == 'gridss' or caller == 'svaba':
                for bp in bps:
                    contig = bps[bp][0]
                    sv_filter = bps[bp][1]
                    sv_start = bps[bp][2]
                    mate_id = bps[bp][3]
                    sv_type = bps[bp][4]
                    #sv_end = bps[mate_id][2]
                    #sv_id = bp + '-' + mate_id
                    if mate_id == 'NA':
                        sv_id = bp
                        sv_end = sv_start + 1
                    else:
                        if mate_id + '-' + bp in variants:
                            pass
                        else:
                            if mate_id in bps:
                                sv_id = bp + '-' + mate_id
                                sv_end = bps[mate_id][2]
                            else:
                                sv_id = bp
                                sv_end = sv_start + 1
                    if trio:
                        variants[sv_id] = (contig, sv_filter, sv_start, sv_end, sv_type, child_gt, father_gt, mother_gt)
                    else:
                        variants[sv_id] = (contig, sv_filter, sv_start, sv_end, sv_type, 'NA', 'NA', 'NA')
        return variants


def get_genotype(vcf):
    """Return 
    {sv_id : (son_gt, father_gt, mother_gt)"""

    genotypes = dict()
    gtlst = []; unique_gt = set()
    with open(vcf, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    samples = line.rstrip('\n').split('FORMAT\t')[-1].split('\t')
                    n_sample = len(samples)
            else:
                col = line.split('\t')
                sv_id = col[2]
                #gtlst = (col[i-n_sample].split(':')[0] for i in range(n_sample))
                for i in range(n_sample):
                    gt = col[i-n_sample].split(':')[0]
                    unique_gt.add(gt)
                    gtlst.append(gt)


                #gtlst= [col[i-n_sample] for i in range(n_sample)]
                #for i in range(n_sample):
                    #gt = line.split('\t')[i-n_sample].split(':')[0]
                genotypes[sv_id] = gtlst
                gtlst = []
    print('==> gathering genotype information')
    print('types of gt in your vcf', unique_gt)
    return genotypes
'''
        #for line in infile:
        #    if line.startswith('#CRHOM'):
        #        samples = line.split('FORMAT\t')[-1].split('\t')
                 
        #        print(samples)
        #samples = [x for x in infile if x.startswith('#CHROM')]
        #samples = x.rstrip('\n').split('FORMAT\t')[-1].split('\t') for x in infile if x.startswith('#CHROM')]
        #samples = samples[0].rstrip('\n').split('FORMAT\t')[-1]
        #samples = samples.split('\t')
        #samples = list(map(lambda x: x.rstrip('\n').split('FORMAT\t')[-1], samples[0]))
        #samples = list(map(lambda sample : sample.split('\t'), samples))

        #n_sample = len(samples)
        #print(n_sample)
        #lines = [variant for variant in infile if not variant.startswith('#')]
        #lines = [variant for variant in infile if not variant.startswith('#')]
        #print(lines)
        for line in infile:
            col = line.split('\t')
            print(line)
            sv_id = col[2]
            print(sv_id)
            #gtlst = (col[i-n_sample].split(':')[0] for i in range(n_sample))
            for i in range(n_sample):
                print(col[i-n_sample])
                
            gtlst= [col[i-n_sample] for i in range(n_sample)]
            #print(gtlst)
            #for i in range(n_sample):
                #gt = line.split('\t')[i-n_sample].split(':')[0]
            genotypes[sv_id] = gtlst
            gtlst = ()
            #samples[i] : line.split('\t')[i-n_sample]

    return genotypes
'''

def filter_genotype(vcf, target_class, caller):
    """Return a filtered genotype_dict with the given argument, target_class : germline / somatic / de novo"""
    genotypes = get_genotype(vcf)
    variants = get_variants(vcf, caller)
    #variants = filter_variant_size_quality(vcf, 'PASS', 50, caller)
    filter_genotypes = dict()

    for sv in genotypes:
        sv_filter = variants[sv][1]
        if sv_filter == 'PASS':
            gts = genotypes[sv]
        
            if gts[0] == '0/1' or gts[0] == '1/1':
                if gts[1] == '0/0' and gts[2] == '0/0':# none of the parents have this variant
                    variant_class = 'de_novo'
                else: # at leat one parent has this variant
                    variant_class = 'germline'
            else:# when child has 0/0 gt
                if gts[1] == '0/1' or gts[1] == '1/1' or gts[2] == '0/1' or gts[2] == '1/1':
                    variant_class = 'somatic'
                else:# when all trio has 0/0 - what is this?
                    variant_class = 'NA'
            if variant_class == target_class:
                filter_genotypes[sv] = gts

    print(len(filter_genotypes.keys()), 'are {} variants'.format(target_class))

    return filter_genotypes
    
def genotype_stats(vcf, target_class, caller):
    genotypes = filter_genotype(vcf, target_class, caller)
    variants = filter_variant_size_quality(vcf, 'PASS', 50, caller)
    stats = dict()

    for sv in variants:
        sv_type = variants[sv][4]
        if not sv_type in stats:
            stats[sv_type] = 0
        if sv in genotypes.keys():
            stats[sv_type] += 1
    #print SV type statistics using pprint
    print("*"*50)
    print("dnSVs found in your VCF")
    pp(stats)
    print("*"*50)
    return None


def filter_variant_size_quality(vcf, filter_criteria='PASS', cutoff=50, caller='foo'):
    """Return a filtered variant dictionary containing variants(except BND type) larger than cutoff and variant quality filtered
    in the form of {variant_id : (chromosome, sv_filter_falg, sv_start, sv_end, sv_type)}
    *all BND variants are included"""

    variants = get_variants(vcf, caller)
    filt_variants = variants.copy()

    #qual_variants = [each_variant for each_variant in variants if variants[each_variant][1] == filter_criteria]
    for each_variant in variants:
        info = variants[each_variant]
        filter_flag, sv_start, sv_end, sv_type = info[1], info[2], info[3], info[4]
        if filter_flag == filter_criteria:
            if sv_type == "BND":
                sv_len = 50
            else:
                sv_len = sv_end - sv_start
            # filter out small SVs
            if sv_len < cutoff:
                del filt_variants[each_variant]
        else:
            del filt_variants[each_variant]
        
    return filt_variants

def variant_stats(vcf, cutoff, caller):
    """Return statistics(# of variants by each SV type) from the input variant dictionary"""
    
    variant_dict = filter_variant_size_quality(vcf, 'PASS', cutoff, caller)
    nsv = (len(variant_dict.keys()))
    stat = dict()
    for sv_id in variant_dict:
        sv_type = variant_dict[sv_id][4]
        if not sv_type in stat:
            stat[sv_type] = 0

        stat[sv_type] += 1

    #print SV type statistics using pprint
    print("*"*50)
    print(f"{nsv} SVs found in your VCF")
    pp(stat)
    print("*"*50)
    return stat


def write_filtered_vcf(vcf, vcf_out, caller):
    """Write a filtered vcf file based on filter_variant_size_quality arguments"""
    vcf_filt = filter_variant_size_quality(vcf, 'PASS', 50, caller)
    
    with open(vcf, 'r') as infile, open(vcf_out, 'w') as outfile:
        write_list = [line for line in infile if line.startswith('#') or line.split('\t')[2] in vcf_filt]
        for each in write_list:
            outfile.write(each)
    print("Your VCF is written at '{}'".format(vcf_out))
    return 0

def write_dnSV_vcf(vcf, vcf_out, target_type, caller):
    vcf_filt = filter_variant_size_quality(vcf, 'PASS', 50, caller)
    genotype = filter_genotype(vcf, target_type, caller)

    with open(vcf, 'r') as infile, open(vcf_out, 'w') as outfile:
        write_list = [line for line in infile if line.startswith('#') or (line.split('\t')[2] in vcf_filt and line.split('\t')[2] in genotype)]
        for each in write_list:
            outfile.write(each)
    print("Your VCF is written at '{}'".format(vcf_out))
    return None


def compare_variants(vcf1, vcf2):
    variants_base = get_variants(vcf1)
    variants_base_unique = variants_base.copy()
    variants_compare = get_variants(vcf2)

    common_variants = 0
    for each in variants_base:
        if variants_base[each] in variants_compare.values():
            del variants_base_unique[each]
            common_variants += 1
    print(common_variants, "common variants")
    return variants_base_unique



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for processing VCF and calculating performance ")
    parser.add_argument("-i", "--input_file", help = " input VCF file ", required = True)
    parser.add_argument("-o", "--output_file", help = " output variant vcf ", required = False)
    parser.add_argument("-d", "--dnsv_output", help = " output de novo variant vcf ", required = False)
    parser.add_argument("-c", "--length_cutoff", help = " SV length cutoff(>=50) ", type = int, default = 50, required = False)
    parser.add_argument("-g", "--genotype", help = " target genotype (germline, somtaic, de_novo ", required = False)
    parser.add_argument("-s", "--caller", help = " name of SV caller : etching/delly/lumpy/gridss/svaba ", required = False, default = 'foo')
    args = parser.parse_args()

    import sys, os, re
    from pprint import pprint as pp

    variant_stats(args.input_file, args.length_cutoff, args.caller)
    if args.output_file:
        write_filtered_vcf(args.input_file, args.output_file, args.caller)

    if args.genotype:
        filter_genotype(args.input_file, args.genotype, args.caller)
        genotype_stats(args.input_file, args.genotype, args.caller)
    if args.dnsv_output:
        write_dnSV_vcf(args.input_file, args.dnsv_output, args.genotype, args.caller)
