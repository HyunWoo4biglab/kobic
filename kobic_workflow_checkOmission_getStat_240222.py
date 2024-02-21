

def get_meta_data(metafile):
    """A function to get sample, family, path information from meta file (ped file)

    Assuming he metafile composed as
    sampleID\tfailyID\tversion\trelationship\tdatapath
    """

    path_dict = dict()
    with open(metafile) as infile:
        lines = [line for line in infile if not line.startswith('paternal')]
        for i, line in enumerate(lines):
            col = line.strip().split()
            #family_num = i // 2
            family_num = re.search('/(\d+)\.Family_(\d+)/', col[1]).group(2)
            family_id = f'family_{family_num}'
            ####
            #pat_id, pat_path, mat_id, mat_path, child_id, child_path, child_sex = col ### change this with the given metadata format
            ####
            if not family_id in path_dict:
                path_dict[family_id] = dict()
            for j in range(0, len(col)-1, 2):
                relation = j // 2 #0 : father, 1: mother, 2 : child
                sample_id = col[j]
                path = col[j+1]
                if not sample_id in path_dict[family_id]:
                    path_dict[family_id][sample_id] = (relation, path)

    return path_dict

def get_files_by_id(sampleid, path):
    """
    | A function to get target fastq files under a given path by its sample_id
    """
    sample_prefix = path.split('/')[-1]
    files = os.listdir(path)
    files = list(filter(lambda x: re.search(sampleid, x), files))
    files = list(filter(lambda x: x.endswith('fastq.gz') or x.endswith('fq.gz'), files))
    list(map(lambda x: path + x, files))

    return files


def check_omit_make_command(path_dict, output, result, threads=30):
    """
    | A function to check whether KMC, KMC union, ETCHING has run properly. If not, rerun.
    """

    if not output.endswith('/'):
        output += '/'
    outputfile = output + f'check_omit_rerun_etching.sh'
    working_dir = output + 'svcall/'
    

    sickle = config.sickle_path
    #samtools = config.samtools_path
    kmc = config.kmc_path
    kmc_tools = config.kmc_tools_path
    etching = config.etching_path
    genome = config.reference_genome
    pgk = config.pgk_path

    family_list = sorted(path_dict.keys())
    print(family_list)
    with open(outputfile, 'w') as outfile:
        outfile.write(f'printf "Sample\\tTotal\\tDEL\\tDUP\\tINV\\tBND\\n" >> {result}\n')
        for family in family_list:

            family_dir = working_dir + family + '/'
            if not os.path.exists(family_dir):
                os.makedirs(family_dir)

            etching_dir = family_dir + 'etching_call/'
            if not os.path.exists(etching_dir):
                os.makedirs(etching_dir)

            trim_dir = family_dir + '1.preprocess/'
            if not os.path.exists(trim_dir):
                os.makedirs(trim_dir)
            parentfile = family_dir + f'{family}_sample.list'
            parent_db = f'{family}_parent'
            pgkp_db = family_dir + f'{family}_pgkp'
            parentfile_list = list(); etching_cmd_d = dict()

            for j, sample in enumerate(path_dict[family]):
                relation, path = path_dict[family][sample]
                if not path.endswith('/'):
                    path += '/'
                files = get_files_by_id(sample, path)
                files.sort()
                trim_files = list(map(lambda x: trim_dir + re.sub('f.*q.gz', 'trim.fastq', x), files))
                single_output = re.sub('1.trim.fastq', 'single.trim.fastq', trim_files[0])
                trim_files.append(single_output)
                files = list(map(lambda x : path + x, files))
                sickle_cmd = f'{sickle} pe -f {files[0]} -r {files[-1]} -o {trim_files[0]} -p {trim_files[1]} -s {trim_files[2]} -t sanger -q 20 -l 70'
                if (not os.path.exists(trim_dir)) or (not os.path.exists(trim_files[0])):
                    outfile.write(f'{sickle_cmd}\n')
                if relation < 2:
                    parentfile_list.extend(trim_files)
                else:
                    etching_prefix = f'{family}_sample_{sample}'
                    etching_check_file = f'{etching_dir}{etching_prefix}.scored.filtered.typed.vcf'
                    etching_output_dir = etching_dir + etching_prefix + '/'
                    if not os.path.exists(etching_output_dir):
                        os.makedirs(etching_output_dir)

                    etching_cmd = f'{etching} -1 {trim_files[0]} -2 {trim_files[1]} -g {genome} -f {pgkp_db} -t {threads} --keep-kmc -o {etching_prefix} --output-dir {etching_output_dir} --work-dir {etching_output_dir}'

                    etching_cmd_d[sample] = {'exist' : os.path.exists(etching_check_file), 'cmd' : etching_cmd}

                if j == len(path_dict[family]) -1:
                    if not os.path.exists(parentfile):
                        parentfile_list = list(filter(lambda x : not 'single.trim' in x, parentfile_list))
                        pfile = open(parentfile, 'w')
                        for p in parentfile_list:
                            pfile.write(f'{p}\n')
                        pfile.close()

                    if not os.path.exists(f'{pgkp_db}.kmc_suf'):
                        kmc_cmd = f'{kmc} -v -k31 -m12 -cx10000 -ci5 -t{threads} -fq @{parentfile} {family_dir}{parent_db} {family_dir}'
                        outfile.write(f'{kmc_cmd}\n')
                        parent_db = family_dir + parent_db
                        kmc_union_cmd = f'{kmc_tools} simple {parent_db} -ci5 {pgk} union {pgkp_db}'
                        outfile.write(f'{kmc_union_cmd}\n')

                    for s in etching_cmd_d:
                        if etching_cmd_d[s]['exist']:
                            get_result(f'{etching_check_file}', sample, outfile, result)
                            outfile.write(f'mv {etching_dir}{etching_prefix}.* {etching_output_dir}\n')
                        else:
                            outfile.write(f'{etching_cmd_d[s]["cmd"]}\n')



#            else:
#                etching_dir = family_dir + 'etching_call/'
#                for j, sample in enumerate(path_dict[family]):
#                if not os.path.exists(etching_dir):
#
#                else:
#                    for j, sample in enumerate(path_dict[family]):
#                        relation, path = path_dict[family][sample]
#                        if relation >= 2:
#                            etching_prefix = f'{family}_sample_{sample}'
#                            etching_check_file = f'{etching_prefix}.scored.filtered.typed.vcf'
#                            etching_output_dir = etching_dir + etching_prefix + '/'
#                            if not os.path.exists(etching_output_dir):
#                                os.makedirs(etching_output_dir)
#
#                            if os.path.exists(etching_check_file):
#                                get_result(f'{etching_dir}{etching_check_file}', result)
#                                outfile.write(f'mv {etching_dir}{etching_prefix}.* {etching_output_dir}\n')
#                            else:
#                                etching_cmd = f'{etching} -1 {trim_files[0]} -2 {trim_files[1]} -g {genome} -f {pgkp_db} -t {threads} --keep-kmc -o {etching_prefix} --output-dir {etching_output_dir} --work-dir {etching_output_dir}'
#                                outfile.write(f'{etching_cmd\n')
#
#                if os.path.exists(etching_dir):
#                    if os.path.exists(
#                    files = os.listdir(etching_dir)
#                    files = list(filter(lambda x: x.endswith('.scored.filtered.typed.vcf'), files))
#                    if len(files) > 0:
#                        
#                
#
#
#    outfile.write(f'printf "Total\tDEL\tDUP\tINV\tBND\n" >> {result}')
#    kmc_p_path = f'{family_dir}/{family}_parent.kmc_suf'
#    kmc_pbkp_path = f'{family_dir}/{family}_pgkp.kmc_suf'
#    if (os.path.exists(kmc_p_path) and os.path.exists(kmc_pgkp_path)):
#        pass
#    else:
#        outfile.write(kmc, kmc_union)
#
#    etching_path = f'{etching_dir}/{etching_prefix}.scored.filtered.typed.vcf'
#    if os.path.exists(path):
#        get_result(etching_path, outputresult)
#
#    else:
        

def get_result(inputfile, sample, outfile, result):

    outfile.write(f'total=$(grep -v \'^#\' {inputfile} -c)\n')
    outfile.write(f'del=$(grep -v \'^#\' {inputfile} | grep \'SVTYPE=DEL\' -c)\n')
    outfile.write(f'dup=$(grep -v \'^#\' {inputfile} | grep \'SVTYPE=DUP\' -c)\n')
    outfile.write(f'inv=$(grep -v \'^#\' {inputfile} | grep \'SVTYPE=INV\' -c)\n')
    outfile.write(f'bnd=$(grep -v \'^#\' {inputfile} | grep \'SVTYPE=BND\' -c)\n')
    outfile.write(f'printf "{sample}\\t')
    outfile.write('${total}\\t${del}\\t${dup}\\t${inv}\\t${bnd}\\n" >> ')
    outfile.write(f'{result}\n')

    return None

if __name__ == '__main__':

    import os, sys, re, argparse, datetime
    import config

    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-i', '--input_ped', help = 'input pedigree fil in .ped format', required = True)
    parser.add_argument('-o', '--output', help = 'output directory', required = False)
    #parser.add_argument('-n', '--ngroup', help = 'number of groups to split the intervals', type = int, required = False, default = 3)
    parser.add_argument('-r', '--sv_result', help = 'output result for sv stats', required = True)
    parser.add_argument('-t', '--num_threads', help = 'number of threads to run programs', type = int, required = False, default = 5)
    args = parser.parse_args()


    path_dict = get_meta_data(args.input_ped)
    check_omit_make_command(path_dict, args.output, args.sv_result, args.num_threads)
