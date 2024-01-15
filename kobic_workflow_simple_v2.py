"""A script to generate deNovo-ETCHING shell script on KOBIC rare-disease cohort"""
# 11-24-2023
# Written by Hyun Woo Kim

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
            family_num = i // 2
            family_id = f'family{family_num}'
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


def get_files(path, trimbool):
    """A function to get the target fastq files present under the given path"""

    #path = '/'.join(path.split('/')[:-1])
    sample_prefix = path.split('/')[-1]
    files = os.listdir(path)
    files = list(filter(lambda x: x.endswith('fastq.gz') or x.endswith('fq.gz'), files))
    if trimbool == 'trim':
        files = list(filter(lambda x: x.endswith('trim.fastq.ga'), files))
    list(map(lambda x: path + x, files))

    return files

def get_files_by_id(sampleid, path):
    """A function to get target fastq files under a given path by its sample_id"""
    sample_prefix = path.split('/')[-1]
    files = os.listdir(path)
    files = list(filter(lambda x: re.search(sampleid, x), files))
    files = list(filter(lambda x: x.endswith('fastq.gz') or x.endswith('fq.gz'), files))
    list(map(lambda x: path + x, files))

    return files

def make_command(path_dict, interval, output, threads=30):
    """A function to write the commands for each sample into a shell script with a given interval"""

    if not output.endswith('/'):
        output += '/'
    start, end = interval
    family_list = sorted(path_dict.keys())
    outputfile = output + f'run_family{start}_{end}.sh'
    working_dir = output + 'svcall/'

    sickle = config.sickle_path
    #samtools = config.samtools_path
    kmc = config.kmc_path
    kmc_tools = config.kmc_tools_path
    etching = config.etching_path
    genome = config.reference_genome
    pgk = config.pgk_path

    with open(outputfile, 'w') as outfile:
        for family in family_list[start : end]:
            family_dir = working_dir + family + '/'
            if not os.path.exists(family_dir):
                os.makedirs(family_dir)
            etching_dir = family_dir + 'etching_call/'
            if not os.path.exists(etching_dir):
                os.makedirs(etching_dir)
 
            parentfile = family_dir + f'{family}_sample.list'
            parent_db = f'{family}_parent'
            pgkp_db = family_dir + f'{family}_pgkp'
            parentfile_list = list(); etching_cmd_list = list()
            trim_dir = family_dir + '1.preprocess/'
            if not os.path.exists(trim_dir):
                os.makedirs(trim_dir)
            for j, sample in enumerate(path_dict[family]):
                relation, path = path_dict[family][sample]
                if not path.endswith('/'):
                    path += '/'
                #files = get_files(path, 'untrim')
                files = get_files_by_id(sample, path)
                files.sort()
                trim_files = list(map(lambda x: trim_dir + re.sub('f.*q.gz', 'trim.fastq.gz', x), files))
                #r1_output = re.sub('f.*q.gz', 'trim.fastq.gz', files[0])
                #r2_output = re.sub('f.*q.gz', 'trim.fastq.gz', files[1])
                single_output = re.sub('1.trim.fastq.gz', 'single.trim.fastq.gz', trim_files[0])
                trim_files.append(single_output)
                files = list(map(lambda x : path + x, files))
                print(files)
                print(trim_files)
                #sickle_cmd = f'{sickle} pe -f {files[0]} -r {files[-1]} -o {sickle_dir + r1_output} -p {sickle_dir + r2_output} -s {sickle_dir + single_output} -t sanger -q 20 -l 70 -g'
                sickle_cmd = f'{sickle} pe -f {files[0]} -r {files[-1]} -o {trim_files[0]} -p {trim_files[1]} -s {trim_files[2]} -t sanger -q 20 -l 70'
                outfile.write(f'{sickle_cmd}\n')

                #files = get_files(output, 'trim')
                #files = list(map(lambda x : sickle_output + x, files))
                if relation < 2:
                    parentfile_list.extend(trim_files)
                else:
                    etching_prefix = f'{family}_sample_{sample}'
                    etching_cmd = f'{etching} -1 {trim_files[0]} -2 {trim_files[1]} -g {genome} -f {pgkp_db} -t {threads} --keep-kmc -o {etching_prefix} --output-dir {etching_dir} --work-dir {etching_dir}'
                    etching_cmd_list.append(etching_cmd)
                if j == len(path_dict[family]) -1: # the last person in the quartet
                    parentfile_list = list(filter(lambda x : not 'single.trim' in x, parentfile_list))
                    print(parentfile_list)
                    pfile = open(parentfile, 'w')
                    for p in parentfile_list:
                        pfile.write(f'{p}\n')
                    pfile.close()

                    kmc_cmd = f'{kmc} -v -k31 -m12 -cx10000 -ci5 -t{threads} -fq @{parentfile} {family_dir}{parent_db} {family_dir}'
                    outfile.write(f'{kmc_cmd}\n')
                    parent_db = family_dir + parent_db
                    kmc_union_cmd = f'{kmc_tools} simple {parent_db} -ci5 {pgk} union {pgkp_db}'
                    outfile.write(f'{kmc_union_cmd}\n')
                    for cmd in etching_cmd_list:
                        outfile.write(f'{cmd}\n')

# Main execution
if __name__ == '__main__':

    import os, sys, re, argparse, datetime
    import config

    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-i', '--input_ped', help = 'input pedigree fil in .ped format', required = True)
    parser.add_argument('-o', '--output', help = 'output directory', required = False)
    parser.add_argument('-n', '--ngroup', help = 'number of groups to split the intervals', type = int, required = False, default = 3)
    parser.add_argument('-t', '--num_threads', help = 'number of threads to run programs', type = int, required = False, default = 5)
    args = parser.parse_args()

    path_dict = get_meta_data(args.input_ped)
    print(path_dict)

    interval_len = len(path_dict) / args.ngroup
    for i in range(args.ngroup):
        interval = (int(i*interval_len), int((i+1)*interval_len))
        make_command(path_dict, interval, args.output, args.num_threads)
