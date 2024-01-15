"""A script to run deNovo-ETCHING automatically on KOBIC rare-disease cohort"""
# 31st July 2023
# Written by Hyun Woo Kim



def getPedigree(ped):
    """ A function to get pedigree information for target trios
        Assuming ped file is composed as the following format

        FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\n
    """

    ped_info = dict(); path_info = dict()
    with open(ped, 'r') as infile:
        lines = [line for line in infile if not line.startswith('family') and not line.startswith('Family')] # in case ped file contains header
        for line in lines:
            #col = line.strip().split('\t')
            family_id, sample_id, paternal_id, maternal_id, sex, phenotype, sample_path, paternal_path, maternal_path = line.strip().split('\t')
            if  paternal_id != 0 and maternal_id != 0:
                if not sample_id in ped_info:
                    ped_info[sample_id] = dict()
                ped_info[sample_id]['paternal'] = paternal_id
                ped_info[sample_id]['maternal'] = maternal_id
                ped_info[sample_id]['familyID'] = family_id
                if ';' in sample_path:
                    
                    sample_path, sample_long_path = sample_path.split(';')
                    paternal_path, paternal_long_path = paternal_path.split(';')
                    maternal_path, maternal_long_path = maternal_path.split(';')

                if not sample_path.endswith('/'):
                    sample_path += '/'
                if not paternal_path.endswith('/'):
                    paternal_path += '/'
                if not maternal_path.endswith('/'):
                    maternal_path += '/'

                path_info[sample_id] = sample_path.strip()
                path_info[paternal_id] = paternal_path.strip()
                path_info[maternal_id] = maternal_path.strip()

    return ped_info, path_info

def get_sample_path(path_file):
    """
    A function to get file path for each target sample(trio) from an argument file
    Assuming path file is composed as the following format

    IndividualID\tfile_path(or directory_path)\n
    """
    path_info = dict()
    with open(path_file, 'r') as infile:
        lines = [line for line in infile if not line.startswith('sample') and not line.startswith('Sample')]
        for line in lines:
            sample_id, sample_path = line.strip().split('\t')
            if ';' in sample_path:
                sample_path = sample_path.split(';')[0]
                longread_path = sample_path.split(';')[-1][1:]
            path_info[sample_id] = sample_path

    return path_info


def get_fastq_paths(path_info):
    """
    A script to get target fastq file path for each individual given a sample's directory
    Assuming the file is stored as the following architecture

    ./sample/R1.fastq
    ./sample/R2.fastq

    This function assumes the path argument 'path_info' is given as './sample/'
    """

    if not path_info.endswith('/'):
        path_info += '/'
    files = os.listdir(path_info)
    #files = list(filter(lambda x : 'R1' in x or 'R2' in x and x.endswith('.fastq'), files))
    files = list(filter(lambda x : x.endswith('.fq.gz') or x.endswith('.fastq.gz'), files))
    #files = list(filter(lambda x : x.endswith('.fastq'), files))
    files.sort()

    assert len(files) == 2, "target fastq files are more than 2!"

    return files

def make_kmc_cmd(family_id, parental_files, output_dir, k=31, m=10, cx=10000, ci=2, n_thread=10):
    ref_genome = config.reference_genome
    kmc_path = config.kmc_path
    if not output_dir.endswith('/'):
        output_dir += '/'
    output_dir += 'kmer_db/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    sample_file = output_dir + 'sample.list'
    output_prefix = output_dir + f'family{family_id}_parents'
    with open(sample_file, 'w') as outfile:
        for parent in parental_files:
            outfile.write(f'{parent}\n')
    kmc_cmd = f'{kmc_path} -v -k31 -m{m} -cx{cx} -ci{ci} -t{n_thread} -fq @{sample_file} {output_prefix} {output_dir}'
    return kmc_cmd


def make_kmc_union_cmd(family_id, output_path, ci=5):
    kmc_tools_path = config.kmc_tools_path
    pgk = config.pgk_path

    if not output_path.endswith('/'):
        output_path += '/'
    output_path += 'kmer_db/'
    output_db = output_path + f'family{family_id}_pgkp'
    parentk = f'{output_path}family{family_id}_parents'

    kmc_union_cmd = f'{kmc_tools_path} simple {parentk} -ci{ci} {pgk} union {output_db}'
    return kmc_union_cmd




def make_etching_command(family_id, target_dir, sample_files, output_path, n_thread):
    etching_path = config.etching_path
    ref_genome = config.reference_genome

    if not output_path.endswith('/'):
        output_path += '/'
    kmer_db = output_path + f'/kmer_db/family{family_id}_pgkp'
    output_dir = output_path + 'etching_variant_call/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    r1 = sample_files[0]
    r2 = sample_files[1]

    output_prefix = f'family{family_id}_' + r1.split('/')[-1].split('_1.fq')[0]
#    etching_cmd = '{program} -1 {r1} -2 {r2} -g {genome} -f {db} -t {threads} -o {output} --output-dir {output_dir} --work-dir {output}'.format(program=etching_path, r1=r1, r2=r2, genome=ref_genome, db=kmer_db, threads=n_thread, output=output_prefix, output_dir=output_dir)
    etching_cmd = '{program} -1 {r1} -2 {r2} -g {genome} -f {db} -t {threads} --keep-kmc -o {output} --output-dir {output_dir} --work-dir {working_dir}'.format(program=etching_path, r1=target_dir+sample_files[0], r2=target_dir+sample_files[1], genome=ref_genome, db=kmer_db, threads=n_thread, output=output_prefix, output_dir=output_dir, working_dir=output_dir)
    return etching_cmd


def make_cmds_on_trios(ped_info, path_info, output_dir):
    """
    A function to run KMC - KMC_TOOLS - ETCHING in a linear manner for each target trio
    
    """
    families_cmd = dict(); run_etching_cmd = dict()
    for sample in ped_info:
        cmds = dict()

        family_id = ped_info[sample]['familyID']

        if family_id in families_cmd:
            continue

        pat = ped_info[sample]['paternal']
        mat = ped_info[sample]['maternal']

        if not output_dir.endswith('/'):
            output_dir += '/'
        output_path = output_dir + f'family{family_id}/'

        # get directory path
        sample_dir = path_info[sample]
        pat_dir = path_info[pat]
        mat_dir = path_info[mat]

        # get fastq file paths
        sample_files = get_fastq_paths(sample_dir)
        pat_files = get_fastq_paths(pat_dir)
        pat_files = list(map(lambda x: pat_dir + x, pat_files))
        mat_files = get_fastq_paths(mat_dir)
        mat_files = list(map(lambda x: mat_dir + x, mat_files))
        pat_files.extend(mat_files)
        parental_files = pat_files

        kmc_cmd = make_kmc_cmd(family_id, parental_files, output_path, 31, 8, 10000, 2, 10)
        kmc_union_cmd = make_kmc_union_cmd(family_id, output_path, 5)

        etching_cmd = make_etching_command(family_id, sample_dir, sample_files, output_path, 10)

        cmds['kmc_cmd'] = kmc_cmd
        cmds['kmc_union_cmd'] = kmc_union_cmd
        #cmds['etching_cmd'] = etching_cmd

        families_cmd[f'family{family_id}'] = cmds
        run_etching_cmd[sample_id] = etching_cmd
        #run_programs_linearly(family, cmds, running_process)
        #run_programs_linearly(family, kmc_cmd, kmc_union_cmd, etching_cmd, running_process)

    return families_cmd, run_etching_cmd


#def run_multiple_trios(families_cmd, n_process):

#    families = list(families_cmd.keys())
#    #for i in range(max(int(len(families)/n_process), 1)):
#        
#    for family in families_cmd:
#        
#        run_programs_linearly(families_cmd[family])
#    
#commands = ['cmd1', 'cmd2', 'cmd3', 'cmd4', 'cmd5'] 
#n = 2 #the number of parallel processes you want
#for j in range(max(int(len(commands)/n), 1)):
#    procs = [subprocess.Popen(i, shell=True) for i in commands[j*n: min((j+1)*n, len(commands))] ]
#    for p in procs:
#        p.wait()

        #result = pool.apply_async(run_programs_linearly, args = (family, kmc_cmd, kmc_union_cmd, etching_cmd, running_process))
        #results.append(result)

        #while len(runing_process.keys()) > 10e
        #    time.sleep(600)
    #for result in results:
    #    result.get()

#    print('All trios are completed!')
#
#    return None

def run_programs_linearly(family, kmc_cmd, kmc_union_cmd, etching_cmd, running_process):
        # Run each programs lienarly

        #process1 = subprocess.Popen(kmc_cmd, shell=True)
        #running_process[family] = process1
        #process1.wait()
        ## RUN KMC
        trial = 1
        parent_db = kmc_cmd.split(' ')[-2] + '.kmc_suf'
        pgkp_db = kmc_union_cmd.split(' ')[-1] + '.kmc_suf'

        print('here is the db', parent_db)
        print('here is the pgkp db', pgkp_db)

        if len(running_process) > 0:
            #if family in run.args:
            time.sleep(10)
        elif len(running_process) > 1:
            time.sleep(30)

        for pid, run in running_process.items():
            while re.search(f"(kmc)(.*)({family})(.*)", run.args):
                print(f'same family {family} in process..')
                time.sleep(300)

        if not os.path.exists(parent_db):
            process1 = run_monitor_process(kmc_cmd, family, trial, running_process)
        else:
            process1 = True

        if process1:
            for pid, run in running_process.items():
                while re.search(f"(kmc)(.*)({family})(.*)", run):
                    time.sleep(300)
            ## RUN KMC_TOOLS UNION if process1 (run kmc) is completed without any error
            trial = 1
            if not os.path.exists(pgkp_db):
                process2 = run_monitor_process(kmc_union_cmd, family, trial, running_process)
            else:
                process2 = True

            if process2:
                ## RUN ETCHING if process2 (run kmc_tools - union) is completed without any error
                trial = 1
                process3 = run_monitor_process(etching_cmd, family, trial, running_process)

        if process3:
            print(f'Family number {family} analysis is comleted!')


def make_pgkp(family, kmc_cmd, kmc_union_cmd, running_process):
    
    pass


def run_etching_w_pgkp(sample, etching_cmd, running_process):
    pass
        
def run_monitor_process(cmd, sample, trial, running_process):
    """
    A function to run a given command 'cmd' an check for the return_code
    If the return_code insists an error, it re-tries the command 4 times more.
    If the error occurs more than 5 times, it quits the job and return an error messge
    """
    if trial <= 5:
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
        running_process[process.pid] = process
        print(process.args)
        print(running_process)
        process.wait()

        while process.poll() is None:
            print(f'Process {process.pid} not complete yet, going to sleep until it\'s done')
            time.sleep(300)
        if process.poll() is not None:
            exit_code = process.returncode
            if exit_code == 0:
                del running_process[process.pid]
                return True
            else:
                # re-run the command
                trial += 1
                run_monitor_process(cmd, sample, trial, running_process)
    else:
        print(f'skipping {sample} due to over 5 times of failure')
        return False


# Main execution
if __name__ == '__main__':

    import os, sys, argparse, datetime
    import subprocess
    import multiprocessing
    #import resources
    import psutil
    import time
    import config
    import re

    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-p', '--input_ped', help = 'input pedigree fil in .ped format', required = True)
    #parser.add_argument('-a', '--input_arg', help = 'input argument file containing sample path', required = True)
    parser.add_argument('-o', '--output', help = 'output directory', required = False, default = 'kobic_rare_etching')
    #parser.add_argument('-k', '--kmer_db', help = 'PGK kmer database', required = True)
    parser.add_argument('-n', '--num_process', help = 'number of process to run parallely', type = int, required = False, default = 5)
    args = parser.parse_args()


    #family_list = [(kmc_cmd, kmc_union_cmd, etching_cmd), (), (), ()]


    ped_info, path_info = getPedigree(args.input_ped)
    family_dict = make_cmds_on_trios(ped_info, path_info, args.output)
    target_list = ['family61', 'family62', 'family79', 'family--']
    running_process = dict()
    manager = multiprocessing.Manager()
    results_d = manager.dict()
    with multiprocessing.Pool(processes = args.num_process) as pool:
        #manager = multiprocessing.Manager()
        #results_dict = manager.dict()
        runs = [pool.apply_async(run_programs_linearly, args=(family, family_dict[family]['kmc_cmd'], family_dict[family]['kmc_union_cmd'], family_dict[family]['etching_cmd'], running_process), callback = lambda res, sid : results_d.update({sid : res}), error_callback = print) for family in family_dict if family in target_list]
#        runs = [pool.apply_async(run_programs_linearly, args=(family, family_dict[family]['kmc_cmd'], family_dict[family]['kmc_union_cmd'], family_dict[family]['etching_cmd'], running_process), callback = print, error_callback = print) for family in family_dict if family in target_list]

        for job in runs:
            job.get()

    for r in results_d:
        print(f'{r}\n{results_d[r]}\n')
    #pool.close()
    #pool.join()


    #pool.apply_async(run_programs_on_trios, args = (args.input_ped, args.input_path, args.output))
    #pool.close()
    #pool.join()

