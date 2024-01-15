"""A script to run ETCHING using python subprocess"""
# 07.21.2023
## added multiprocessing and job monitoring for automated running


def get_samples(target_dir, target_list):
    samples = os.listdir(target_dir)
    samples.sort()
    samples = list(filter(lambda sample : sample in target_list), samples)



    samples = list(filter(lambda files : files.endswith('.fastq'), samples))
    samples = list(filter(lambda files : 'R1' in files or 'R2' in files, samples))

    samples = list(map(lambda f : target_dir + f, samples))

    return samples

def get_fastq_files(target_dir, ):
    samples = os.listdir(target_dir)


def run_etching(input_dir, ):

    import config as cf

    target_dir = input_dir
    samples = get_samples(target_dir)
    assert len(samples) == 2, "target fastq samples are more than 2!"

    etching_path = cf.etching_path
    kmer_db = args.kmer_db
    n_thread = args.n_thread
    output_prefix = args.output
    #output_dir = args.output_dir
    ref_genome = args.ref_genome

    etching_cmd = '{program} -1 {fwd} -2 {rev} -g {genome} -f {db} -t {threads} -o {output}'.format(program=etching_path, fwd=target_dir+samples[0], rev=target_dir+samples[1], genome=ref_genome, db=kmer_db, threads=n_thread, output=output_prefix)

    print('-'*10, 'your command')
    print(etching_cmd)
    process = Popen(etching_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    process.wait()


def make_kmc_cmd(ref_genome, parents, n_thread, output_dir, k=31, m=10, cx=10000, ci=2):
    kmc_path = config.kmc_path
    output_dir += '/kmer_db/'
    if not os.path.exists(output_dir):
        os.makedir(output_dir)
    sample_file = output_dir + 'sample.list'
    output_prefix = output_dir + 'parents'
    with open(sample_file, 'w') as outfile:
        for parent in parents:
            outfile.write(f{'parent}\n')
    kmc_cmd = f'{kmc_path} -v -k31 -m{m} -cx{cx} -ci{ci} -t{n_thread} -fa @{sample_file} {output_prefix} {output_dir}'

    return kmc_cmd

def make_kmc_union_cmd(pgk, parentk, output_db, ci=5):
    kmc_tools_path = config.kmc_tools_path
    #output_db = output_dir + 'pgkp'
    kmc_union_cmd = f'{kmc_tools_path} simple {parentk} -ci{ci} {pgk} union {output_db}'

    return kmc_union_cmd


def make_etching_command(target_dir, kmer_db, ref_genome, n_thread, output_dir):
    etching_path = config.etching_path
    files = get_samples(target_dir)
    files.sort()
    assert len(files) == 2, "target fastq files are more than 2!"

    output_prefix = fwd.split('.fastq')[0]

    etching_cmd = '{program_path} -1 {r1} -2 {r2} -g {genome} -f {kmer_db} -t {threads} -o {output}'.format(program=etching_path, r1=target_dir+files[0], r2=target_dir+files[1], genome=ref_genome, db=kmer_db, threads=n_thread, output=output_prefix)

    return etching_cmd

def run_programs_linearly(trio_dir, output_dir):
    """Run kmc - kmc_tools union - ETCHING in a linear manner(step by step on one trio at a time)"""
    trios = sys.listdir(trio_dir)

    for target_trio in trios:
        parent = list(); child = list()
        samples = sys.listdir(trio_dir + target_trio)
        for sample in samples:
            if "Father" in sample or "Mother" in sample:
                parents.append(f'{trio_dir}{target_trio}/{sample}/')
            else:
                child.append(f'{trio_dir}{target_trio}/{sample}/')
        target_output_dir = output_dir + target_trio + '/'

        # RUN KMC with parent files
        kmc_cmd = make_kmc_cmd(ref_genome, parents, n_trhead, output_dir, 31, 10, 10000, 2)
        process1 = subprocess.Popen(kmc_cmd, shell=True)
        process1.wait()
        #process1 = run_process_check_output(kmc_cmd)

        # RUN KMC_TOOLS to unionize parent kmerDB with PGK
        while process1.poll() is None:
            print('Process not completed yet, going into sleep...')
            time.sleep(300)
        if process.poll() is not None:
            exit_code = process.returncode

        if exit_code = 0:#The process has completed
            del running_processes[data]
            
        #if process.poll() is not None:
        parentk = target_output_dir + '/kmer_db/parents'
        pgkp = target_output_dir + '/kmer_db/pgkp'
        kmc_union_cmd = make_kmc_union_cmd(pgk, parentk, pgkp, 5)
        process2 = subprocess.Popen(kmc_union_cmd, shell=True)
        process2.wait()
    
        # RUN ETCHING with PGKP
        while process2.poll() is None:
            time.sleep(300)
        #if process2.poll() is not None:

        etching_cmd = make_etching_cmd(target_dir, kmer_db, ref_genome, n_thread, output_dir)
        process3 = subprocess.Popen(etching_cmd)
        process3.wait()

        if process3.poll() is not None:
            if os.path.is_file(output_vcf)
                print(f"Trio {target_trio} is done!")


def run_program_parallely(trio_dir, output_dir):
    """Run each program(kmc, kmc_tools union, ETCHING) on multiple trios in a parallel manenr"""

    process = subprocess.Popen(command, shell=True)
    running_processes[data] = process



# Function to monitor running processes
def monitor_process_returncode(max_process):
    while True:
    while len(running_processes.keys()) >= max_process:
        for data, process in list(running_processes.items()):
            if process.poll() is not None:
                exit_code = process.returncode

                if exit_code = 0:
                    # The process has completed
                    del running_processes[data]
                    print(f"Your process for data number {data} has completed successfully.")
                else:
                    print(f"Your process for data number {data} has encountered an error (Exit code: {exit_code}).")
                    # check if the desired final output is generated. If not,

                    output_file_path = f"data_{data}/{final_output}"
                    if not os.path.exists(output_file_path):
                        print(f"Rerunning Program A for data number {data}...")
                        run_program_A(data)
                    else:
                        print(f"Desired output not achieved, but 'final_output.txt' exists for data number {data}.")
                
                del running_processes[data]
            
    
        time.sleep(600)  # Adjust the sleep duration as per your requirement


def run_process_check_output(cmd):
    """Run program and monitor its error. If so, check the desired final output and rerun the program recursively"""
    try:
        process = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        running_processes[] = process
        print(f"Program A for data number {data} has completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Program A for data number {data} has encountered an error (Exit code: {e.returncode}).")

        # Rerun the process if desired output not achieved
        output_file_path = f"data_{data}/{final_output}"
        if not os.path.exists(output_file_path):
            print(f"Rerunning Program A for data number {data}...")
            run_process_check_output(cmd)
        else:
            print(f"Desired output not achieved, but '{final_output}' exists for data number {data}.")


"""
    while True:
        for data, process in list(running_processes.items()):
            if process.poll() is not None:
                exit_code = process.returncode

                if exit_code = 0:
                    # The process has completed
                    del running_processes[data]
                    print(f"Your process for data number {data} has completed successfully.")
                else:
                    print(f"Your process for data number {data} has encountered an error (Exit code: {exit_code}).")
                    # check if the desired final output is generated. If not,

                    output_file_path = f"data_{data}/{final_output}"
                    if not os.path.exists(output_file_path):
                        print(f"Rerunning Program A for data number {data}...")
                        run_program_A(data)
                    else:
                        print(f"Desired output not achieved, but 'final_output.txt' exists for data number {data}.")

                del running_processes[data]


        time.sleep(600)  # Adjust the sleep duration as per your requirement
"""


def check_my_resources(running_processes):
    """A function to check system resources for your processes"""

    my_process_ids = [process.pid for process in list(running_processes.values())]

    # Calculate total CPU and memory usage for your processes
    my_cpu_percent = sum(psutil.Process(pid).cpu_percent(interval=None) for pid in my_process_ids)
    my_memory_percent = sum(psutil.Process(pid).memory_percent() for pid in my_process_ids)

    # Get the number of your running processes
    job_count = len(running_processes)

    # Modify the conditions as per your requirements
    if my_cpu_percent < 80 and my_memory_percent < 80 and job_count < num_processes:
        return True
    else:
        return False

def check_resources():
    cpu_percent = psultil.cpu_percent()
    memory_percent = psultil.virtual_memory().percent

     # Modify the conditions as per your requirements
    if cpu_percent < 80 and memory_percent < 80 and job_count < num_processes:
        return True
    else:
        return False
    
    
# Main execution
if __name__ == '__main__':

    import os, argparse, datetime
    from subprocess import *
    import multiprocessing
    import resources
    import psutil
    import time


    # get target trio samples 
    ## This assumes samples are stored in each trio directory
    target_samples = get_samples(target_dir)

    # run KMC- KMC_TOOLS UNION - ETCHING 
    for i in range(len(target_samples)):
        run_programs_linearly(i, output_dir)
        running_process[i] = 
        while check_running_jobs():
            date.sleep(600)

        #catch = target_samples[i:i+10]
        #for j in catch:
        #    run_programs_linearly(j, output_dir)

        run_programs_linearly(trio_dir, output_dir)



    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    
    # Start the process monitoring in the background
    monitor_process = multiprocessing.Process(target=monitor_processes)
    monitor_process.start()
    
    for data in range(100):
        while not check_resources():
            pass  # Wait until resources are available
        
        # Start a new process to run program A
        result = pool.apply_async(run_program_A, args=(data,))
        results.append(result)
    
    # Wait for all processes to complete
    for result in results:
        result.get()
    
    # Wait for the process monitoring to finish
    monitor_process.join()

    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-i', '--input_dir', help = 'input directory path', required = True)
    parser.add_argument('-g', '--ref_genome', help = 'reference genome file', required = True)
    parser.add_argument('-t', '--n_thread', help = 'number of threads to run etching (default = 20)', required = False, default = 20)
    parser.add_argument('-k', '--kmer_db', help = 'KMC kmer db ', required = True)
    parser.add_argument('-o', '--output', help = 'output prefix', required = False, default = 'ETCHING_run.{}'.format(date_time))
    #parser.add_argument('-q', '--apply_bqsr', type = bool, help = 'sickle option True/False', required = False, default = True)
    args = parser.parse_args()
    #run_etching(args)
    #__main__(args)
