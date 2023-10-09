import os
from tqdm import tqdm
from collections import Counter

#生成所有k-mer序列
dna_list = ['A', 'G', 'T', 'C']

def all_possible_seq(k):
    
    def generate_seq(k, chars, current_string, results):
        if k == 0:
            results.append(current_string)
            return
        for char in chars:
            generate_seq(k - 1, chars, current_string + char, results)

    results = []
    generate_seq(k, dna_list, '', results)
    return results


ambiguity_symbol_dict = {
    'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
    'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 
    'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'], 
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}
# e.g. 'W':['A', 'T'] --> 'A':0.5, 'T':0.5
#      'D':['A', 'G', 'T'] --> 'A':1/3, 'G':1/3, 'T':1/3
#      'A':['A'] --> 'A':1/1

#print(ambiguity_symbol_dict)

def nucleic_acid_expand(dna_seq):
    results = ['']
    for base in dna_seq:
        results = [ seq + new_seq 
                  for seq in results 
                    for new_seq in ambiguity_symbol_dict[base]]

    return results


import itertools
# from line_profiler import profile

expand_dict = {}

def nucleic_acid_expand_product(dna_seq):
    if dna_seq in expand_dict:
        return expand_dict[dna_seq]
    else:
        bases = [ambiguity_symbol_dict[base] for base in dna_seq]
        combinations = [''.join(comb) for comb in itertools.product(*bases)]
        expand_dict[dna_seq] = combinations
    return combinations


# Test 
# print(nucleic_acid_count('AWCBN'))
# print(nucleic_acid_expand('AWCBN'))

#创建kmer序列与对应tf的字典
# @profile
def get_dict(seq):
    curr_dict = {}
    for k in range(2, 5):
        dna_segs = all_possible_seq(k)                #k-mer序列
        new_dict = {seg:0 for seg in dna_segs}        #创建一个初始字典，key为k-mer序列，value为0
        curr_dict.update(new_dict)                         #在k值更新之后更新字典
        all_occur = len(seq) - k + 1                  #k-mer序列的数量
        for x in range(all_occur):                    #计算每个k-mer序列的出现次数
            key = seq[x: x+k]                         
            # curr_dict[key] += 1
            # curr_dict.update({ key: curr_dict[key]+1 })
            # curr_dict.update({ key1: curr_dict[key]+0.5, key2: curr_dict[key]+0.5 })
            # curr_dict.update({ lambda N: keyN: curr_dict[key]+1/N })
            expand_result = nucleic_acid_expand_product(dna_seq=key)
            for key in expand_result:
                curr_dict[key] += 1/len(expand_result)
        
        # curr_dict = {k:v for k,v in curr_dict.items()}#字典储存在curr_dict 

    for key in curr_dict.keys():                      
        l = len(seq) - len(key) + 1
        curr_dict[key] = curr_dict[key] / l           #计算tf
    return curr_dict

def feat2str(int_float_dict):
    feat_str = ''
    for index, value in int_float_dict.items():
        feat_str += str(index) + ':' + str(value) + ' '  
    #feat_str = str(label) + '\t' + feat_str 
    return feat_str

folder_path_0 = './dataset/single_seq_nonprobiotic'
label_0 = 0
dirs_0 = [path for path in os.listdir(folder_path_0) if path.endswith('txt')]
seq_list_0 = []
filename_list_0 =[]
for input_file in dirs_0[:10]:
# for input_file in dirs:
    with open(os.path.join(folder_path_0, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list_0.append(dna_seq)
        filename_list_0.append(input_file)

folder_path_1 = './dataset/single_seq_probiotic'
label_1 = 1
dirs_1 = [path for path in os.listdir(folder_path_1) if path.endswith('txt')]
seq_list_1 = []
filename_list_1 =[]
for input_file in dirs_1:
# for input_file in dirs:
    with open(os.path.join(folder_path_1, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list_1.append(dna_seq)
        filename_list_1.append(input_file)


import multiprocessing as mp
from multiprocessing import Pool, Queue, Manager, Process
# num_process = mp.cpu_count()
num_process = 2
print(f'using {num_process} cores')


def get_feat_str_list(seq_list, feat_str_list, queue):
    for seq in seq_list:
        feat_dict = get_dict(seq)
    # #给kmer序列标号并创建字典
    # index_ATCG_map = {ATCG: int(i) for i, ATCG in enumerate(feat_dict.keys())}
    # #创建index和tf的字典
    # feat_dict = {index_ATCG_map[ATCG]: value for i, (ATCG, value) in enumerate(feat_dict.items())}
        feat_dict = {i: value for i, (ATCG, value) in enumerate(feat_dict.items())}
        feat_str = feat2str(feat_dict)
        feat_str_list.append(feat_str)
        queue.put(1)

def tqdm_listener(total:int, queue:Queue):
    pbar = tqdm(total=total)
    while True:
        if not queue.empty():
            k = queue.get()
            if k == 1:
                pbar.update(1)
            else:
                break
    pbar.close()

def process_feat(seq_list, filename_list):

    manager = Manager()
    # all_feat_str_list = manager.list()
    feat_str_list = manager.list()
    tqdm_queue = manager.Queue()
    # iter_count = manager.Value('i', 0)
    len(seq_list)

    batch_size = len(seq_list) // num_process
    batches = [seq_list[i:i+batch_size] for i in range(0, len(seq_list), batch_size)]

    tqdm_process = Process(target=tqdm_listener, args=(len(seq_list), tqdm_queue))
    tqdm_process.start()

    # start all process
    processes = []
    for batch in batches:
        p = Process(target=get_feat_str_list, args=(batch, feat_str_list, tqdm_queue))
        processes.append(p)
        p.start()

    # wait all finish
    for p in processes:
        p.join()

    tqdm_queue.put(-1)
    tqdm_process.join()

    # with Pool(processes=num_process) as pool:
    #     results = []
    #     for dna_seq in seq_list:
    #         result = pool.apply_async(get_feat_str, args=(dna_seq, shared_expand_dict))
    #         results.append(result)

    #     with tqdm(total=len(seq_list)) as pbar:
    #         for result in results:
    #             feat_str_list.append(result.get())
    #             pbar.update(1)

            # for result in pool.imap(process_seq, seq_list):
            
    # pool.close()
    # pool.join()

    return list(feat_str_list)


nonprobio_feat_list = process_feat(seq_list_0, filename_list_0)
probio_feat_list = process_feat(seq_list_1, filename_list_1)

with open('result_nonprobio.tsv', 'w') as f:
    for feat_str in nonprobio_feat_list:
        f.write(str(label_0) + '\t' + feat_str + '\n')


with open('result_probio.tsv', 'w') as f:
    for feat_str in probio_feat_list:
        f.write(str(label_1) + '\t' + feat_str + '\n')
