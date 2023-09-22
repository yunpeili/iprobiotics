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

# Test 
# print(nucleic_acid_count('AWCBN'))
# print(nucleic_acid_expand('AWCBN'))

#创建kmer序列与对应tf的字典
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
            expand_result = nucleic_acid_expand(dna_seq=key)
            for key in expand_result:
                curr_dict[key] += 1/len(expand_result)
        
        # curr_dict = {k:v for k,v in curr_dict.items()}#字典储存在curr_dict 

    for key in curr_dict.keys():                      
        l = len(seq) - len(key) + 1
        curr_dict[key] = curr_dict[key] / l           #计算tf
    return(curr_dict)

def feat2str(int_float_dict):
    feat_str = ''
    for index, value in int_float_dict.items():
        feat_str += str(index) + ':' + str(value) + ' '  
    #feat_str = str(label) + '\t' + feat_str 
    return feat_str

folder_path_0 = '/home/nusri/下载/dna_feat/dataset/single_seq_nonprobiotic'
label_0 = 0
dirs_0 = [path for path in os.listdir(folder_path_0) if path.endswith('txt')]
seq_list_0 = []
filename_list_0 =[]
for input_file in dirs_0:
# for input_file in dirs:
    with open(os.path.join(folder_path_0, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list_0.append(dna_seq)
        filename_list_0.append(input_file)

folder_path_1 = '/home/nusri/下载/dna_feat/dataset/single_seq_probiotic'
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

def get_feat_str_list(seq_list, filename_list):
    feat_str_list = []
    for seq, filename in tqdm(zip(seq_list, filename_list), total=len(seq_list)):
        feat_str = ''
        try:
            feat_dict = get_dict(seq)
        except Exception as e:
            print('ERROR: ', filename)
            raise('ERROR')
        #给kmer序列标号并创建字典
        index_ATCG_map = {ATCG: int(i) for i, ATCG in enumerate(feat_dict.keys())}
        #创建index和tf的字典
        feat_dict = {index_ATCG_map[ATCG]: value for i, (ATCG, value) in enumerate(feat_dict.items())}
        feat_str += feat2str(feat_dict)
        feat_str_list.append(feat_str)
    return feat_str_list

with open('result_nonprobio.tsv', 'w') as f:
    for feat_str in get_feat_str_list(seq_list_0, filename_list_0):
        f.write(str(label_0) + '\t' + feat_str + '\n')
    for feat_str in get_feat_str_list(seq_list_1, filename_list_1):
        f.write(str(label_1) + '\t' + feat_str + '\n')
