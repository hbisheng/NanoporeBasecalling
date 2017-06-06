from os import listdir
from os.path import isfile, join
"""
files = ["nanocall7.3.error_profile.txt", \
        "nanocall9.error_profile.txt", \
        "deepnano_metrichor_r73.error_profile.txt", \
        "deepnano_r9.error_profile.txt", \
        "metrichor_r9.error_profile.txt"]
"""
files = ["nanocall7.3.new_error_profile.txt", \
        "nanocall9.new_error_profile.txt", \
        "deepnano_metrichor_r73.new_error_profile.txt", \
        "deepnano_r9.new_error_profile.txt", \
        "metrichor_r9.new_error_profile.txt"]
""" 
error_profile file format: 
    [0]QUERY_NAME
    [1]READ_LENGTH
    [2]ALIGNMENT_LENGTH   
    [3]HARD    
    [4]SOFT    
    [5]EQUAL   
    [6]DIFF    
    [7]INS     
    [8]DEL
"""

nanocall_r7 = {}
nanocall_r9 = {}
deepnano_r7_1d = {}
deepnano_r7_2d = {}
deepnano_r9 = {}
metrichor_r7_1d = {}
metrichor_r7_2d = {}
metrichor_r9 = {}
id2filename = {}

def add_to_dict(query_name, items):
    if query_name.count(":") > 0:
        # Nanocall output
        file_name = query_name.split(":")[1]
        if file_name.startswith("LomanLabz"):
            # R7.3
            seq = query_name.split(":")[2]
            if int(seq) == 0:
                file_name = file_name + "_template"
            else:
                file_name = file_name + "_complement"
            nanocall_r7[file_name] = items
        else:
            # R9
            id = query_name.split(":")[0]
            id2filename[id] = file_name
            nanocall_r9[file_name] = items
    elif query_name.count("/") > 0:
        # Deepnano R9
        filename = query_name.split("/")[-1].split(".")[0]
        deepnano_r9[filename] = items
    elif query_name.count(".") > 0:
        # Deepnano R7 or Metrichor R7
        filename = query_name.split(".")[0]
        suffix = query_name.split(".")[-1]
        if suffix.endswith("2d_rnn"):
            # Deepnano R7 2D
            deepnano_r7_2d[filename] = items
        elif suffix.endswith("template_rnn"):
            # Deepnano R7 1D template
            deepnano_r7_1d[filename + "_template"] = items
        elif suffix.endswith("complement_rnn"):
            # Deepnano R7 1D complement
            deepnano_r7_1d[filename + "_complement"] = items
        elif suffix.endswith("2d"):
            metrichor_r7_2d[filename] = items
        elif suffix.endswith("template"):
            metrichor_r7_1d[filename + "_template"] = items
        elif suffix.endswith("complement"):
            metrichor_r7_1d[filename + "_complement"] = items
    else:
        # Metrichor R9
        id = query_name.split("_")[0]
        if id in id2filename:
            metrichor_r9[id2filename[id]] = items


"""
    [0]READ_LENGTH [1]ALIGNMENT_LENGTH [2]HARD [3]SOFT [4]EQUAL [5]DIFF [6]INS [7]DEL
"""

def cal_identity(items):
    if items[1] == 0:
        return 0
    return 1.0 * items[4] / (items[4] + items[5] + items[6] + items[7])

def cal_fraction_aligned(items):
    if items[1] == 0:
        return 0
    return 1.0 * items[0] / (items[0] + items[2] + items[3])

def get_avg_statistics(dic):
    ident_sum = 0.0
    fraction_aligned_sum = 0.0
    read_length_sum = 0
    align_length_sum = 0
    insertions_sum = 0
    deletions_sum = 0
    aligned_num = 0
    for key in dic:
        items = dic[key]
        if items[1] > 0: # ignore those unmapped to the reference
            aligned_num += 1
            ident_sum += 1.0 * items[4] / (items[4] + items[5] + items[6] + items[7])
            fraction_aligned_sum += 1.0 * items[0] / (items[0] + items[2] + items[3])
            read_length_sum += items[0]
            align_length_sum += items[1]
            insertions_sum += 1.0 * items[6] / (items[4] + items[5] + items[6] + items[7])
            deletions_sum += 1.0 * items[7] / (items[4] + items[5] + items[6] + items[7])
    print "\tReads number", len(dic)
    print "\tAligned percentage", 1.0 * aligned_num / len(dic)
    print "\tAvg read length", 1.0 * read_length_sum / aligned_num
    print "\tAvg alignment length", 1.0 * align_length_sum / aligned_num
    print "\tAvg identity", ident_sum / aligned_num
    print "\tAvg fraction aligned", fraction_aligned_sum / aligned_num
    print "\tAvg insertion percentage", 1.0 * insertions_sum / aligned_num
    print "\tAvg deletion percentage", 1.0 * deletions_sum / aligned_num

for file in files:
    with open(file, 'r') as fi:
        for line in fi:
            items = line[:-1].split("\t")
            add_to_dict(items[0], [float(val) for val in items[1:]])

print "NANOCALL_R7"
get_avg_statistics(nanocall_r7)
print "NANOCALL_R9"
get_avg_statistics(nanocall_r9)
print "DEEPNANO_R7_1D"
get_avg_statistics(deepnano_r7_1d)
print "DEEPNANO_R7_2D"
get_avg_statistics(deepnano_r7_2d)
print "DEEPNANO_R9"
get_avg_statistics(deepnano_r9)
print "METRICHOR_R7_1D"
get_avg_statistics(metrichor_r7_1d)
print "METRICHOR_R7_2D"
get_avg_statistics(metrichor_r7_2d)
print "METRICHOR_R9"
get_avg_statistics(metrichor_r9)


f = open("r7_1d_res_new", "w")
for key in deepnano_r7_1d:
    if key not in nanocall_r7:
        f.write("0, 0, 0, 0")
    else:
        items = nanocall_r7[key]
        f.write("%d, %d, %f, %f" % (items[0],\
            items[1], \
            cal_fraction_aligned(items), 
            cal_identity(items)))
    items = deepnano_r7_1d[key]
    f.write(",%d, %d, %f, %f" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))

    items = metrichor_r7_1d[key]
    f.write(",%d, %d, %f, %f\n" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))
f.close()

f = open("r7_2d_res_new", "w")
for key in deepnano_r7_2d:
    items = deepnano_r7_2d[key]
    f.write("%d, %d, %f, %f" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))

    items = metrichor_r7_2d[key]
    f.write(",%d, %d, %f, %f\n" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))
f.close()

f = open("r9_res_new", "w")
for key in deepnano_r9:
    
    items = nanocall_r9[key]
    f.write("%d, %d, %f, %f" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))
    
    items = deepnano_r9[key]
    f.write(",%d, %d, %f, %f" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))

    items = metrichor_r9[key]
    f.write(",%d, %d, %f, %f\n" % (items[0],\
        items[1], \
        cal_fraction_aligned(items), 
        cal_identity(items)))
f.close()
