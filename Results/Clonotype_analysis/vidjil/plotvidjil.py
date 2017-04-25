
from collections import Counter
from collections import defaultdict
from bisect import bisect_left

def readFastaMul(nomFi):
	listofall=[]
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	for l in lines:
		listofall.append(int(l.split("\t")[4][:-1]))
	return listofall





def count_intervals(sequence, intervals):
    count = defaultdict(int)
    intervals.sort()
    for item in sequence:
        pos = bisect_left(intervals, item)
        if pos == len(intervals):
            count[None] += 1
        else:
            count[intervals[pos]] += 1
    return count


listequal=readFastaMul("outvidjil.txt")
data=[1,1,1,1,1,10,199998]
#print Counter(listequal)
print count_intervals(listequal, [10,100,1000,10000,100000,80229])
#print count_intervals(data, [10])
