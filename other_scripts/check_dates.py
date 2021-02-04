import glob
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

files15 = glob.glob("TIMPS_2015/*.nc")
datesT15 = [f.split("_")[3][4:8] for f in files15]
count15 = Counter(datesT15)
print(np.mean(list(count15.values())))

files16 = glob.glob("TIMPS_2016/*.nc")
datesT16 = [f.split("_")[3][4:8] for f in files16]
count16 = Counter(datesT16)
print(np.mean(list(count16.values())))

files17 = glob.glob("TIMPS_2017/*.nc")
datesT17 = [f.split("_")[3][4:8] for f in files17]
count17 = Counter(datesT17)
print(np.mean(list(count17.values())))

files18 = glob.glob("TIMPS_2018/*.nc")
datesT18 = [f.split("_")[3][4:8] for f in files18]
count18 = Counter(datesT18)
print(np.mean(list(count18.values())))

files19 = glob.glob("TIMPS_2019/*.nc")
datesT19 = [f.split("_")[3][4:8] for f in files19]
count19 = Counter(datesT19)
print(np.mean(list(count19.values())))

#print("Plotting")
#plt.bar(count17.keys(),count17.values())
#plt.show()
