import sys
from collections import *

def main():
  counter = Counter()
  common_words = defaultdict(float)
  norm = 0.
  for line in open(sys.argv[1]):
    line = line.strip()
    line.replace(" ","")
    if line == '':
      continue
    words = line.split('_')
    for word in words:
      word = word.strip()
      if word  == "":
        continue
      counter[word] +=1 
  for word,count in counter.most_common(5000):
    #print word
    norm += count
    #print counter[word]
    common_words[word] = counter[word]
  for word in common_words:
    common_words[word] /= norm;
    print "%s\t%s"%(word,repr(common_words[word]))
if __name__=="__main__":
  main()
